#include "PyIMAQ.h"
#include "pch.h"
#include <windows.h>
#include <mmsystem.h>
#include <stdio.h>
#define _NIWIN  
#include "niimaq.h"
#include "fftw3.h"
#include "time.h"
#include <algorithm>
#include <complex>


// IMAQ handles
static SESSION_ID session_id;
static BUFLIST_ID buflist_id;
static INTERFACE_ID interface_id;
static Int8** imaq_buffers;

// Dimensions of acquisition
static Int32 acqWinWidth;
static Int32 acqWinHeight;
static Int32 bytesPerPixel;

// Size of each buffer/frame
static unsigned int buffer_size;
static unsigned int number_of_bscans;
static unsigned int bwidth;  // Equal to acqWinHeight / number_of_bscans

// Master buffer number since acquisiton start
static uInt32 buffer_number = NULL;

// User-supplied buffer count
static int number_of_buffers;


// Flag that turns linear interpolation of lambda spectrum on and off
static bool interp_enabled;

// Flag that turns apod on and off
static bool apod_enabled;
static float* apod_window; // Window like background spectrum or 

// 32 bit fft plan for FAST-SDOCT
static fftwf_plan fft_plan;

static uint16_t* examined_addr = NULL; // Pointer to the examined IMAQ buffer
static uInt32* examined_count = NULL;  // The buffer returned by ExamineBuffer

static float* src_frame;  // Input array for FFT

// Wavelength to wavenumber interpolation params
static float* k;  // Linear-in-wavenumber array
static float* lam; // Linear-in-wavelength array
static float d_lam; // Wavelength spacing
static int* nn;  // LUT of nearest-neighbor indices for each wavenumber in k

// Motion parameters and memory
static fftwf_complex *t0;  // Reference frame
static fftwf_complex *tn;  // Delta frame
static fftwf_complex *r; // Phase correlation output

static std::complex<float> i_t0;
static std::complex<float> i_tn;
static std::complex<float> pc;
static std::complex<float> pcnorm;

static int* px_dx_lut;  // LUT to convert from FFT pixels to pixel displacements
static unsigned int pc_roi;  // The offset into the top of the A-scan from which to crop the bwidth x bwidth ROI
static unsigned int bscan_idx;  // Index of current B-scan being processed

float* fourier_filter;  // Filter applied to the 2D spatial spectral images prior to phase correlation
float* apod_filter; // Pixmap of attenuation values applied to 2D spatial images prior to forward FFT 

// 32 bit 2D fft plans for phase correlation
static fftwf_plan pc_t0_fft_plan;
static fftwf_plan pc_tn_fft_plan;
static fftwf_plan pc_r_fft_plan;

std::complex<float> fftwf_conj(std::complex<float> imag)
{
	return std::conj(imag);
}

extern "C"
{

	__declspec(dllexport) int interfaceReset()
	{

		return imgInterfaceReset(interface_id);

	}

	__declspec(dllexport) int rebootCamera()
	{

		char* cmdbf = "REBOOT\r";
		uInt32 size = sizeof(cmdbf);

		int err = 0;
		err = imgSessionSerialFlush(session_id);
		err = imgSessionSerialWrite(session_id, cmdbf, &size, 3000);

		return err;

	}


	__declspec(dllexport) int open(char* name)
	{
		int error = 0;

		error = imgInterfaceOpen(name, &interface_id);
		error = imgSessionOpen(interface_id, &session_id);

		return error;
	}


	__declspec(dllexport) int close()
	{
		int error = 0;

		printf("Unlocking buffer memory...\n");
		error = imgMemUnlock(buflist_id);

		printf("Disposing of buffers...\n");
		for (int i = 0; i < number_of_buffers; i++)
		{
			if (imaq_buffers[i] != NULL)
			{
				error = imgDisposeBuffer(imaq_buffers[i]);
			}
		}

		printf("Disposing of buffer list...\n");
		error = imgDisposeBufList(buflist_id, FALSE);
		printf("Closing session...\n");
		error = imgClose(session_id, TRUE);
		printf("Closing interface...\n");
		error = imgClose(interface_id, TRUE);

		return error;
	}


	__declspec(dllexport) int defineBufferSize()
	{
		int error = 0;

		error = imgGetAttribute(session_id, IMG_ATTR_ROI_WIDTH, &acqWinWidth);
		error = imgGetAttribute(session_id, IMG_ATTR_ROI_HEIGHT, &acqWinHeight);
		error = imgGetAttribute(session_id, IMG_ATTR_BYTESPERPIXEL, &bytesPerPixel);
		printf("defineBufferSize: width %i, height %i, bytes %i\n", acqWinWidth, acqWinHeight, bytesPerPixel);

		buffer_size = acqWinWidth * acqWinHeight * bytesPerPixel;

		return error;
	}


	__declspec(dllexport) int getBufferSize()
	{

		int error = 0;

		error = imgGetAttribute(session_id, IMG_ATTR_ROI_WIDTH, &acqWinWidth);
		error = imgGetAttribute(session_id, IMG_ATTR_ROI_HEIGHT, &acqWinHeight);
		error = imgGetAttribute(session_id, IMG_ATTR_BYTESPERPIXEL, &bytesPerPixel);

		buffer_size = acqWinWidth * acqWinHeight * bytesPerPixel;

		return buffer_size;

	}


	__declspec(dllexport) int getFrameSize()
	{

		int error = 0;

		error = imgGetAttribute(session_id, IMG_ATTR_ROI_WIDTH, &acqWinWidth);
		error = imgGetAttribute(session_id, IMG_ATTR_ROI_HEIGHT, &acqWinHeight);

		return acqWinWidth * acqWinHeight;

	}


	__declspec(dllexport) int sessionSerialWrite(char* msg)
	{
		uInt32 msg_length = strlen(msg);
		imgSessionSerialFlush(session_id);
		return imgSessionSerialWrite(session_id, msg, &msg_length, 3000);
	}

	__declspec(dllexport) int sessionSerialRead(char* buff, uInt32* length)
	{
		return imgSessionSerialRead(session_id, buff, length, 3000);
	}


	__declspec(dllexport) int getBufferDim(int* x, int* y, int* bytes)
	{
		int error = 0;

		error = imgGetAttribute(session_id, IMG_ATTR_ROI_WIDTH, x);
		error = imgGetAttribute(session_id, IMG_ATTR_ROI_HEIGHT, y);
		error = imgGetAttribute(session_id, IMG_ATTR_BYTESPERPIXEL, bytes);

		return error;

	}


	__declspec(dllexport) int initBuffers(int buffer_no)
	{
		unsigned int bufCmd;

		number_of_buffers = buffer_no;

		int error = 0;

		error = defineBufferSize();

		error = imgCreateBufList(number_of_buffers, &buflist_id);

		imaq_buffers = new Int8 * [number_of_buffers];

		for (int i = 0; i < number_of_buffers; i++)
		{
			error = imgCreateBuffer(session_id, FALSE, buffer_size, (void**)&imaq_buffers[i]);
			error = imgSetBufferElement2(buflist_id, i, IMG_BUFF_ADDRESS, imaq_buffers[i]);
			error = imgSetBufferElement2(buflist_id, i, IMG_BUFF_SIZE, buffer_size);
			bufCmd = (i == (number_of_buffers - 1)) ? IMG_CMD_LOOP : IMG_CMD_NEXT;
			error = imgSetBufferElement2(buflist_id, i, IMG_BUFF_COMMAND, bufCmd);
		}

		error = imgMemLock(buflist_id);

		error = imgSessionConfigure(session_id, buflist_id);

		return error;

	}


	__declspec(dllexport) int startAcq()
	{
		return imgSessionAcquire(session_id, TRUE, NULL);
	}


	__declspec(dllexport) int stopAcq()
	{
		return imgSessionStopAcquisition(session_id);
	}


	__declspec(dllexport) int abortAcq()
	{
		return imgSessionAbort(session_id, &buffer_number);
	}


	__declspec(dllexport) int sessionConfigure()
	{
		return imgSessionConfigure(session_id, buflist_id);
	}


	__declspec(dllexport) int getFrame(int frame_no, UINT16* frame_dst, int* buffer_number_out)
	{
		int error = 0;

		void* buffer_addr = NULL;

		error = imgSessionExamineBuffer2(session_id, frame_no, &buffer_number, &buffer_addr);

		memcpy(frame_dst, buffer_addr, buffer_size);  // Copy buffer to output

		error = imgSessionReleaseBuffer(session_id);  // Release the buffer

		return error;

	}


	__declspec(dllexport) int getDroppedFrames(int* dst)
	{

		return imgGetAttribute(session_id, IMG_ATTR_LOST_FRAMES, dst);

	}


	__declspec(dllexport) int getCurrentFrame(UINT16* frame_dst, int* buffer_number_out)
	{
		int error = 0;

		void* buffer_addr = NULL;

		error = imgSessionExamineBuffer2(session_id, IMG_CURRENT_BUFFER, &buffer_number, &buffer_addr);

		memcpy(frame_dst, buffer_addr, buffer_size);  // Copy buffer to output

		error = imgSessionReleaseBuffer(session_id);  // Release the buffer

		return error;

	}


	__declspec(dllexport) int copyBuffer(int frame_no, UINT16* frame_dst)
	{

		return imgSessionCopyBuffer(session_id, frame_no, (uInt8*)frame_dst, 0);

	}

	// -- FAST SD-OCT FUNCTIONS -------------------------------------------------------------------


	__declspec(dllexport) int motion_plan(int zstart, int bscans, float* apodization_filter, float* fourier_domain_filter)
	{

		examined_count = new uInt32[1];

		fourier_filter = fourier_domain_filter;
		apod_filter = apodization_filter;

		number_of_bscans = bscans;
		bwidth = (int)(acqWinHeight / number_of_bscans);

		printf("Planning phase correlation motion tracking for %i B-scans of width %i\n", number_of_bscans, bwidth);

		pc_roi = zstart;

		t0 = fftwf_alloc_complex((int)((acqWinWidth * acqWinHeight) / 2));  // A subsection of the array bwidth x bwidth will actually be used for the transform
		tn = fftwf_alloc_complex((int)((acqWinWidth * acqWinHeight) / 2));
		r = fftwf_alloc_complex(bwidth * bwidth);

		printf("Planning FFTWF transforms with complex source %i x %i\n", bwidth, bwidth);
		pc_t0_fft_plan = fftwf_plan_dft_2d(bwidth, bwidth, t0, t0, FFTW_FORWARD, FFTW_PATIENT);
		pc_tn_fft_plan = fftwf_plan_dft_2d(bwidth, bwidth, tn, tn, FFTW_FORWARD, FFTW_PATIENT);
		pc_r_fft_plan = fftwf_plan_dft_2d(bwidth, bwidth, r, r, FFTW_BACKWARD, FFTW_PATIENT);

		return 0;
	}



	__declspec(dllexport) int SDOCT_plan(bool interp_flag, bool apod_flag, float* spectrometer_lambda, float* apod)
	{
		
		if (interp_flag)
		{
			lam = spectrometer_lambda;
			interp_enabled = true;
		}
		else
		{
			interp_enabled = false;
		}

		if (apod_flag)
		{
			apod_window = apod;
			apod_enabled = true;
		}
		else
		{
			apod_enabled = false;
		}
		
		src_frame = fftwf_alloc_real(acqWinHeight * acqWinWidth);
		fftwf_complex* dummy_out = fftwf_alloc_complex((int)((acqWinWidth * acqWinHeight) / 2));

		// FFTW "many" plan
		int n[] = { acqWinWidth };
		int idist = acqWinWidth;
		int odist = (int)(acqWinWidth / 2);
		int istride = 1;
		int ostride = 1;
		int* inembed = n;
		int* onembed = &odist;
		printf("Planning FFTWF many r2c transform: %i transforms of size %i...\n", acqWinHeight, acqWinWidth);
		fft_plan = fftwf_plan_many_dft_r2c(1, n, acqWinHeight, src_frame, inembed, istride, idist, dummy_out, onembed, ostride, odist, FFTW_PATIENT);

		fftwf_free(dummy_out);

		if (interp_enabled)
		{
			printf("Planning lambda to k interpolation...\n");
			// Populate linear-in-wavenumber array
			float lam_max = lam[acqWinWidth - 1];
			float lam_min = lam[0];
			d_lam = lam_max - lam_min;

			float d_k = (1 / lam_min - 1 / lam_max) / 2048;

			// TODO cleanup
			k = new float[acqWinWidth];
			nn = new int[2 * acqWinWidth];
			printf("Allocated nn at %p with elements %i\n", nn, sizeof(nn));

			for (int i = acqWinWidth - 1; i > -1; i -= 1)
			{
				k[i] = 1 / ((1 / lam_max) + d_k * i);
			}

			//(Naively, but only once) find nearest upper and lower indices for linear interpolation
			for (int i = 0; i < acqWinWidth; i++)
			{
				// printf("Defining k array at k = %i\n", i);
				float nearestj;  // Holds index of nearest-yet wavelength
				float nearest = acqWinWidth;  // Default to width of entire spectrometer, it's larger than any possible distance
				for (int j = 1; j < acqWinWidth - 1; j++)  // Can exclude first and last lambda-- all k are within these
				{
					float distance = lam[j] - k[i];
					if (std::abs(distance) < std::abs(nearest))  // If new nearest neighbor found
					{
						nearest = distance;
						nearestj = j;
					}
				}
				// printf("Writing to %p\n", nn[i * 1]);
				if (nearest >= 0)  // Distance greater than 0 means closest point is to the right
				{
					nn[i * 1] = nearestj - 1;  // nn holds upper and lower interpolation bounds for all k
					nn[i * 2] = nearestj;
				}
				else
				{
					nn[i * 1] = nearestj;
					nn[i * 2] = nearestj + 1;
				}
			}
		}
	
		return 0;
	}

	#define MOTION_ERROR 2147483647

	__declspec(dllexport) float SDOCT_motion(int dt, float* output)
	{

		int halfwidth = acqWinWidth / 2;

		// Reference frame grab ----------------------------------------------------------------------------------------------------------------------------------
		examined_addr = NULL;
		if ((imgSessionExamineBuffer2(session_id, 1, examined_count, (void**)&examined_addr) == -1) || (examined_addr == NULL) || (*examined_count == -1))
		{
			return MOTION_ERROR;
		}
		// Copy array from IMAQ buffer and convert it to float
		std::copy(examined_addr, examined_addr + acqWinWidth * acqWinHeight, src_frame);

		if (apod_enabled)
		{
			for (int i = 0; i < acqWinHeight; i++)  // For each A-line
			{
				for (int j = 0; j < acqWinWidth; j++)
				{
					src_frame[i * acqWinWidth + j] *= apod_window[j];
				}
			}
		}

		fftwf_execute_dft_r2c(fft_plan, src_frame, t0);
		imgSessionReleaseBuffer(session_id);
		printf("Acquired reference image\n");

		// Delayed frame grab ----------------------------------------------------------------------------------------------------------------------------------
		examined_addr = NULL;
		if ((imgSessionExamineBuffer2(session_id, 2, examined_count, (void**)&examined_addr) == -1) || (examined_addr == NULL) || (*examined_count == -1))
		{
			return MOTION_ERROR;
		}

		std::copy(examined_addr, examined_addr + acqWinWidth * acqWinHeight, src_frame);

		if (apod_enabled)
		{
			for (int i = 0; i < acqWinHeight; i++)  // For each A-line
			{
				for (int j = 0; j < acqWinWidth; j++)
				{
					src_frame[i * acqWinWidth + j] *= apod_window[j];
				}
			}
		}

		fftwf_execute_dft_r2c(fft_plan, src_frame, tn);
		imgSessionReleaseBuffer(session_id);
		printf("Acquired dt image\n");


		// Phase correlation  ----------------------------------------------------------------------------------------------------------------------------------

		for (int bscan_idx = 0; bscan_idx < number_of_bscans; bscan_idx++)
		{

			printf("Phase correlation of B-scan %i, t0 at %p and tn at %p\n", bscan_idx + 1, t0, tn);

			fftwf_execute_dft(pc_t0_fft_plan, t0, t0);
			fftwf_execute_dft(pc_tn_fft_plan, tn, tn);

			printf("FFTs successful...\n");

			for (int i = 0; i < bwidth; i++)
			{
				// printf("B-scan %i / %i, A-scan %i / %i\n", bscan_idx + 1, number_of_bscans, i + 1, bwidth);
				for (int j = 0; j < bwidth; j++)
				{

					// printf("Defining R matrix element [%i][%i], %i...\n", i, j, (i * bwidth + j));

					// Convert to floats
					memcpy(&i_tn, tn + (bscan_idx * halfwidth * bwidth) + (i * halfwidth + pc_roi + j), sizeof(fftwf_complex));
					memcpy(&i_t0, t0 + (bscan_idx * halfwidth * bwidth) + (i * halfwidth + pc_roi + j), sizeof(fftwf_complex));

					// Correlation
					pc = (i_t0 * std::conj(i_tn));
					pcnorm = pc / std::abs(pc);

					// printf("r = %f at [%i, %i]\n", pcnorm, i, j);

					// Copy phase corr result to R array
					memcpy(r + (i * bwidth + j), &pcnorm, sizeof(fftwf_complex));

				}
			}

			// TODO apply filter to R array before FFT

			fftwf_complex* R = fftwf_alloc_complex(bwidth * bwidth);

			// Convert R back to spatial domain
			fftwf_execute_dft(pc_r_fft_plan, r, R);

			// TODO move these declarations out of the loop
			int maxi = -1;
			int maxj = -1;
			float maxval = -1;
			std::complex<float> tmp;
			// Find peak
			for (int i = 0; i < bwidth; i++)
			{
				for (int j = 0; j < bwidth; j++)
				{
					memcpy(&tmp, R[i * bwidth + j], sizeof(fftwf_complex));
					float mag = std::abs(tmp);
					// printf("R = %f at [%i, %i]\n", mag, i, j);
					if (mag > maxval)
					{
						maxval = mag;
						maxi = i;
						maxj = j;
					}
				}
			}

			printf("Maxval = %f at [%i %i]\n", maxval, maxi, maxj);

		}

		return 0;

	}

	__declspec(dllexport) int SDOCT_GetBuffer(int frame_no, uInt32* count, fftwf_complex* dst_frame)
	{
		// Grab IMAQ buffer
		examined_addr = NULL;
		if ((imgSessionExamineBuffer2(session_id, frame_no, count, (void**)&examined_addr) == -1) || (examined_addr == 0))
		{
			return -1;
		}
		if (*count == -1)
		{ 
			return -1;
		}

		if (interp_enabled)
		{
			uint16_t* lambda_aline;
			for (int i = 0; i < acqWinHeight; i++)  // For each A-line
			{
				lambda_aline = examined_addr + i * acqWinWidth;
				for (int j = 0; j < acqWinWidth; j++)  // For each element of each A-line
				{
					// Casts to float
					float y1 = lambda_aline[nn[j * 1]];  // y-values from neighbors in spectrum
					float y2 = lambda_aline[nn[j * 2]];

					// Lambda and K values
					float x1 = lam[nn[j * 1]];  // corresponding initial wavelength
					float x = k[j];  // linear-in-wavenumber interpolation point

					if (y1 == y2)
					{
						src_frame[acqWinWidth * i + j] = y1;
					}
					else
					{
						src_frame[acqWinWidth * i + j] = y1 + (x - x1) / (y2 - y1) * d_lam;
					}

				}

			}
		}
		else
		{
			std::copy(examined_addr, examined_addr + acqWinWidth * acqWinHeight, src_frame);

		}
		if (apod_enabled)
		{

			for (int i = 0; i < acqWinHeight; i++)  // For each A-line
			{
				for (int j = 0; j < acqWinWidth; j++)  // For each element of each A-line
				{
					src_frame[i * acqWinWidth + j] = apod_window[j];
				}
			}
		}

		// FFT

		fftwf_execute_dft_r2c(fft_plan, src_frame, dst_frame);

		// Release the IMAQ buffer back into the ring
		imgSessionReleaseBuffer(session_id);
		return 0;

	}


	__declspec(dllexport) int SDOCT_cleanup()
	{

		fftwf_destroy_plan(fft_plan);
		fftwf_free(src_frame);
		delete k;
		delete nn;
		return 0;

	}


	__declspec(dllexport) int copyCurrentBuffer(UINT16* frame_dst)
	{

		return imgSessionCopyBuffer(session_id, IMG_CURRENT_BUFFER, (uInt8*)frame_dst, 0);

	}


	__declspec(dllexport) int externalLineTrigConfigure(int trigger_line, bool trigger_rising)
	{
		
		UINT32 triggerPolarity;
		if (trigger_rising)
		{
			triggerPolarity = IMG_TRIG_POLAR_ACTIVEH;
		}
		else
		{
			triggerPolarity = IMG_TRIG_POLAR_ACTIVEL;
		}

		UINT32 triggerLine;
		switch (trigger_line)
		{
			case 0:
				triggerLine = IMG_EXT_TRIG0;
			case 1:
				triggerLine = IMG_EXT_TRIG1;
			case 2:
				triggerLine = IMG_EXT_TRIG2;
			case 3:
				triggerLine = IMG_EXT_TRIG3;
			case 4:
				triggerLine = IMG_EXT_RTSI0;
			case 5:
				triggerLine = IMG_EXT_RTSI1;
			case 6:
				triggerLine = IMG_EXT_RTSI2;
			case 7:
				triggerLine = IMG_EXT_RTSI3;
			case 8:
				triggerLine = IMG_EXT_RTSI4;
			case 9:
				triggerLine = IMG_EXT_RTSI5;
			case 10:
				triggerLine = IMG_EXT_RTSI6;

		}

		return imgSessionLineTrigSource2(session_id, IMG_SIGNAL_EXTERNAL, triggerLine, triggerPolarity, 0);
	}


	__declspec(dllexport) int configureROI(UINT32 top, UINT32 left, UINT32 height, UINT32 width)
	{

		return imgSessionConfigureROI(session_id, top, left, height, width);

	}


	__declspec(dllexport) int setAttributeROI(UINT32 top, UINT32 left, UINT32 height, UINT32 width)
	{
		int err = imgSetAttribute2(session_id, IMG_ATTR_ACQWINDOW_TOP, top);
		err = imgSetAttribute2(session_id, IMG_ATTR_ACQWINDOW_LEFT, left);
		err = imgSetAttribute2(session_id, IMG_ATTR_ACQWINDOW_HEIGHT, height);
		err = imgSetAttribute2(session_id, IMG_ATTR_ACQWINDOW_WIDTH, width);
		return err;
	}


	__declspec(dllexport) int setCameraAttributeString(char* attribute, char* value)
	{

		return imgSetCameraAttributeString(session_id, attribute, value);

	}


	__declspec(dllexport) int setCameraAttributeNumeric(char* attribute, double value)
	{

		return imgSetCameraAttributeNumeric(session_id, attribute, value);

	}


	__declspec(dllexport) void showError(int error_code)
	{
		Int8 error_buffer[256];
		imgShowError(error_code, error_buffer);
		printf(error_buffer);
		printf("\n");
	}


}
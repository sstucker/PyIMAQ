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
#include "string.h"

  
// Define everything at global module level

#define MAX_FRAMESTAMP 4096

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
static int bwidth;  // Equal to acqWinHeight / number_of_bscans
static uInt32 halfwidth;  // The size of the spatial A-line-- not actually half, but N / 2 + 1

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
static uInt32* examined_count_ref = NULL;
static uInt32* examined_count_dt = NULL;

static float* src_frame;  // Input array for FFT
static unsigned int* frame_stamp;  // Array of frame numbers
static unsigned int grabbed;  // Number of times grab function has been called since planning

// Wavelength to wavenumber interpolation params
static float* k;  // Linear-in-wavenumber array
static float* lam; // Linear-in-wavelength array
static float d_lam; // Wavelength spacing
static int* nn;  // LUT of nearest-neighbor indices for each wavenumber in k

// Motion parameters and memory
static int dt;
static int rpts;

static float* mot_output;  // Buffer for result of motion quant

static fftwf_complex* t0;  // Reference frame
static fftwf_complex* t0_roi;  // Reference frame ROI for phase corr
static fftwf_complex* tn;  // Delta frame
static fftwf_complex* tn_roi;  // Delta frame ROI for phase corr
static fftwf_complex* r; // Phase correlation output

static std::complex<float> i_t0;
static std::complex<float> i_tn;
static std::complex<float> pc;
static std::complex<float> pcnorm;

static int* px_dx_lut;  // LUT to convert from FFT pixels to pixel displacements
static unsigned int pc_roi;  // The offset into the top of the A-scan from which to crop the bwidth x bwidth ROI
static unsigned int bscan_idx;  // Index of current B-scan being processed
static UINT16* ref_stamp; // Frame stamps of reference A-lines

static float* fourier_filter;  // Filter applied to the 2D spatial spectral images prior to phase correlation
static float* apod_filter; // Pixmap of attenuation values applied to 2D spatial images prior to forward FFT 

// 32 bit 2D fft plans for phase correlation
static fftwf_plan pc_t0_fft_plan;
static fftwf_plan pc_tn_fft_plan;
static fftwf_plan pc_r_fft_plan;

// Maxima finding buffers
static fftwf_complex* R;
static int maxi;
static int maxj;
static int npeak;
static float dx;
static float dx_norm;
static float dy;
static float dy_norm;
static std::complex<float>* maxima;
static float* avg_corr_x;
static float* avg_corr_y;
static float* avg_pos_x;
static float* avg_pos_y;
static float maxval;
static float mag;
static std::complex<float> tmp;
static int* fftshift;

static float* aline_dc;
static unsigned int aline_dc_count;

int mod(int a, int b)
{
	int r = a % b;
	return r < 0 ? r + b : r;
}

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

		printf("Initializing %i buffers of size %i\n", number_of_buffers, buffer_size);

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

	__declspec(dllexport) int getSourceFFT(int index, float* dst)
	{
		// Returns the A-line at index from src_buffer, the input buffer for the most recent OCT processing FFT
		if (index >= acqWinHeight)
		{
			// Index out of bounds
			return -1;
		}
		memcpy(dst, src_frame + index * acqWinWidth, sizeof(float) * acqWinWidth);
		return 0;
	}

	__declspec(dllexport) int getDC(float* dst)
	{
		// Returns the reference spectrum computed via mean
		memcpy(dst, aline_dc, sizeof(float) * acqWinWidth);
		return 0;
	}

	__declspec(dllexport) int SDOCT_get_pc_buffers(fftwf_complex* t0_dst, fftwf_complex* tn_dst)
	{

		for (int i = 0; i < bwidth; i++)
		{
			for (int j = 0; j < bwidth; j++)
			{
				memcpy(tn_dst + bwidth * i + j, tn_roi + bwidth * i + j, sizeof(fftwf_complex));
				memcpy(t0_dst + bwidth * i + j, t0_roi + bwidth * i + j, sizeof(fftwf_complex));
			}
		}

		return 0;
	}

	__declspec(dllexport) int SDOCT_get_R(fftwf_complex* R_dst)
	{

		for (int i = 0; i < bwidth; i++)
		{
			for (int j = 0; j < bwidth; j++)
			{
				memcpy(R_dst + bwidth * i + j, R + bwidth * i + j, sizeof(fftwf_complex));
			}
		}

		return 0;
	}

	__declspec(dllexport) int SDOCT_get_t0(fftwf_complex* t0_dst)
	{

		for (int i = 0; i < acqWinHeight; i++)
		{
			for (int j = 0; j < halfwidth; j++)
			{
				memcpy(t0_dst + halfwidth * i + j, t0 + halfwidth * i + j, sizeof(fftwf_complex));
			}
		}

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
			printf("Apodization enabled...\n");
			apod_window = new float[acqWinWidth];
			memcpy(apod_window, apod, acqWinWidth * sizeof(float));
			apod_enabled = true;
		}
		else
		{
			apod_enabled = false;
		}

		src_frame = fftwf_alloc_real(acqWinHeight * acqWinWidth);
		halfwidth = acqWinWidth / 2 + 1;
		frame_stamp = new unsigned int[acqWinHeight];
		fftwf_complex* dummy_out = fftwf_alloc_complex(halfwidth * acqWinHeight);
		grabbed = 0;

		// Initialize mean aline matrix
		aline_dc = new float[acqWinWidth];
		std::memset(aline_dc, 0, acqWinWidth * sizeof(float));
		aline_dc_count = 0;

		// FFTW "many" plan
		int n[] = { acqWinWidth };
		int idist = acqWinWidth;
		int odist = halfwidth;
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
				float nearestj = -1;  // Holds index of nearest-yet wavelength
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

	__declspec(dllexport) int motion_plan(int zstart, int weighted_avg_npeak, int repeat, int bscans, float* apodization_filter, float* fourier_domain_filter)
	{
		if (repeat < 1)
		{
			return -1; // Repeat must be a positive integer
		}
		else
		{
			rpts = repeat;
		}
		mot_output = new float[3 * bscans];

		examined_count = new uInt32[1];
		examined_count_ref = new uInt32[1];
		examined_count_dt = new uInt32[1];
		ref_stamp = new UINT16[acqWinHeight];

		number_of_bscans = bscans;
		bwidth = acqWinHeight / number_of_bscans;  // The number of A-scans in each independent B-scan to calculate motion along, as well as the height of the ROI

		fourier_filter = new float[bwidth * bwidth];
		apod_filter = new float[bwidth * bwidth];

		// Copy these because the user may lose track of these arrays after calling "plan"
		memcpy(fourier_filter, fourier_domain_filter, sizeof(float) * bwidth * bwidth);
		memcpy(apod_filter, apodization_filter, sizeof(float) * bwidth * bwidth);

		printf("Planning phase correlation motion tracking for %i B-scans of width %i\n", number_of_bscans, bwidth);

		pc_roi = zstart;

		t0 = fftwf_alloc_complex(halfwidth * acqWinHeight);  // A subsection of the array bwidth x bwidth will actually be used for the transform
		tn = fftwf_alloc_complex(halfwidth * acqWinHeight);
		
		t0_roi = fftwf_alloc_complex(bwidth * bwidth);
		tn_roi = fftwf_alloc_complex(bwidth * bwidth);
		
		r = fftwf_alloc_complex(bwidth * bwidth);
		R = fftwf_alloc_complex(bwidth * bwidth);

		printf("Planning FFTWF transforms with complex source %i x %i\n", bwidth, bwidth);
		pc_t0_fft_plan = fftwf_plan_dft_2d(bwidth, bwidth, t0, t0, FFTW_FORWARD, FFTW_PATIENT);
		pc_tn_fft_plan = fftwf_plan_dft_2d(bwidth, bwidth, tn, tn, FFTW_FORWARD, FFTW_PATIENT);
		pc_r_fft_plan = fftwf_plan_dft_2d(bwidth, bwidth, r, R, FFTW_BACKWARD, FFTW_PATIENT);

		// Populate fftshift matrix
		fftshift = new int[bwidth];
		int idx = 0;
		for (int i = 0; i < bwidth / 2; i++)
		{
			fftshift[idx] = i;
			idx++;
		}
		for (int i = 0 - bwidth / 2; i < 0; i++)
		{
			fftshift[idx] = i;
			idx++;
		}

		// Allocate memory for weighted peak finding
		npeak = weighted_avg_npeak;
		maxima = new std::complex<float>[2 * npeak + 1 * 2 * npeak + 1];
		avg_corr_x = new float[2 * npeak + 1];
		avg_corr_y = new float[2 * npeak + 1];
		avg_pos_x = new float[2 * npeak + 1];
		avg_pos_y = new float[2 * npeak + 1];

		return 0;
	}

	__declspec(dllexport) int SDOCT_motion(int n, int* time_lags, float* output)
	{
		// int n: length of dt
		// int* dt: array of time lags to acquire and average together. Length n
		// float* output: pointer to array length number_of_bscans * 3 where resultant maximum correlation values, and dx and dy estimates are copied

		memset(mot_output, 0, sizeof(float) * 3 * number_of_bscans);

		for (int rpt = 0; rpt < rpts; rpt++)
		{
			// Reference frame grab ----------------------------------------------------------------------------------------------------------------------------------
			examined_addr = NULL;
			if ((imgSessionExamineBuffer2(session_id, IMG_LAST_BUFFER, examined_count_ref, (void**)&examined_addr) == -1) || (examined_addr == NULL) || (*examined_count == -1))
			{
				printf("Reference frame grab failed!\n");
				return -1;  // If frame grab fails
			}

			// Get frame stamp (only of first A-line)
			memcpy(ref_stamp, examined_addr, sizeof(UINT16));

			// Copy array from IMAQ buffer and convert it to float
			for (int i = 0; i < acqWinHeight; i++)  // For each A-line
			{
				src_frame[acqWinWidth * i] = 0; // Zero the stamp value so it doesn't affect the image
				for (int j = 1; j < acqWinWidth; j++)  // For each element of each A-line
				{
					src_frame[acqWinWidth * i + j] = examined_addr[acqWinWidth * i + j];
				}
			}

			// Apodization of reference frame
			if (apod_enabled)
			{

				// Compute mean spectra. TODO add option to precompute this
				std::memset(aline_dc, 0, acqWinWidth * sizeof(float));  // Zero the DC A-line array
				for (int i = 0; i < acqWinHeight; i++)  // For each A-line
				{
					for (int j = 0; j < acqWinWidth; j++)  // For each element of each A-line
					{
						aline_dc[j] += src_frame[acqWinWidth * i + j] / acqWinHeight;
					}
				}

				for (int i = 0; i < acqWinHeight; i++)  // For each A-line
				{
					for (int j = 0; j < acqWinWidth; j++)  // For each element of each A-line
					{
						src_frame[acqWinWidth * i + j] = src_frame[acqWinWidth * i + j] - aline_dc[j];  // Subtract mean A-line / DC
						src_frame[acqWinWidth * i + j] = src_frame[acqWinWidth * i + j] * apod_window[j];  // Multiply by apodization window
					}
				}
			}

			fftwf_execute_dft_r2c(fft_plan, src_frame, t0);  // Reference frame FFT
			imgSessionReleaseBuffer(session_id);

			// TODO move reference frame FFT and ROI here

			for (int dt_idx = 0; dt_idx < n; dt_idx++)  // For each frame dt
			{

				dt = time_lags[dt_idx];

				// Delayed frame grab ----------------------------------------------------------------------------------------------------------------------------------
				examined_addr = NULL;
				if ((imgSessionExamineBuffer2(session_id, *examined_count_ref - dt, examined_count_dt, (void**)&examined_addr) == -1) || (examined_addr == NULL) || (*examined_count == -1))
				{
					printf("dt frame grab failed!\n");
					return -1;
				}

				// Ensure A-lines are actually with delay dt with respect to reference frame
				if (mod(ref_stamp[0] - examined_addr[0], MAX_FRAMESTAMP) != acqWinHeight * dt)
				{
					printf("Tried to compare the wrong frames: got %i with dt = %i!\n", *examined_count_dt, mod(ref_stamp[0] - examined_addr[0], MAX_FRAMESTAMP));
					imgSessionReleaseBuffer(session_id);
					return -1;
				}

				// printf("Repeat %i, t0 buffer number = %i, grabbed tn frame with dt = %i, buffer number %i\n", rpt, *examined_count_ref, dt, *examined_count_dt);

				// Copy array from IMAQ buffer and convert it to float
				for (int i = 0; i < acqWinHeight; i++)  // For each A-line
				{
					src_frame[acqWinWidth * i] = 0; // Zero the stamp value so it doesn't affect the image
					for (int j = 1; j < acqWinWidth; j++)  // For each element of each A-line
					{
						src_frame[acqWinWidth * i + j] = examined_addr[acqWinWidth * i + j];
					}
				}

				// Apodization
				if (apod_enabled)
				{

					// Use the DC spectrum calculated from reference frame

					for (int i = 0; i < acqWinHeight; i++)  // For each A-line
					{
						for (int j = 0; j < acqWinWidth; j++)  // For each element of each A-line
						{
							src_frame[acqWinWidth * i + j] = src_frame[acqWinWidth * i + j] - aline_dc[j];  // Subtract mean A-line / DC
							src_frame[acqWinWidth * i + j] = src_frame[acqWinWidth * i + j] * apod_window[j];  // Multiply by apodization window
						}
					}
				}

				fftwf_execute_dft_r2c(fft_plan, src_frame, tn);
				imgSessionReleaseBuffer(session_id);

				// Phase correlation  ----------------------------------------------------------------------------------------------------------------------------------

				for (int bscan_idx = 0; bscan_idx < number_of_bscans; bscan_idx++)  // For each B-scan in the frame
				{

					// Copy processing FFT output to PC buffers and apodize

					for (int i = 0; i < bwidth; i++)  // For each row of ROI
					{
						for (int j = 0; j < bwidth; j++)  // For each col of ROI
						{
							memcpy(t0_roi + i * bwidth + j, t0 + bscan_idx * bwidth * halfwidth + i * halfwidth + pc_roi + j, sizeof(fftwf_complex));
							memcpy(tn_roi + i * bwidth + j, tn + bscan_idx * bwidth * halfwidth + i * halfwidth + pc_roi + j, sizeof(fftwf_complex));
							t0_roi[i * bwidth + j][0] *= apod_filter[i * bwidth + j];
							tn_roi[i * bwidth + j][0] *= apod_filter[i * bwidth + j];
							t0_roi[i * bwidth + j][1] *= apod_filter[i * bwidth + j];
							tn_roi[i * bwidth + j][1] *= apod_filter[i * bwidth + j];
						}
					}

					// Convert phase correlation ROIs to Fourier domain
					fftwf_execute_dft(pc_t0_fft_plan, t0_roi, t0_roi);
					fftwf_execute_dft(pc_tn_fft_plan, tn_roi, tn_roi);

					// Calculate phase corr matrix R
					for (int i = 0; i < bwidth; i++)
					{
						for (int j = 0; j < bwidth; j++)
						{

							// Convert to floats
							memcpy(&i_tn, tn_roi + i * bwidth + j, sizeof(fftwf_complex));
							memcpy(&i_t0, t0_roi + i * bwidth + j, sizeof(fftwf_complex));

							// Correlation
							pc = (i_t0 * std::conj(i_tn));

							if (std::abs(pc) == 0)
							{
								// Zero division case results in 0, not NaN
								pcnorm = 0;
							}
							else
							{
								pcnorm = pc / std::abs(pc);
							}

							// Copy phase corr result to R array
							memcpy(r + (i * bwidth + j), &pcnorm, sizeof(fftwf_complex));

						}
					}

					// TODO apply fourier domain filter. Not sure that this actually improves performance, however

					// Convert R back to spatial domain
					fftwf_execute_dft(pc_r_fft_plan, r, R);

					maxval = -1;
					maxi = -1;
					maxj = -1;

					// Naive search for max corr pixel

					for (int i = 0; i < bwidth; i++)
					{
						for (int j = 0; j < bwidth; j++)
						{
							memcpy(&tmp, R + i * bwidth + j, sizeof(fftwf_complex));
							mag = std::abs(tmp);

							if (mag > maxval)
							{
								maxval = mag;
								maxi = fftshift[i];
								maxj = fftshift[j];
							}
						}
					}

					if (maxval == -1)
					{
						return -1;
					}

					// Weighted average of npeak pixels around max
					if (npeak > 0)
					{
						// Populate matrix of corr values

						int i = 0;
						for (int di = -npeak; di < npeak + 1; di++)
						{
							int j = 0;
							for (int dj = -npeak; dj < npeak + 1; dj++)
							{
								memcpy(maxima + i * 2 * npeak + 1 + j, R + bwidth * mod(maxi + di, bwidth) + mod(maxj + dj, bwidth), sizeof(fftwf_complex));
								j++;
							}
							i++;
						}

						// Take averages along axes

						memset(avg_corr_x, 0, (2 * npeak + 1) * sizeof(float));
						memset(avg_corr_y, 0, (2 * npeak + 1) * sizeof(float));
						memset(avg_pos_x, 0, (2 * npeak + 1) * sizeof(float));
						memset(avg_pos_x, 0, (2 * npeak + 1) * sizeof(float));

						i = 0;
						for (int di = -npeak; di < npeak + 1; di++)
						{
							int j = 0;
							for (int dj = -npeak; dj < npeak + 1; dj++)
							{
								avg_corr_x[i] += std::abs(maxima[2 * npeak + 1 * i + j]);
								avg_pos_x[i] += maxi + di;
								avg_corr_y[j] += std::abs(maxima[2 * npeak + 1 * i + j]);
								avg_pos_y[j] += maxj + dj;
								j++;
							}
							i++;
						}

						// Find weighted dx, dy

						dx = 0;
						dy = 0;
						dx_norm = 0;
						dy_norm = 0;

						for (int i = 0; i < 2 * npeak + 1; i++)
						{
							dx += (avg_corr_x[i] / (2 * npeak + 1)) * (avg_pos_x[i] / (2 * npeak + 1));
							dy += (avg_corr_y[i] / (2 * npeak + 1)) * (avg_pos_y[i] / (2 * npeak + 1));
							dx_norm += (avg_corr_x[i] / (2 * npeak + 1));
							dy_norm += (avg_corr_y[i] / (2 * npeak + 1));
						}

						dx = dx / dx_norm;
						dy = dy / dy_norm;

					}
					else
					{
						dx = maxi;
						dy = maxj;
					}

					// Normalize FFT result
					maxval = maxval / halfwidth;

					// printf("Rpt %i, dt %i, B-scan %i\n", rpt, dt, bscan_idx);
					// printf("Adding to mot_output[%i], output[%i], output[%i]\n", 3 * bscan_idx, 3 * bscan_idx + 1, 3 * bscan_idx + 2);
					mot_output[3 * bscan_idx + 0] += maxval;
					mot_output[3 * bscan_idx + 1] += dx / dt;
					mot_output[3 * bscan_idx + 2] += dy / dt;

				}  // For each B-scan

			}  // For each dt

		}  // For each rpt

		// Copy results to output array
		for (int i = 0; i < (3 * number_of_bscans); i++)
		{
			mot_output[i] = mot_output[i] / (n * rpts);  // Normalize average of all results for n time lags, rpts repeats
			memcpy(output + i, mot_output + i, sizeof(float));
		}

		return 0;

	}

	__declspec(dllexport) int SDOCT_GetBuffer(int frame_no, uInt32* count, fftwf_complex* dst_frame, UINT16* stamp)
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

		// Copy frame stamp array
		for (int i = 0; i < acqWinHeight; i++)
		{
			// printf("A-line %i, %i stamped %i\n", *count, i, examined_addr[i * acqWinWidth]);
			memcpy(stamp + i, examined_addr + acqWinWidth * i, sizeof(UINT16));  // Frame stamp in first index of each line
			*(examined_addr + acqWinWidth * i) = 0;
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

		// Apodization
		if (apod_enabled)
		{

			std::memset(aline_dc, 0, acqWinWidth * sizeof(float));

			// Compute mean spectra
			for (int i = 0; i < acqWinHeight; i++)  // For each A-line
			{
				for (int j = 0; j < acqWinWidth; j++)  // For each element of each A-line
				{
					aline_dc[j] += src_frame[acqWinWidth * i + j] / acqWinHeight;
				}
			}

			for (int i = 0; i < acqWinHeight; i++)  // For each A-line
			{
				for (int j = 0; j < acqWinWidth; j++)  // For each element of each A-line
				{
					src_frame[acqWinWidth * i + j] = src_frame[acqWinWidth * i + j] - aline_dc[j];  // Subtract mean A-line / DC
					src_frame[acqWinWidth * i + j] = src_frame[acqWinWidth * i + j] * apod_window[j];  // Multiply by apodization window
				}
			}
		}


		// FFT

		fftwf_execute_dft_r2c(fft_plan, src_frame, dst_frame);

		// Release the IMAQ buffer back into the ring
		imgSessionReleaseBuffer(session_id);

		return *count;

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

	__declspec(dllexport) int configureTriggerBufferWithTTL1(int timeout)
	{
		return imgSessionTriggerConfigure2(session_id, IMG_SIGNAL_EXTERNAL, IMG_EXT_TRIG1, IMG_TRIG_POLAR_ACTIVEH, timeout, IMG_TRIG_ACTION_BUFFER);
	}

	__declspec(dllexport) int configureROI(UINT32 top, UINT32 left, UINT32 height, UINT32 width)
	{

		int err =  imgSessionConfigureROI(session_id, top, left, height, width);
		
		err = imgSetAttribute2(session_id, IMG_ATTR_ROI_TOP, top);
		err = imgSetAttribute2(session_id, IMG_ATTR_ROI_LEFT, left);
		err = imgSetAttribute2(session_id, IMG_ATTR_ROI_HEIGHT, height);
		err = imgSetAttribute2(session_id, IMG_ATTR_ROI_WIDTH, width);

		err = imgSetAttribute2(session_id, IMG_ATTR_ACQWINDOW_TOP, top);
		err = imgSetAttribute2(session_id, IMG_ATTR_ACQWINDOW_LEFT, left);
		err = imgSetAttribute2(session_id, IMG_ATTR_ACQWINDOW_HEIGHT, height);
		err = imgSetAttribute2(session_id, IMG_ATTR_ACQWINDOW_WIDTH, width);

		return err;

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
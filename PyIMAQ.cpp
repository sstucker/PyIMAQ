#include "PyIMAQ.h"
#include "pch.h"
#include <windows.h>
#include <mmsystem.h>
#include <stdio.h>
#define _NIWIN  
#include "niimaq.h"
#include "fftw3.h"


// 
static SESSION_ID session_id;
static BUFLIST_ID buflist_id;
static INTERFACE_ID interface_id;
static Int8** imaq_buffers;

// Dimensions of acquisition
static Int32 acqWinWidth;
static Int32 acqWinHeight;
static Int32 bytesPerPixel;

// Size of each buffer/frame
unsigned int buffer_size;

// Master buffer number since acquisiton start
static uInt32 buffer_number = NULL;

// User-supplied buffer count
static int number_of_buffers;

// 32 bit fft plan for FAST-SDOCT
static fftwf_plan fft_plan;
static int aline_n;

extern "C"
{
	
	__declspec(dllexport) int interfaceReset()
	{

		return imgInterfaceReset(interface_id);

	}

	__declspec(dllexport) int rebootCamera()
	{

		char* cmdbf = "REBOOT";
		uInt32 size = sizeof(cmdbf);

		int err = 0;
		err = imgSessionSerialFlush(session_id);
		err = imgSessionSerialWrite(session_id, cmdbf, &size, 10);

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

		error = imgMemUnlock(buflist_id);

		for (int i = 0; i < number_of_buffers; i++)
		{
			if (imaq_buffers[i] != NULL)
			{
				error = imgDisposeBuffer(imaq_buffers[i]);
			}
		}

		error = imgDisposeBufList(buflist_id, FALSE);
		error = imgClose(session_id, TRUE);
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
		
		return buffer_size;

	}


	__declspec(dllexport) int getFrameSize()
	{

		return acqWinWidth * acqWinHeight;

	}

	__declspec(dllexport) int sessionSerialWrite(char* msg)
	{
		uInt32 msg_length = strlen(msg);
		printf(msg);
		printf(" size of %i\n", strlen(msg));
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

		imaq_buffers = new Int8*[number_of_buffers];

		for (int i = 0; i < number_of_buffers; i++)
		{
			printf("Configuring buffer %i...\n", i);
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

	__declspec(dllexport) int SDOCT_plan_fft(int aline_size)
	{
		aline_n = aline_size;
		void* dummy = fftwf_alloc_real(aline_size);
		fft_plan = fftwf_plan_dft_r2c_1d(aline_size, (float*)dummy, (fftwf_complex*)dummy, FFTW_PATIENT);
		fftwf_free(dummy);
		return 0;

	}

	__declspec(dllexport) int SDOCT_cleanup()
	{

		fftwf_destroy_plan(fft_plan);
		return 0;

	}

	__declspec(dllexport) int SDOCT_copyBuffer_FFT_ONLY(int frame_no, fftwf_complex* frame_dst)
	{

		imgSessionCopyBuffer(session_id, frame_no, (uInt8*)frame_dst, 0);
		for (int i = 0; i < 1; i++)
		{
			printf("%i", i);
			printf("%i", i * acqWinWidth);
			//fftwf_execute_dft_r2c(fft_plan, (float*)frame_dst[i * acqWinWidth], (fftwf_complex*)frame_dst[i*acqWinWidth]);
		}
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
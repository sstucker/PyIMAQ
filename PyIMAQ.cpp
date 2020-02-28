#include "PyIMAQ.h"
#include "pch.h"
#include <windows.h>
#include <mmsystem.h>
#include <stdio.h>
#define _NIWIN  
#include "niimaq.h"  


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

extern "C"
{
	
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

		error = imgSessionAbort(session_id, NULL);
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

	__declspec(dllexport) int abortAcq()
	{
		return imgSessionAbort(session_id, &buffer_number);
	}

	__declspec(dllexport) int getFrame(int frame_no, UINT16* frame_dst, int* buffer_number_out)
	{
		int error = 0;

		void* buffer_addr = NULL;

		printf("Examining buffer...\n");
		error = imgSessionExamineBuffer2(session_id, frame_no, &buffer_number, &buffer_addr);
		Int8 errbuf[128];
		printf("%i\n",error);
		imgShowError(error, errbuf);
		printf(errbuf); 
		printf("\n");
		fflush(stdout);
		printf("About to memcpy...\n");
		printf("buffer_addr: %i\n", buffer_addr);
		printf("dest addr: %i\n", frame_dst);
		for (int i = 0; i < buffer_size/bytesPerPixel; i+=(buffer_size/bytesPerPixel)/8)
		{
			printf("Buffer index %i / %i\n",i,buffer_size/bytesPerPixel);
			printf("UINT16 value at pixel %i\n", ((UINT16*)buffer_addr)[i]);
			fflush(stdout);
		}
		memcpy(frame_dst, buffer_addr, buffer_size);  // Copy buffer to output

		error = imgSessionReleaseBuffer(session_id);  // Release the buffer

		return error;

	}

	__declspec(dllexport) void showError(int error_code)
	{
		Int8 error_buffer[256];
		imgShowError(error_code, error_buffer);
		printf(error_buffer);
		printf("\n");

	}

}
"""
Convenience functions that allow IMAQ device and ring buffer to be accessed from Python.

No handles to be dealt with; they live as globals as part of the library. Thus only one
interface can be used at a time.

"""

import ctypes as c
import numpy as np
from numpy.ctypeslib import ndpointer
import matplotlib.pyplot as plt

path_to_dll = 'x64/Debug/PyIMAQ.dll'
img = c.CDLL(path_to_dll)

bool_p = ndpointer(dtype=np.bool,ndim=1,flags='C_CONTIGUOUS')
int_p = ndpointer(dtype=np.int,ndim=1,flags='C_CONTIGUOUS')
bool_p = ndpointer(dtype=np.bool,ndim=1,flags='C_CONTIGUOUS')
uint16_p = ndpointer(dtype=np.uint16,ndim=1,flags='C_CONTIGUOUS')
float_p = ndpointer(dtype=np.float32,ndim=1,flags='C_CONTIGUOUS')
double_p = ndpointer(dtype=np.float64,ndim=1,flags='C_CONTIGUOUS')

img.open.argtypes = [c.c_char_p]
img.open.restype = c.c_int
def imgOpen(name):
	name = name.encode('utf-8') 
	return img.open(name)

img.initBuffers.argtypes = [c.c_int]
img.initBuffers.restype = c.c_int
def imgInitBuffer(number_of_buffers):
	return img.initBuffers(number_of_buffers)

img.getBufferSize.restype = c.c_int
def imgGetBufferSize():
	return img.getBufferSize()

img.close.restype = c.c_int
def imgClose():
	return img.close()

img.startAcq.restype = c.c_int
def imgStartAcq():
	return img.startAcq()

img.abortAcq.restype = c.c_int
def imgAbortAcq():
	return img.abortAcq()

img.showError.argtypes = [c.c_int]
def imgShowErrorMsg(error_code):
	img.showError(error_code)

img.getFrame.argtypes = [c.c_int, uint16_p, int_p]
img.getFrame.restype = c.c_int
def imgGetFrame(frame_number, dst, buffnumber=np.empty(1,dtype=np.int)):
	return img.getFrame(frame_number, dst, buffnumber)

# Test

print("imgOpen err ")
imgShowErrorMsg(imgOpen('img0'))

print("imgInitBuffers")
imgShowErrorMsg(imgInitBuffer(3))

print("imgStartAcq")
imgShowErrorMsg(imgStartAcq())

print(imgGetBufferSize(),"bytes/frame")

all = []

for i in range(100):

	fbuff = np.empty(imgGetBufferSize(),dtype=(np.uint16))

	print("getFrame")
	imgShowErrorMsg(imgGetFrame(i,fbuff))
	print("Image range",min(fbuff),max(fbuff))
	all.append(fbuff)

print("imgAbortAcq")
imgShowErrorMsg(imgAbortAcq())

print("imgClose")
imgShowErrorMsg(imgClose())


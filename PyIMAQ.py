"""
Wrapper functions for interface with IMAQ device and a ring buffer

No handles to API objects to be dealt with: are globals as part of the library.
Thus only one interface can be used at a time.
"""
import pathlib
import ctypes as c
import numpy as np
from numpy.ctypeslib import ndpointer

path_to_dll = str(pathlib.Path(__file__).parent.absolute())+'/bin/PyIMAQ.dll'

try:
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
    	return img.open(name.encode('utf-8'))

    img.initBuffers.argtypes = [c.c_int]
    img.initBuffers.restype = c.c_int
    def imgInitBuffer(number_of_buffers):
    	return img.initBuffers(number_of_buffers)

    img.getBufferSize.restype = c.c_int
    def imgGetBufferSize():
    	return img.getBufferSize()

    img.getFrameSize.restype = c.c_int
    def imgGetFrameSize():
    	return img.getFrameSize()

    img.sessionSerialWrite.argtypes = [c.c_char_p, c.c_uint32]
    img.sessionSerialWrite.restype = c.c_int
    def imgSessionSerialWrite(msg):
        msg = msg.encode('utf-8')
        return img.sessionSerialWrite(msg, len(msg))

    img.close.restype = c.c_int
    def imgClose():
    	return img.close()

    img.startAcq.restype = c.c_int
    def imgStartAcq():
    	return img.startAcq()

    img.startAcq.restype = c.c_int
    def imgStopAcq():
    	return img.stopAcq()

    img.abortAcq.restype = c.c_int
    def imgAbortAcq():
    	return img.abortAcq()

    img.showError.argtypes = [c.c_int]
    def imgShowErrorMsg(error_code):
    	img.showError(error_code)

    img.abortAcq.restype = c.c_int
    def imgSessionConfigure():
    	return img.sessionConfigure()

    img.getFrame.argtypes = [c.c_int, uint16_p, int_p]
    img.getFrame.restype = c.c_int
    def imgGetFrame(frame_number, dst, buffnumber=np.empty(1,dtype=np.int)):
    	return img.getFrame(frame_number, dst, buffnumber)

    img.getCurrentFrame.argtypes = [uint16_p, int_p]
    img.getCurrentFrame.restype = c.c_int
    def imgGetCurrentFrame(dst, buffnumber=np.empty(1,dtype=np.int)):
    	return img.getCurrentFrame(dst, buffnumber)

    img.copyBuffer.argtypes = [c.c_int, uint16_p]
    img.copyBuffer.restype = c.c_int
    def imgCopyBuffer(frame_number, dst):
        return img.copyBuffer(frame_number, dst)

    img.copyCurrentBuffer.argtypes = [uint16_p]
    img.copyCurrentBuffer.restype = c.c_int
    def imgCopyCurrentBuffer(dst):
        return img.copyCurrentBuffer(dst)

    img.getDroppedFrames.argtypes = [int_p]
    img.getDroppedFrames.restype = c.c_int
    def imgGetDroppedFrames():
        dropped = np.empty(1,dtype=np.int)
        img.getDroppedFrames(dropped)
        return int(dropped[0])

    img.externalLineTrigConfigure.argtypes = [c.c_int, c.c_bool]
    img.externalLineTrigConfigure.restype = c.c_int
    def imgExternalLineTrigConfigure(trigger_line, trigger_rising):
        return img.externalLineTrigConfigure(trigger_line, trigger_rising)

    img.configureROI.argtypes = [c.c_int, c.c_int, c.c_int, c.c_int]
    img.configureROI.restype = c.c_int
    def imgConfigureROI(top, left, height, width):
        return img.configureROI(top, left, height, width)

    img.setAttributeROI.argtypes = [c.c_int, c.c_int, c.c_int, c.c_int]
    img.setAttributeROI.restype = c.c_int
    def imgSetAttributeROI(top, left, height, width):
        return img.setAttributeROI(top, left, height, width)

    img.setCameraAttributeString.argtypes = [c.c_char_p, c.c_char_p]
    img.setCameraAttributeString.restype = c.c_int
    def imgSetCameraAttributeString(attribute, value):
        attribute = attribute.encode('utf-8')
        value = value.encode('utf-8')
        return img.setCameraAttributeString(attribute, value)

    img.setCameraAttributeNumeric.argtypes = [c.c_char_p, c.c_double]
    img.setCameraAttributeNumeric.restype = c.c_int
    def imgSetCameraAttributeNumeric(attribute, value):
        attribute = attribute.encode('utf-8')
        return img.setCameraAttributeNumeric(attribute, value)

except OSError:
    print('PyIMAQ: PyIMAQ DLL load failed.')
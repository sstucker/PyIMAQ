"""
Wrapper functions for interface with IMAQ device and a ring buffer

No handles to API objects to be dealt with: are globals as part of the library.
Thus only one interface can be used at a time.
"""
import sys
import pathlib
import ctypes as c
import numpy as np
from numpy.ctypeslib import ndpointer

path_to_dll = str(pathlib.Path(__file__).parent.absolute())+'\\bin\\PyIMAQ.dll'

try:
    img = c.CDLL(path_to_dll)

    bool_p = ndpointer(dtype=np.bool,ndim=1,flags='C_CONTIGUOUS')
    int_p = ndpointer(dtype=np.int,ndim=1,flags='C_CONTIGUOUS')
    bool_p = ndpointer(dtype=np.bool,ndim=1,flags='C_CONTIGUOUS')
    uint16_p = ndpointer(dtype=np.uint16,ndim=1,flags='C_CONTIGUOUS')
    uint32_p = ndpointer(dtype=np.uint32,ndim=1,flags='C_CONTIGUOUS')
    float_p = ndpointer(dtype=np.float32,ndim=1,flags='C_CONTIGUOUS')
    double_p = ndpointer(dtype=np.float64,ndim=1,flags='C_CONTIGUOUS')
    complex64_p = ndpointer(dtype=np.complex64,ndim=1,flags='C_CONTIGUOUS')


    img.rebootCamera.restype = c.c_int
    def imgRebootCamera():
        return img.rebootCamera()

    img.getBufferSize.restype = c.c_int
    def imgInterfaceReset():
        return img.interfaceReset()

    img.open.argtypes = [c.c_char_p]
    img.open.restype = c.c_int
    def imgOpen(name):
        try:
            img = c.CDLL(path_to_dll)
        except OSError:
            print('PyIMAQ: failed to open DLL', path_to_dll)
            return -1
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

    img.sessionSerialWrite.argtypes = [c.c_char_p]
    img.sessionSerialWrite.restype = c.c_int
    def imgSessionSerialWrite(msg):
        return img.sessionSerialWrite(c.c_char_p(msg.encode('utf-8')))

    img.sessionSerialRead.argtypes = [c.c_char_p, uint32_p]
    img.sessionSerialRead.restype = c.c_int
    def imgSessionSerialRead(buffer_size=128):
        msgbuff = c.create_string_buffer(buffer_size)
        msglen = np.empty(1, dtype=np.uint32)
        err = img.sessionSerialRead(msgbuff, msglen)
        return str(msgbuff)[0:int(msglen)]

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

    img.configureTriggerBufferWithTTL1.argtypes = [c.c_int]
    img.configureROI.restype = c.c_int
    def imgConfigTrigBufferWithTTL1(timeout):
        return img.configureTriggerBufferWithTTL1(timeout)

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

    # FAST SD-OCT FUNCTIONS -------------------------------------------------------


    img.SDOCT_get_t0.argtypes = [complex64_p]
    img.SDOCT_get_t0.restype = c.c_int
    def octCopyReferenceFrame(dst):
        return img.SDOCT_get_t0(dst)

    img.SDOCT_get_R.argtypes = [complex64_p]
    img.SDOCT_get_R.restype = c.c_int
    def octCopyR(dst):
        return img.SDOCT_get_R(dst)

    img.SDOCT_get_pc_buffers.argtypes = [complex64_p, complex64_p]
    img.SDOCT_get_pc_buffers.restype = c.c_int
    def octCopyPhaseCorrBuffers(t0_dst, tn_dst):
        return img.SDOCT_get_pc_buffers(t0_dst, tn_dst)

    img.getDC.argtypes = [float_p]
    img.getDC.restype = c.c_int
    def octCopyDCSpectrum():
        buffer = np.zeros(2048, dtype=np.float32)
        img.getDC(buffer)
        return buffer

    img.getSourceFFT.argtypes = [c.c_int, float_p]
    img.getSourceFFT.restype = c.c_int
    def octCopySpectrum(index):
        buffer = np.zeros(2048, dtype=np.float32)
        img.getSourceFFT(index, buffer)
        return buffer

    img.motion_plan.argtypes = [c.c_int, c.c_int, c.c_int, c.c_int, float_p, float_p]
    img.motion_plan.restype = c.c_int
    def octMotionPlan(zstart, npeak, number_of_bscans, repeat=1, apod_filter_2d=np.ones(32*32, dtype=np.float32), fourier_filter_2d=np.ones(32*32, dtype=np.float32)):
        return img.motion_plan(zstart, npeak, repeat, number_of_bscans, apod_filter_2d, fourier_filter_2d)
    
    img.SDOCT_motion.argtypes = [c.c_int, int_p, float_p]
    img.SDOCT_motion.restype = c.c_int
    def octMotion(n, dt, output):
        return img.SDOCT_motion(n, dt, output)

    img.SDOCT_plan.argtypes = [c.c_bool, c.c_bool, float_p, float_p]
    img.SDOCT_plan.restype = c.c_int
    def octPlan(lam=None, apod=None):
        if lam is None:
            lambda_arr = np.zeros(1, dtype=np.float32)
            interp_flag = False
        else:
            lambda_arr = lam.astype(np.float32)
            interp_flag = True
        if apod is None:
            apod_arr = np.zeros(1, dtype=np.float32)
            apod_flag = False
        else:
            apod_arr = apod.astype(np.float32)
            apod_flag = True
        return img.SDOCT_plan(interp_flag, apod_flag, lambda_arr, apod_arr)

    def octCleanup():
        return img.SDOCT_cleanup()

    img.SDOCT_GetBuffer.argtypes = [c.c_int, int_p, complex64_p, uint16_p]
    img.SDOCT_GetBuffer.restype = c.c_int
    def octCopyBuffer(frame_number, dst):
        c = np.empty(1, dtype=np.int)
        err = img.SDOCT_GetBuffer(frame_number, c, dst, frame_stamp)
        if err is -1:
            return -1
        else:
            return c

except OSError:
    print('PyIMAQ: failed to open DLL', path_to_dll)
    sys.exit()
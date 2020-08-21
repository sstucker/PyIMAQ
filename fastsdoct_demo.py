# -*- coding: utf-8 -*-

from PyIMAQ import *
import numpy as np
import multiprocessing as mp
from threading import Thread
import matplotlib.pyplot as plt
import matplotlib
import time

# -----------------------------------------------------------------------------

ALINE_SIZE = 2048
ALINES_PER_BUFFER = 10

# -----------------------------------------------------------------------------

camera_name = "img1"
buffer_size = 12

# ----------------------------------------------------------------------------- 

# Open the camera interface
print("imgOpen")
imgShowErrorMsg(imgOpen('img1'))

print('Camera reboot...')
# imgShowErrorMsg(imgRebootCamera())
imgSessionSerialWrite("REBOOT\\r")
time.sleep(2)

octPlanFFT(ALINE_SIZE)

imgSetAttributeROI(0, 0, ALINES_PER_BUFFER, ALINE_SIZE)

# Initialize a ring buffer w/ [buffer_size] buffers
print("imgInitBuffers... ",buffer_size,"buffers in ring...")
imgShowErrorMsg(imgInitBuffer(buffer_size))

# Start camera acquisition
print("imgStartAcq")
imgShowErrorMsg(imgStartAcq())


print(imgGetBufferSize(),"bytes/frame...")
print(imgGetFrameSize(),"A-lines/frame...")

buffer = np.zeros(ALINE_SIZE, dtype=np.uint16)

print('Acquiring buffer')
imgCopyBuffer(0, buffer);
print('Buffer acquired. Mean', np.mean(buffer))

print('Acquiring preprocessed buffer...')
buffer2 = np.zeros(ALINE_SIZE, dtype=np.complex64)
octCopyBufferFFTOnly(0, buffer2)

print("imgStopAcq", flush=True)
imgShowErrorMsg(imgStopAcq())

print("imgClose")
imgShowErrorMsg(imgClose())
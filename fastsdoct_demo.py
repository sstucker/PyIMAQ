# -*- coding: utf-8 -*-

from PyIMAQ import *
import numpy as np
import multiprocessing as mp
from threading import Thread
import matplotlib.pyplot as plt
import matplotlib
import time
import pyfftw
import os

# -----------------------------------------------------------------------------

ALINE_SIZE = 2048
ALINES_PER_BUFFER = 64

# -----------------------------------------------------------------------------

camera_name = "img1"
buffer_size = 10

# -----------------------------------------------------------------------------


# Open the camera interface
print("imgOpen")
imgOpen('img0')

imgSetAttributeROI(0, 0, ALINES_PER_BUFFER, ALINE_SIZE)

# Initialize a ring buffer w/ [buffer_size] buffers
print("imgInitBuffers... ", buffer_size, "buffers in ring...")
imgInitBuffer(buffer_size)

lam = np.linspace(1235, 1385, ALINE_SIZE, dtype=np.float32)
octPlan()
octMotionPlan(20, 2)

# Start camera acquisition
print("imgStartAcq")
imgStartAcq()

print(imgGetBufferSize(), "bytes/frame...")
print(imgGetFrameSize(), "UINT16/frame...")

N_BUFFERS = 10
start = time.time()
for i in range(N_BUFFERS):
    buffer = np.zeros(int(ALINE_SIZE / 2) * ALINES_PER_BUFFER, dtype=np.complex64)
    # print(i, 'returned', octGetBuffer(i, buffer))
    octGetBuffer(i, buffer)
elapsed = time.time() - start
print('Processed', N_BUFFERS * ALINES_PER_BUFFER, 'A-lines in', elapsed)
print('Virtual A-line rate:', (N_BUFFERS * ALINES_PER_BUFFER) / elapsed, 'hz')

result = np.empty(2, dtype=np.float32)
for i in range(10):
    octMotion(1, result)

print("imgStopAcq", flush=True)
imgStopAcq()

octCleanup()

imgClose()

plt.plot(np.abs(buffer))
plt.show()
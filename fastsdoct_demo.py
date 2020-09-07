# -*- coding: utf-8 -*-
from src.main.python.PyScanPattern.ScanPattern import Figure8ScanPattern
from PyIMAQ import *
import numpy as np
import multiprocessing as mp
from threading import Thread
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import matplotlib
import time
import pyfftw
import os
import nidaqmx
from nidaqmx.stream_writers import AnalogMultiChannelWriter
from nidaqmx.constants import LineGrouping, Edge, AcquisitionType

# -----------------------------------------------------------------------------

ALINE_SIZE = 2048
ALINES_PER_BUFFER = 64

# -----------------------------------------------------------------------------

camera_name = "img1"  # img0: free running trigger, img1: manual trigger
buffer_size = 8

# -----------------------------------------------------------------------------

# Create the DAQmx task
scan_task = nidaqmx.Task()
timing = nidaqmx.task.Timing(scan_task)
daq_name = 'Dev1'
daq_channel_ids = [
            daq_name + '/' + 'ao0',
            daq_name + '/' + 'ao1',
            daq_name + '/' + 'ao2'
            ]
for ch_name in daq_channel_ids:
    scan_task.ao_channels.add_ao_voltage_chan(ch_name)

pat = Figure8ScanPattern(debug=True)
pat.generate(0.0625, ALINES_PER_BUFFER / 2, rotation_rad=np.pi/4)
scan_task.timing.cfg_samp_clk_timing(pat.get_sample_rate(),
                                     source="",
                                     active_edge=Edge.RISING,
                                     sample_mode=AcquisitionType.CONTINUOUS)
scan_writer = AnalogMultiChannelWriter(scan_task.out_stream)

scan_signals = pat.get_signals()
scan_signals[0] = scan_signals[0] * 3

scan_writer.write_many_sample(np.array(scan_signals))

# img = plt.imshow(np.zeros([int(ALINE_SIZE / 2), ALINES_PER_BUFFER]))
# plt.show()

# Open the camera interface
print("imgOpen")
imgOpen(camera_name)

imgSetAttributeROI(0, 0, ALINES_PER_BUFFER, ALINE_SIZE)

# Initialize a ring buffer w/ [buffer_size] buffers
print("imgInitBuffers... ", buffer_size, "buffers in ring...")
imgInitBuffer(buffer_size)

lam = np.linspace(1235, 1385, ALINE_SIZE, dtype=np.float32)
octPlan()
octMotionPlan(20, 2)

scan_task.start()

# Start camera acquisition
print("imgStartAcq")
imgStartAcq()

print(imgGetBufferSize(), "bytes/frame...")
print(imgGetFrameSize(), "UINT16/frame...")

acquired = []
N_BUFFERS = 0
start = time.time()

fig = plt.figure()
im = plt.imshow([np.zeros(int(ALINE_SIZE / 2), ALINES_PER_BUFFER), animated=True)

def animate(i):
    buffer = np.zeros(int(ALINE_SIZE / 2) * ALINES_PER_BUFFER, dtype=np.complex64)
    val = octCopyBuffer(i, buffer)
    if val is -1:
        return []
    b = np.transpose(np.array(np.split(buffer, ALINES_PER_BUFFER)))
    try:
        im.set_array(20*np.log10(np.abs(b[0:200, :])))
    except ZeroDivisionError:
        im.set_array(np.abs(b[0:200, :]))
    i += 1
    N_BUFFERS += 1
    return im

ani = FuncAnimation(fig, animate, interval=10, blit=False)
plt.show()

elapsed = time.time() - start
print('Processed', N_BUFFERS * ALINES_PER_BUFFER, 'A-lines in', elapsed)
print('Virtual A-line rate:', (N_BUFFERS * ALINES_PER_BUFFER) / elapsed, 'hz')

print("imgStopAcq", flush=True)
imgStopAcq()

scan_task.stop()
scan_task.close()

octCleanup()

imgClose()
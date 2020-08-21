# -*- coding: utf-8 -*-

from PyIMAQ import *
import numpy as np
import multiprocessing as mp
from threading import Thread
import matplotlib.pyplot as plt
import matplotlib
import time

# -----------------------------------------------------------------------------

camera_name = "img0"
buffer_size = 24

# -----------------------------------------------------------------------------

print("Frame grab demo from camera",camera_name,"...")

# Open the camera interface
print("imgOpen")
imgShowErrorMsg(imgOpen('img0'))


# Initialize a ring buffer w/ [buffer_size] buffers
print("imgInitBuffers... ",buffer_size,"buffers in ring...")
imgShowErrorMsg(imgInitBuffer(buffer_size))

# Start camera acquisition
print("imgStartAcq")
imgShowErrorMsg(imgStartAcq())

time.sleep(5)

print(imgGetBufferSize(),"bytes/frame...")

fqueue = mp.Queue(buffer_size)

plt.figure(1)
plt.ion()
line, = plt.plot([])
plt.xlim([0, 2048])
plt.ylim([0, 200])
plt.show()

while True:
    
    start = time.time()
    
    fbuff = np.empty(imgGetFrameSize(),dtype=np.uint16)
    
    imgGetCurrentFrame(fbuff,np.empty(1,dtype=np.int))
    
    f = fbuff[2048*100:2048*101]
    del fbuff        
    print(str(1000 / (time.time() - start))[0:6],"hz")
    print(imgGetDroppedFrames(),"dropped frames")
    
    plt.pause(0.001)
    line.set_ydata(f)
    line.set_xdata(range(len(f)))
    
    plt.draw()   

plt.close()

print("imgAbortAcq")
imgShowErrorMsg(imgAbortAcq())

print("imgClose")
imgShowErrorMsg(imgClose())

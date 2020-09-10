import pandas as pd
import numpy as np
import scipy as sp
import pylab as pl
from collections import *
from copy import copy
import os
import subprocess
import numpy.linalg as LA

from statsmodels import robust 
import seaborn as sns
from matplotlib import cm


def smooth(x,w=7):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """
    x=x.as_matrix()
    window_len=w
    window='hamming'

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    
    ## make sure we get the right length of smoothed signal
    
    #Imran's approach:
    #print 'y len:',y.size
    #print 'x len:',x.size
    #print 'windows/2:',window_len/2
    #print 'end:',y.size-x.size-window_len/2
    #print y
    return y[(window_len/2):-(y.size-x.size-window_len/2)]
    
#     # Modefied SciPy Cookbook:
#     if window_len%2==0:
#         #return y[(window_len/2-1):-(window_len/2)]
#         # to be consistent with Imran's starting point
#         return y[(window_len/2):-(window_len/2-1)]
#     else:
#         return y[(window_len/2-1):-(window_len/2+1)]
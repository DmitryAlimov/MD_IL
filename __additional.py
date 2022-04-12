import os, sys, time
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
from scipy.interpolate import interp1d
from scipy import interpolate
from scipy.fft import fft, fftfreq


def smooth(y, d):   ### d: steps for smoothing y
    y_smoothed = np.array(y)
    for i in range(0, len(y)):
        y_smoothed[i] = sum(y[max(0, i - d//2):min(len(y), i + d//2)]) / (-max(0, i - d//2) + min(len(y), i + d//2))
    return y_smoothed
def interp(x, y, x_new):
    s = interpolate.InterpolatedUnivariateSpline(x, y)
    y_new = s(x_new)
    return y_new
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt

def butter_lowpass(cutoff, fs, order=5):
    """
    Create a Butterworth lowpass filter.
    
    Parameters:
    cutoff : float
        The cutoff frequency of the filter.
    fs : float
        The sampling frequency of the data.
    order : int, optional
        The order of the filter (default is 5).
        
    Returns:
    b, a : ndarray, ndarray
        Numerator (b) and denominator (a) polynomials of the IIR filter.
    """
    nyq = 0.5 * fs              # Nyquist Frequency
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    """
    Apply a Butterworth lowpass filter to the data.
    
    Parameters:
    data : array_like
        The input data to be filtered.
    cutoff : float
        The cutoff frequency of the filter.
    fs : float
        The sampling frequency of the data.
    order : int, optional
        The order of the filter (default is 5).
        
    Returns:
    y : ndarray
        The filtered data.
    """
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = filtfilt(b, a, data)
    return y
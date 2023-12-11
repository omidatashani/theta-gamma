function [S] = bandPassFSignal(signal, Fs, Bandpass)
% The aim of this function is to filter the signal into specific
% frequencies
% Input:
% ======
% signal --> the signal which is needed to be filtered
% Fs --> sampling frequency 
% Badpass --> a 1 x 2 vector indicating the low and high values of the
% filter 
% Output:
% =======
% S --> the filtered signal 
% -------------------------------------------------------------------------

Wn_filt                 = [Bandpass(1)/(Fs/2) Bandpass(2)/(Fs/2)]; % normalized frequencies to Nyquist frequency 
[b,a]                   = butter(3,Wn_filt);                       % a third order Butterworth filter 
S                       = filtfilt(b,a,signal);

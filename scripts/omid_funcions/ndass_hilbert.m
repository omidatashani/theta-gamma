% Matlab code for calculation hilbert transform
% 
% This code accompanies the book, titled "Analyzing Neural Time Series Data" 
% (MIT Press). Using the code without following the book may lead to confusion, 
% incorrect data analyses, and misinterpretations of results. 
% 
% inputs:
% sig : time domain of signal for transformation it should be N trials *
% time points
% 
% npref: response signal prefered condition, spike counts in specific time interval correspond to
% each trials
% freqsbands: frequency band for hilbert transform
% fs : sampling frequency 
% output:
% analytic_sig : the complex value analytic signal   
% freqsbands :  frequenny band were used for hilbert transform
% Cohen M.
% Edited by Omid Amir Atashani

function [analytic_sig,freqsbands] = ndass_hilbert(sig,freqsbands,fs)

padding_pnt = 500; %% zero Padding

if isempty(fs)
fs = 1000; % Sampling frequency
end

sig_h = [];
for trii = 1:size(sig,1)
    sig_h = [sig_h [ zeros(1,padding_pnt) (sig(trii,:)) zeros(1,padding_pnt)]];  
end


%% FIR filter  setting

% specify Nyquist freuqency
nyquist = 1000/2;

% filter frequency band
if isempty(freqsbands)
freqsbands=[4 8]; %% Frequency band can be defined based on your focuse in Hertz here
end
% transition width
trans_width = 0.2; % fraction of 1, thus 20%

% filter order
filt_order = 700;%750round(3*(fs/filtbound(1)));700

% frequency vector (as fraction of Nyquist
ffrequencies  = [ 0 (1-trans_width)*freqsbands(1) freqsbands (1+trans_width)*freqsbands(2) nyquist ]/nyquist;

% shape of filter (must be the same number of elements as frequency vector
idealresponse = [ 0 0 1 1 0 0 ];

% get filter weights
filterweights = firls(filt_order,ffrequencies,idealresponse);


%% 
pnts=size(sig,2)+1000;trials=size(sig,1);

s_1_tot=filtfilt(filterweights,1,sig_h);

s_1_tot= hilbert(s_1_tot);

analytic_sig =[]; 

analytic_sig(:,:) = (reshape(s_1_tot,pnts,trials))';

analytic_sig=analytic_sig(:,501:end-500);

end

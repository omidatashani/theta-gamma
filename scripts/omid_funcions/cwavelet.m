% Matlab code for calculation wavelet transform
% inputs:
% sig : time domain of signal for transformation it should be N trials *
% time points
% freqs2use: center frequencies for wavelent transform
% fs : sampling frequency 
% output:
% analytic_sig : the complex value analytic signal   
% freqs2use : center frequencies were used for wavelent transform
% Edited by Omid Amir Atashani

function [analytic_sig, freqs2use] = cwavelet(sig, freqs2use, fs)

    if isempty(freqs2use)
        freqs2use=[4:30, 33:3:130];
    end

    pnts = length(sig);
    time = -1:1/fs:1;
    half_wavelet = (length(time)-1)/2;

    num_cycles    = 7 * ones(1, length(freqs2use));
    n_wavelet     = length(time);
    
    analytic_sig_h = zeros(pnts, length(freqs2use)); % Pre-allocate memory for speed

    %% core
    for fi = 1:length(freqs2use)
        % create wavelet
        s = num_cycles(fi) / (2 * pi * freqs2use(fi));
        wavelet_time = exp(2 * 1i * pi * freqs2use(fi) .* time) .* exp(-time.^2./(2*(s^2)));
        
        % phase angles via convolution
        convolution_result_time = conv(sig, wavelet_time(end:-1:1), 'same');  % note 'same' to keep the output size equal to input
        
        % Store the result in analytic_sig_h
        analytic_sig_h(:, fi) = convolution_result_time;
    end % end frequency loop

    analytic_sig = analytic_sig_h';
end

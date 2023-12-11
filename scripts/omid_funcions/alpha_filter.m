function filtered = alpha_filter(sig,fs)

nyquist = fs/2;
freqsbands=[8 13];

trans_width = 0.2;

filt_order = 950;

ffrequencies  = [ 0 (1-trans_width)*freqsbands(1) freqsbands (1+trans_width)*freqsbands(2) nyquist ]/nyquist;

idealresponse = [ 0 0 1 1 0 0 ];

filterweights = firls(filt_order,ffrequencies,idealresponse);

filtered=filtfilt(filterweights,1,sig);
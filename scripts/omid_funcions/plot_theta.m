clc
clear
%%
lfp1 = HPC;
lfp2 = PFC;

[analytic_sig1] = theta_hilbert(lfp1);
plot(abs(analytic_sig1(:)))
%%
nyquist = 1000/2;

freqsbands=[6 12]; 

trans_width = 0.2;

filt_order = 1100;

ffrequencies  = [ 0 (1-trans_width)*freqsbands(1) freqsbands (1+trans_width)*freqsbands(2) nyquist ]/nyquist;

idealresponse = [ 0 0 1 1 0 0 ];

filterweights = firls(filt_order,ffrequencies,idealresponse);

lfp1m = [lfp1(1,:),lfp1(2,:)];

lfp1f=filtfilt(filterweights,1,lfp1m);
figure
plot(lfp1f(1,:))
xlim([1 2800])
%%
figure
plot(angle(analytic_sig1(1,:)))
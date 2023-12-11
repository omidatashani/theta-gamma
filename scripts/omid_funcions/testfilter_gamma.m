% This function test FIR filter
%
% Edited  By Omid Amir Atashani
function testfilter_gamma()
nyquist = 500;
lower_filter_bound = 30; % Hz
upper_filter_bound = 150; % Hz
transition_width   = 0.1;
filter_order       =120;

ffrequencies  = [ 0 (1-transition_width)*lower_filter_bound lower_filter_bound upper_filter_bound (1+transition_width)*upper_filter_bound nyquist ]/nyquist;
idealresponse = [ 0 0 1 1 0 0 ];
filterweights = firls(filter_order,ffrequencies,idealresponse);

fft_filtkern  = abs(fft(filterweights));

hz_filtkern   = linspace(0,nyquist,ceil(length(fft_filtkern)/2)); % list of frequencies in Hz corresponding to filter kernel

figure
plot(ffrequencies*nyquist,idealresponse,'r')
hold on

fft_filtkern  = abs(fft(filterweights));
fft_filtkern  = fft_filtkern./max(fft_filtkern); % normalized to 1.0 for visual comparison ease
plot(hz_filtkern,fft_filtkern(1:ceil(length(fft_filtkern)/2)),'b')

set(gca,'ylim',[-.1 1.1],'xlim',[0 250])
legend({'ideal';'best fit'})

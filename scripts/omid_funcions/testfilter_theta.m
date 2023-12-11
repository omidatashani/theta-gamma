% This function test FIR filter
%
% Edited By Omid Amir Atashani
function testfilter_theta()
nyquist = 500;
lower_filter_bound = 6; % Hz
upper_filter_bound = 12; % Hz
transition_width   = 0.2;
filter_order       = 1100;

ffrequencies  = [ 0 (1-transition_width)*lower_filter_bound lower_filter_bound upper_filter_bound (1+transition_width)*upper_filter_bound nyquist ]/nyquist;
idealresponse = [ 0 0 1 1 0 0 ];
filterweights = firls(filter_order,ffrequencies,idealresponse);

% Compute the frequency response of the filter
[H,w] = freqz(filterweights, 1, 'half', nyquist*2); % 'half' computes the response at frequencies around [0, nyquist]
f = w/(2*pi) * nyquist*2;  % Convert frequency points from rad/sample to Hz

% Normalize the actual filter response for visual comparison ease
H_normalized = abs(H)./max(abs(H));

% Create the figure
figure;

% Plot magnitude response
plot(ffrequencies*nyquist, idealresponse, 'r'); % Plot ideal response
hold on;
plot(f, H_normalized, 'b'); % Plot actual filter response
set(gca,'ylim',[-.1 1.1],'xlim',[0 25]);
legend({'Ideal', 'Actual'});
title('Filter Magnitude Response');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Plot phase response
figure
phase_response = rad2deg(unwrap(angle(H)));  % Compute phase response, unwrap it, and convert to degrees
plot(f, phase_response, 'b');
set(gca,'ylim',[-1400 0],'xlim',[lower_filter_bound upper_filter_bound]);
title('Filter Phase Response');
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
grid on;

end



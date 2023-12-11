function plotFilterResponse(signal, Fs, Bandpass)
    % Your bandPassFSignal function
    [S] = bandPassFSignal(signal, Fs, Bandpass);

    % Compute the filter coefficients for the Butterworth filter
    Wn_filt = [Bandpass(1)/(Fs/2) Bandpass(2)/(Fs/2)]; % normalized frequencies to Nyquist frequency 
    [b,a] = butter(3,Wn_filt);  % a third order Butterworth filter 

    % Compute the frequency response of the filter
    [H,w] = freqz(b, a, 'half', Fs); % 'half' computes the response at 4096 frequencies around [0, Fs/2]
    f = w/(2*pi) * Fs;  % Convert frequency points from rad/sample to Hz

    % Compute the ideal response
    lower_filter_bound = Bandpass(1);
    upper_filter_bound = Bandpass(2);
    transition_width = 0.2; % Example transition width, modify as needed
    ffrequencies = [0 (1-transition_width)*lower_filter_bound lower_filter_bound upper_filter_bound (1+transition_width)*upper_filter_bound Fs/2]/(Fs/2);
    idealresponse = [0 0 1 1 0 0];

    % Normalize the actual filter response for visual comparison ease
    H_normalized = abs(H)./max(abs(H));

    % Create the figure
    figure;
    
    % Plot magnitude response
    
    plot(ffrequencies*(Fs/2), idealresponse, 'r'); % Plot ideal response
    hold on;
    plot(f, H_normalized, 'b'); % Plot actual filter response
    set(gca,'ylim',[-.1 1.1],'xlim',[(1) (upper_filter_bound + 45)]); % Adjust the x-limits to show the interval of interest and a bit before and after
    legend({'Ideal', 'Actual'});
    title('Filter Magnitude Response');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    grid on;
    
    % Plot phase response
    figure
    phase_response = rad2deg(angle(H));  % Compute phase response and convert to degrees
    plot(f, phase_response, 'b');
    set(gca,'ylim',[-200 200],'xlim',[lower_filter_bound upper_filter_bound]);
    title('Filter Phase Response');
    xlabel('Frequency (Hz)');
    ylabel('Phase (degrees)');
    grid on;
    
end



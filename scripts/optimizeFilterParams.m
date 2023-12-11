function [bestOrder, bestWidth, minError] = optimizeFilterParams(Fs, Bandpass)
    minError = Inf;  % Initialize minimum error to infinity
    bestOrder = [];
    bestWidth = [];

    % Define ranges of parameters to test
    orderRange = 2:2:10;  % Example: test filter orders from 2 to 10 in steps of 2
    widthRange = 0.1:0.1:0.5;  % Example: test transition widths from 0.1 to 0.5 in steps of 0.1

    for order = orderRange
        for width = widthRange
            % Create Butterworth filter for current parameters
            Wn_filt = [(Bandpass(1)-width)/(Fs/2), (Bandpass(2)+width)/(Fs/2)];
            [b, a] = butter(order, Wn_filt);
            
            % Compute frequency response of current filter
            [H, w] = freqz(b, a, 'half', Fs);
            f = w/(2*pi) * Fs;
            
            % Create ideal response at the same frequency points as the actual response
            idealresponse = zeros(size(f));
            idealresponse(f >= (Bandpass(1)-width) & f <= Bandpass(1)) = 1;
            idealresponse(f >= Bandpass(2) & f <= (Bandpass(2)+width)) = 1;
            
            % Compute mean squared error between actual and ideal responses
            mse = mean((abs(H) - idealresponse).^2);
            
            % Update best parameters if current error is lower than minimum error
            if mse < minError
                minError = mse;
                bestOrder = order;
                bestWidth = width;
            end
        end
    end
end


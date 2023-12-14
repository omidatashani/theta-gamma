clc
clear
%%
cd '/Users/omid/Desktop/scripts'
addpath(genpath("FMAToolbox-master/"))
addpath("omid_funcions/")
load('/Users/omid/Desktop/2021-06-03_13-38-14_posttrial5/HPC_100_CH18_0.continuous.mat')
load('/Users/omid/Desktop/2021-06-03_13-38-14_posttrial5/PFC_100_CH22_0.continuous.mat')
load("/Users/omid/Desktop/2021-06-03_13-38-14_posttrial5/2021-06-03_13-38-14_posttrial5-states_SM.mat")
%%
samplingRate = 2500;
[HPClfpCleaned,HPCnoisyInds] = removeArtefacts(HPC,samplingRate,[4 8], ...
    [1 0.1]);
[PFClfpCleaned,PFCnoisyInds] = removeArtefacts(PFC,samplingRate,[4 8], ...
    [1 0.1]);
%%
thetaband=[6 12];
samplingRate = 1250;
%%
% LFPtheta = bandPassFSignal(HPClfpCleaned,samplingRate,thetaband);
% thetaPhase=angle(hilbert(LFPtheta));
%%
LFPtheta = theta_filter(HPClfpCleaned,samplingRate);
thetaPhaseHPC=angle(hilbert(LFPtheta));
LFPtheta = theta_filter(PFClfpCleaned,samplingRate);
thetaPhasePFC=angle(hilbert(LFPtheta));
%%
% figure
% plot(thetaPhase)
%%
% plot(LFPtheta)
%%
% plotFilterResponse(HPClfpCleaned,samplingRate,thetaband)
% %%
% testfilter_theta()
%%
cuttime = 10801;
tic
[wtH,f] = cwt(HPClfpCleaned(1:cuttime),'amor',samplingRate,'FrequencyLimits',[20 180]);
toc
[wtP,f] = cwt(PFClfpCleaned(1:cuttime),'amor',samplingRate,'FrequencyLimits',[20 180]);
figure
h=pcolor(1:cuttime,f,abs(wtH).^2);
h.EdgeColor = 'none';
colorbar
title('Raw Wavelet Power (cwt)')
xlabel('Time')
ylabel('Frequency (Hz)')
%%
% num_octaves = log2(180/20);
% desired_frequencies = 81; % from 20 to 180 with 2 Hz increment
% voices_per_octave = round(desired_frequencies / num_octaves);
% 
% [wtP, f] = cwt(PFClfpCleaned(1:cuttime), 'amor', samplingRate, 'FrequencyLimits', [20 180], 'VoicesPerOctave', voices_per_octave);
% %%
% % Your obtained wavelet transform and frequency vector
% [wtP, f] = cwt(PFClfpCleaned(1:cuttime), 'amor', samplingRate, 'FrequencyLimits', [20 180]);
% 
% % Desired frequency vector
% desired_frequencies = 20:2:180;
% 
% % Interpolate the wavelet transform coefficients at the desired frequencies
% wtP_interpolated = interp1(f, wtP, desired_frequencies, 'linear', 'extrap');
% 
% % Now, wtP_interpolated contains the wavelet transform coefficients at your desired frequencies.
% 
% figure
% h=pcolor(1:cuttime,desired_frequencies,abs(wtP_interpolated).^2);
% h.EdgeColor = 'none';
% colorbar
% title('Raw Wavelet Power (cwt)')
% xlabel('Time')
% ylabel('Frequency (Hz)')
%%
% % Smoothing parameters
% time_smoothing_length = round(8e-3 * samplingRate) * 2 + 1; % For ±8ms in time
% frequency_smoothing_length = 3;  % For ±2Hz in frequency since you now have 2Hz steps
% 
% % Use the provided smoothing function
% flag_SMOOTH = true;
% wtP_smooth = smoothCFSAbdel(wtP_interpolated, flag_SMOOTH, frequency_smoothing_length, time_smoothing_length);
% 
% figure
% h=pcolor(1:cuttime,desired_frequencies,abs(wtP_smooth).^2);
% h.EdgeColor = 'none';
% colorbar
% title('Raw Wavelet Power (cwt)')
% xlabel('Time')
% ylabel('Frequency (Hz)')
%%
log10_power = log10(abs(wtH).^2 + eps);  % Adding eps to avoid log(0)
%%
power_spectrum = abs(wtH).^2;
mean_power_time = mean(power_spectrum, 2);
std_power_time = std(power_spectrum, 0, 2);
z_scored_power_time = (power_spectrum - mean_power_time) ./ std_power_time;
%%
mean_power_freq = mean(power_spectrum, 1);
std_power_freq = std(power_spectrum, 0, 1);
z_scored_power_freq = (power_spectrum - mean_power_freq) ./ std_power_freq;
%%
total_power = sum(abs(wtH).^2, 'all');
normalized_power = abs(wtH).^2 / total_power;
%% Plotting log10 normalized power spectrum
figure;
h = pcolor(1:cuttime, f, log10_power);
h.EdgeColor = 'none';
colorbar;
title('log10 Normalized Wavelet Power (cwt)');
xlabel('Time');
ylabel('Frequency (Hz)');
%%
figure;
h = pcolor(1:cuttime, f, z_scored_power_time);
h.EdgeColor = 'none';
colorbar;
title('z scored power over t Normalized Wavelet Power (cwt) from paper');
xlabel('Time');
ylabel('Frequency (Hz)');
%%
figure;
h = pcolor(1:cuttime, f, z_scored_power_freq);
h.EdgeColor = 'none';
colorbar;
title('z scored power over f Normalized Wavelet Power (cwt)');
xlabel('Time');
ylabel('Frequency (Hz)');
%%
figure;
h = pcolor(1:cuttime, f, normalized_power);
h.EdgeColor = 'none';
colorbar;
title('Normalized Wavelet Power to sum (cwt)');
xlabel('Time');
ylabel('Frequency (Hz)');
%%
num_bins = 20;  %  20 bins, yielding 18-degree bins
bin_edges = linspace(-pi, pi, num_bins + 1);
[~,~,bin_idx] = histcounts(thetaPhaseHPC, bin_edges);
bin_idx = bin_idx(1:cuttime);
PPC = zeros(num_bins, length(f));
for freq_idx = 1:length(f)
    for bin = 1:num_bins
        idx = bin_idx == bin;
        phase_angles = angle(wtP(freq_idx, idx));  % only take phase angles for the current frequency
        N = sum(idx);  % number of time points in this bin
        ppc_sum = 0;
        for n = 1:N
            for m = n+1:N
                ppc_sum = ppc_sum + cos(phase_angles(n) - phase_angles(m));
            end
        end
        PPC(bin, freq_idx) = ppc_sum / (N * (N - 1));  % store the PPC value for the current frequency
    end
end
phase_bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
figure;
h = pcolor(phase_bin_centers, f, PPC');
h.EdgeColor = 'none';
%shading interp;
colorbar;
xlabel('Theta Phase (radians)');
ylabel('Frequency (Hz)');
title('Pairwise Phase Consistency');
%%
num_bins = 20;  
bin_edges = linspace(-pi, pi, num_bins + 1);
[~,~,bin_idx] = histcounts(thetaPhaseHPC, bin_edges);
bin_idx = bin_idx(1:cuttime);

PPC = zeros(num_bins, length(f));

for freq_idx = 1:length(f)
    for bin = 1:num_bins
        idx = bin_idx == bin;
        phase_angles = angle(wtP(freq_idx, idx));  % only take phase angles for the current frequency
        PPC(bin, freq_idx) = ppc(phase_angles, 2);  
    end
end
phase_bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
figure;
h = pcolor(phase_bin_centers, f, PPC');
h.EdgeColor = 'none';
%shading interp;
colorbar;
xlabel('Theta Phase (radians)');
ylabel('Frequency (Hz)');
title('Pairwise Phase Consistency');

%%
% [wcoh, wcs, f] = wcoherence(HPClfpCleaned(1:cuttime),PFClfpCleaned(1:cuttime), ...
%     samplingRate,'FrequencyLimits',[20 180]);
% figure
% h=pcolor(1:cuttime,f,wcoh);
% h.EdgeColor = 'none';
% colorbar
% title('Wavelet Coherence')
% xlabel('Time')
% ylabel('Frequency (Hz)')
%%
% [sample, ThetaTS] = lu_wcsthetaextract(HPClfpCleaned, PFClfpCleaned, ...
%     LFPtheta, thetaPhase, samplingRate, [1; 1201]);
% %%
% selection = 600;
% 
% figure
% for i=1:10
%     nexttile(i)
%     plot(sample(:, selection+i));
%     title(string(selection+i))
% end
%%
tic
[wtm,fm] = cwavelet(HPClfpCleaned(1:cuttime), 20:2:180, 1250);
toc
wpow = abs(wtm).^2;
figure
h=pcolor(1:cuttime,20:2:180,wpow);
h.EdgeColor = 'none';
colorbar
title('Raw Wavelet Power (mine)')
xlabel('Time')
ylabel('Frequency (Hz)')
%%
log10_powerm = log10(abs(wtm).^2 + eps);  % Adding eps to avoid log(0)
%%
power_spectrum = abs(wtm).^2;
mean_power_time = mean(power_spectrum, 2);
std_power_time = std(power_spectrum, 0, 2);
z_scored_power_timem = (power_spectrum - mean_power_time) ./ std_power_time;
%%
mean_power_freq = mean(power_spectrum, 1);
std_power_freq = std(power_spectrum, 0, 1);
z_scored_power_freqm = (power_spectrum - mean_power_freq) ./ std_power_freq;
%%
total_power = sum(abs(wtm).^2, 'all');
normalized_powerm = abs(wtm).^2 / total_power;
%%
figure;
h = pcolor(1:cuttime, fm, log10_powerm);
h.EdgeColor = 'none';
colorbar;
title('log10 Normalized Wavelet Power (mine)');
xlabel('Time');
ylabel('Frequency (Hz)');
%%
figure;
h = pcolor(1:cuttime, fm, z_scored_power_timem);
h.EdgeColor = 'none';
colorbar;
title('z scored power over t Normalized Wavelet Power (mine) from paper');
xlabel('Time');
ylabel('Frequency (Hz)');
%%
figure;
h = pcolor(1:cuttime, fm, z_scored_power_freqm);
h.EdgeColor = 'none';
colorbar;
title('z scored power over f Normalized Wavelet Power (mine)');
xlabel('Time');
ylabel('Frequency (Hz)');
%%
figure;
h = pcolor(1:cuttime, fm, normalized_powerm);
h.EdgeColor = 'none';
colorbar;
title('Normalized Wavelet Power to sum (mine)');
xlabel('Time');
ylabel('Frequency (Hz)');
%%
% params.Fs=1250; % sampling frequency
% params.fpass=[1 200]; % band of frequencies to be kept
% params.tapers=[2 3];
% params.pad=2; % pad factor for fft
% params.err=[0 0.05];
% params.trialave=0;
% movingwin=[0.2 0.001];
% 
% [C1,phi1,S12,S1,S2,t_mt,f_mt] = cohgramc(HPClfpCleaned(1:cuttime),PFClfpCleaned(1:cuttime) ...
%     ,movingwin,params);
% figure
% h=pcolor(t_mt,f_mt,C1');
% h.EdgeColor = 'none';
% colorbar
% title('Coherence')
% xlabel('Time')
% ylabel('Frequency (Hz)')
%%
rem_states = find(states==5);

% Finding the discontinuities in the REM states indices
discont = find(diff(rem_states) > 1);

start_indices = [rem_states(1), rem_states(discont + 1)];
end_indices = [rem_states(discont), rem_states(end)];

rem_LFPHPC = cell(length(start_indices), 1);
rem_LFPPFC = cell(length(start_indices), 1);

for i = 1:length(start_indices)

    start_idx = (start_indices(i) - 1) * 1250 + 1;
    end_idx = end_indices(i) * 1250;

    rem_LFPPFC{i} = PFClfpCleaned(start_idx:end_idx);
    rem_LFPHPC{i} = HPClfpCleaned(start_idx:end_idx);

end
rem_LFPHPC{2} = rem_LFPHPC{2}(1:76000);
%%
figure
plot(rem_LFPHPC{2})
title('Raw LFP HPC trial 2');
%%
samplingRate = 1250;
num_bins = 20;
bin_edges = linspace(-pi, pi, num_bins + 1);
% Smoothing parameters
time_smoothing_length = round(8e-3 * samplingRate) * 2 + 1; % For ±8ms in time
frequency_smoothing_length = 3;  % For ±2Hz in frequency since you now have 2Hz steps
flag_SMOOTH = true;


PPC = cell(length(rem_LFPHPC),1);  
tic
for i = 1:length(rem_LFPHPC)

    signalHPC = rem_LFPHPC{i};

    LFPtheta = theta_filter(signalHPC, samplingRate);
    
    thetaPhase = angle(hilbert(LFPtheta));

    [~,~,bin_idx] = histcounts(thetaPhase, bin_edges);
    bin_idx = bin_idx(1:length(signalHPC));  % Ensure bin_idx is the same length as the signal
    
    signalPFC = rem_LFPPFC{i};

    [wt, f] = cwt(signalPFC, 'amor', samplingRate, 'FrequencyLimits', [20 180]);

    %[wt,f] = cwavelet(signalPFC, 20:2:180, 1250);
    %wt = smoothCFSAbdel(wt, flag_SMOOTH, frequency_smoothing_length, time_smoothing_length);
    
    PPC{i} = zeros(num_bins, length(f));
    
    for freq_idx = 1:length(f)
        for bin = 1:num_bins
            idx = bin_idx == bin;
            phase_angles = angle(wt(freq_idx, idx));  % Only take phase angles for the current frequency
            PPC{i}(bin, freq_idx) = ppc(phase_angles, 2);  
        end
    end
end
toc
%%
x_labels = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;  % Compute bin centers

for i = 1:length(PPC)
    figure;

    h = pcolor(x_labels, f, PPC{i}');
    h.EdgeColor = 'none';
    
    %shading interp;  % Interpolate grid for smoother color

    colorbar;
    xlabel('HPC Theta Phase (radians)');
    ylabel('PFC Frequency (Hz)');
    title(['Pairwise Phase Consistency for REM trial ' num2str(i) ' (my function)']);
end
%%
ppc2 = PPC{1}+ PPC{2};
ppc2=ppc2/2;
h = pcolor(x_labels, f, ppc2');
h.EdgeColor = 'none';
shading interp;  % Interpolate grid for smoother color
colorbar;
xlabel('HPC Theta Phase (radians)');
ylabel('PFC Frequency (Hz)');
title(['Pairwise Phase Consistency averaged in trials (my function)']);
%%
inds = [1:7, 18:20];
PPC_matrix = [];
for i = 1:length(PPC)
    PPC_matrix = [PPC_matrix; mean(PPC{i}(inds, :), 1)];  % averaging across phase bin indices
end

c1 = 1; 
c2 = 0;
c3 = 0; 

sem_plot(PPC_matrix,f,c1, c2, c3)

xlabel('Frequency (Hz)');
ylabel('Average PPC');
title('Average PPC with SEM');
%%
win = 1; 

semplot(PPC_matrix, f', win, c1, c2, c3);
xlabel('Frequency (Hz)');
ylabel('Average PPC');
title('Average PPC with SEM (PFC)');
grid on;
%%
samplingRate = 1250;
num_bins = 20;
bin_edges = linspace(-pi, pi, num_bins + 1);

Power = cell(2,1);  % Store power values

for i = 1:2
    signalHPC = rem_LFPHPC{i};
    LFPtheta = theta_filter(signalHPC, samplingRate);
    thetaPhase = angle(hilbert(LFPtheta));
    [~,~,bin_idx] = histcounts(thetaPhase, bin_edges);
    bin_idx = bin_idx(1:length(signalHPC));  % Ensure bin_idx is the same length as the signal

    signalPFC = rem_LFPPFC{i};
    %[wt, f] = cwt(signalPFC, 'amor', samplingRate, 'FrequencyLimits', [20 180]);
    [wt,f] = cwavelet(signalPFC, 20:2:180, 1250);
    % Z-score normalization (Time)
    mean_power_time = mean(wt, 1);
    std_power_time = std(wt, 0, 1);
    NormalizedPower_ztime = (wt - mean_power_time) ./ std_power_time;
    
    % Z-score normalization (Frequency)
    mean_power_freq = mean(wt, 2);
    std_power_freq = std(wt, 0, 2);
    NormalizedPower_zfreq = (wt - mean_power_freq) ./ std_power_freq;

    wt = smoothCFSAbdel(NormalizedPower_zfreq, flag_SMOOTH, frequency_smoothing_length, time_smoothing_length);
    
    Power{i} = zeros(num_bins, length(f));  % Initialize power matrix
    
    for freq_idx = 1:length(f)
        for bin = 1:num_bins
            idx = bin_idx == bin;
            % Compute power for the current frequency and bin
            Power{i}(bin, freq_idx) = mean(abs(wt(freq_idx, idx)).^2);
        end
    end
end
%%
x_labels = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;  % Compute bin centers

for i = 1:length(Power)
    figure;
    h = pcolor(x_labels, f, Power{i}');
    h.EdgeColor = 'none';
    shading interp;  % Interpolate grid for smoother color
    colorbar;
    xlabel('HPC Theta Phase (radians)');
    ylabel('PFC Frequency (Hz)');
    title(['Power for REM trial ' num2str(i) ' (cwt)']);
end
%%
samplingRate = 1250;
num_bins = 20;
bin_edges = linspace(-pi, pi, num_bins + 1);

Power = cell(2,1);  
NormalizedPower_log10 = cell(2,1);  
NormalizedPower_zphase = cell(2,1);  
NormalizedPower_zfrequ = cell(2,1);  
NormalizedPower_total = cell(2,1);  

for i = 1:2
    signalHPC = rem_LFPHPC{i};
    LFPtheta = theta_filter(signalHPC, samplingRate);
    thetaPhase = angle(hilbert(LFPtheta));
    [~,~,bin_idx] = histcounts(thetaPhase, bin_edges);
    bin_idx = bin_idx(1:length(signalHPC));  % Ensure bin_idx is the same length as the signal

    signalPFC = rem_LFPPFC{i};
    %[wt, f] = cwt(signalPFC, 'amor', samplingRate, 'FrequencyLimits', [20 180]);
    [wt,f] = cwavelet(signalPFC, 20:2:180, 1250);

    % Z-score normalization (Time)
    mean_power_time = mean(wt, 1);
    std_power_time = std(wt, 0, 1);
    NormalizedPower_ztime = (wt - mean_power_time) ./ std_power_time;
    
    % Z-score normalization (Frequency)
    mean_power_freq = mean(wt, 2);
    std_power_freq = std(wt, 0, 2);
    NormalizedPower_zfreq = (wt - mean_power_freq) ./ std_power_freq;

    wt = smoothCFSAbdel(NormalizedPower_ztime, flag_SMOOTH, frequency_smoothing_length, time_smoothing_length);
    
    Power{i} = zeros(num_bins, length(f));  % Initialize power matrix
    
    for freq_idx = 1:length(f)
        for bin = 1:num_bins
            idx = bin_idx == bin;
            Power{i}(bin, freq_idx) = mean(abs(wt(freq_idx, idx)).^2);
        end
    end
    
    % Log10 Normalization
    NormalizedPower_log10{i} = log10(Power{i} + eps);
    
    % Z-score normalization (Phase)
    mean_power_phase = mean(Power{i}, 1);
    std_power_phase = std(Power{i}, 0, 1);
    NormalizedPower_zphase{i} = (Power{i} - mean_power_phase) ./ std_power_phase;
    
    % Z-score normalization (Frequency)
    mean_power_freq = mean(Power{i}, 2);
    std_power_freq = std(Power{i}, 0, 2);
    NormalizedPower_zfrequ{i} = (Power{i} - mean_power_freq) ./ std_power_freq;
    
    % Total power normalization
    total_power = sum(Power{i}, 'all');
    NormalizedPower_total{i} = Power{i} / total_power;
end

%%
x_labels = (bin_edges(1:end-1) + bin_edges(2:end)) / 2; 

for i = 1:length(Power)
    figure;
    h = pcolor(x_labels, f, NormalizedPower_zphase{i}');  
    h.EdgeColor = 'none';
    %shading interp; 
    colorbar;
    xlabel('HPC Theta Phase (radians)');
    ylabel('PFC Frequency (Hz)');
    title(['z scored across phases Power for REM trial ' num2str(i)]);
end
%%
figure
NormalizedPower_zphase2 = NormalizedPower_zphase{1} + NormalizedPower_zphase{2};
NormalizedPower_zphase2=NormalizedPower_zphase2/2;
h = pcolor(x_labels, f, NormalizedPower_zphase2');
h.EdgeColor = 'none';
shading interp;
colorbar;
xlabel('HPC Theta Phase (radians)');
ylabel('PFC Frequency (Hz)');
title('z scored across phases Power averaged in trials');
%%
% weighted_phases = bsxfun(@times, NormalizedPower_zphase2, exp(1i * x_labels'));
% sum_power = sum(NormalizedPower_zphase2, 2); % Sum power across frequency bins for each phase bin
% gravity_phase_center = angle(sum(weighted_phases, 2) ./ sum_power); % Gravity center for each phase bin
% 
% % Calculate the gravity center for frequency
% weighted_frequencies = NormalizedPower_zphase2 .* f; % Element-wise multiplication
% gravity_frequency_center = sum(weighted_frequencies(:)) / sum(NormalizedPower_zphase2(:)); 
% 
% % Phase standard deviation calculation for each phase bin
% phase_diffs = bsxfun(@minus, x_labels', gravity_phase_center); % Phase differences
% phase_diffs = wrapToPi(phase_diffs); % Correct phase differences for circular data
% phase_squared_diffs = NormalizedPower_zphase2 .* (phase_diffs.^2);
% phase_std_dev = sqrt(sum(phase_squared_diffs, 2) ./ sum_power);
% 
% % Calculate the overall gravity phase center and standard deviation
% overall_gravity_phase_center = angle(mean(exp(1i * gravity_phase_center)));
% overall_phase_std_dev = mean(phase_std_dev);
% 
% % Define the interval for averaging PPC
% interval_start = overall_gravity_phase_center - 7 * overall_phase_std_dev;
% interval_end = overall_gravity_phase_center + overall_phase_std_dev;
% 
% % Normalize the interval to be within -pi to pi
% interval_start = wrapToPi(interval_start);
% interval_end = wrapToPi(interval_end);
% 
% % Display the results
% fprintf('Overall Gravity Phase Center: %.4f radians\n', overall_gravity_phase_center);
% fprintf('Overall Phase Standard Deviation: %.4f radians\n', overall_phase_std_dev);
% fprintf('Gravity Frequency Center: %.4f Hz\n', gravity_frequency_center);
% fprintf('Interval Start: %.4f radians\n', interval_start);
% fprintf('Interval End: %.4f radians\n', interval_end);
% 
% 
% 
% %%
% threshold_percentile = 75;
% 
% for i = 1:2
%     HPC_LFP = rem_LFPHPC{i};
%     PFC_LFP = rem_LFPPFC{i};
% 
%     LFPtheta = theta_filter(HPC_LFP, samplingRate);
%     thetaPhaseHPC = angle(hilbert(LFPtheta)); 
% 
%     [HPC_peaks, HPC_peak_locs] = findpeaks(LFPtheta);
% 
%     power_threshold = prctile(HPC_peaks, 100 - threshold_percentile);
%     selected_HPC_cycles = {};
%     selected_PFC_cycles = {};
% 
%     for k = 1:length(HPC_peaks)
%         if HPC_peaks(k) > power_threshold
%             peak_loc = HPC_peak_locs(k);
%             
%             % Find the phase change just before the peak (cycle start)
%             phase_diff = diff([0; thetaPhaseHPC(1:peak_loc)]);  % Include zero to align with thetaPhaseHPC indices
%             potential_starts = find(phase_diff < -pi);
%             if isempty(potential_starts)
%                 cycle_start_idx = 1;  % If no wrapping occurs, start from the first index
%             else
%                 cycle_start_idx = potential_starts(end) + 1;  % Take the last wrapping point before the peak
%             end
%             
%             % Find the phase change just after the peak (cycle end)
%             phase_diff = diff([thetaPhaseHPC(peak_loc:end); 0]);  % Include zero to align with thetaPhaseHPC indices
%             potential_ends = find(phase_diff < -pi);
%             if isempty(potential_ends)
%                 cycle_end_idx = length(HPC_LFP);  % If no wrapping occurs, end at the last index
%             else
%                 cycle_end_idx = potential_ends(1) + peak_loc;  % Take the first wrapping point after the peak
%             end
% 
%             % Extract the cycle based on the identified start and end indices
%             selected_HPC_cycles{end+1} = HPC_LFP(cycle_start_idx:cycle_end_idx);
%             selected_PFC_cycles{end+1} = PFC_LFP(cycle_start_idx:cycle_end_idx);
%         end
%     end
% 
%     selected_cycles_HPC{i} = selected_HPC_cycles;
%     selected_cycles_PFC{i} = selected_PFC_cycles;
% end
% %%
% samplingRate = 1250;
% num_bins = 20;
% bin_edges = linspace(-pi, pi, num_bins + 1);
% frequencies = 20:2:180; 
% x_labels = (bin_edges(1:end-1) + bin_edges(2:end)) / 2; 
% Bandpass = [6 12];  % Theta band
% 
% NormalizedPower_log10 = cell(2,1);
% NormalizedPower_zphase = cell(2,1);
% NormalizedPower_zfreq = cell(2,1);
% 
% for i = 1:2
%     Power = []; 
%     for cycle_idx = 1:length(selected_cycles_HPC{i})
%         signalHPC = selected_cycles_HPC{i}{cycle_idx};
% %         LFPtheta = theta_filter(signalHPC, samplingRate);
%         LFPtheta = bandPassFSignal(signalHPC, samplingRate, Bandpass);
%         thetaPhase = angle(hilbert(LFPtheta));
%         [~,~,bin_idx] = histcounts(thetaPhase, bin_edges);
%         bin_idx = bin_idx(1:length(signalHPC)); % Ensure bin_idx is the same length as the signal
% 
%         signalPFC = selected_cycles_PFC{i}{cycle_idx};
%         [wt, f] = cwavelet(signalPFC, frequencies, samplingRate);
%         wt = smoothCFSAbdel(wt, flag_SMOOTH, frequency_smoothing_length, time_smoothing_length);
% 
%         Power_cycle = zeros(num_bins, length(f));  
%         for freq_idx = 1:length(f)
%             for bin = 1:num_bins
%                 idx = bin_idx == bin;
%                 Power_cycle(bin, freq_idx) = mean(abs(wt(freq_idx, idx)).^2);
%             end
%         end
%         Power = cat(3, Power, Power_cycle); % Concatenate the power matrix of this cycle to the rest
%     end
%     
%     % Average the power over all cycles for this trial
%     Power_avg = nanmean(Power, 3);
%     
%     % Log10 Normalization
%     NormalizedPower_log10{i} = log10(Power_avg + eps);
%     
%     % Z-score normalization (Phase)
%     mean_power_phase = mean(Power_avg, 1);
%     std_power_phase = std(Power_avg, 0, 1);
%     NormalizedPower_zphase{i} = (Power_avg - mean_power_phase) ./ std_power_phase;
%     
%     % Z-score normalization (Frequency)
%     mean_power_freq = mean(Power_avg, 2);
%     std_power_freq = std(Power_avg, 0, 2);
%     NormalizedPower_zfreq{i} = (Power_avg - mean_power_freq) ./ std_power_freq;
% end
% 
% NormalizedPower_zphase_avg = mean(cat(3, NormalizedPower_zphase{:}), 3);
% 
% figure;
% h = pcolor(x_labels, f, NormalizedPower_zphase_avg');
% h.EdgeColor = 'none';
% %shading interp; 
% colorbar;
% xlabel('HPC Theta Phase (radians)');
% ylabel('PFC Frequency (Hz)');
% title('Z-Scored Power Across Phases in High Amp Cycles');
%%
% tic
% PPC = cell(2,1);  
% 
% for i = 1:2
%     PPC{i} = zeros(num_bins, length(frequencies)); 
% 
%     all_phases = cell(num_bins, length(frequencies));
% 
%     for cycle_idx = 1:length(selected_cycles_HPC{i})
%         signalHPC = selected_cycles_HPC{i}{cycle_idx};
%         %LFPtheta = theta_filter(signalHPC, samplingRate);
%         LFPtheta = bandPassFSignal(signalHPC, samplingRate, Bandpass);
%         thetaPhase = angle(hilbert(LFPtheta));
% 
%         signalPFC = selected_cycles_PFC{i}{cycle_idx};
%         [wt, f] = cwavelet(signalPFC, frequencies, samplingRate);
%         wt = smoothCFSAbdel(wt, flag_SMOOTH, frequency_smoothing_length, time_smoothing_length);
% 
%         [~, ~, bin_idx] = histcounts(thetaPhase, bin_edges);
% 
%         for freq_idx = 1:length(frequencies)
%             for bin = 1:num_bins
%                 idx = bin_idx == bin;
%                 phase_angles = angle(wt(freq_idx, idx));
% 
%                 phase_angles = phase_angles(:);
% 
%                 if isempty(all_phases{bin, freq_idx})
%                     all_phases{bin, freq_idx} = phase_angles;
%                 else
%                     all_phases{bin, freq_idx} = [all_phases{bin, freq_idx}; phase_angles];
%                 end
%             end
%         end
%     end
% 
%     for freq_idx = 1:length(frequencies)
%         for bin = 1:num_bins
%             PPC{i}(bin, freq_idx) = ppc(all_phases{bin, freq_idx},2);
%         end
%     end
% end
% toc
% %%
% figure
% ppc2 = PPC{1}+PPC{2};
% ppc2=ppc2/2;
% h = pcolor(x_labels, f, ppc2');
% h.EdgeColor = 'none';
% shading interp;
% colorbar;
% xlabel('HPC Theta Phase (radians)');
% ylabel('PFC Frequency (Hz)');
% title('Pairwise Phase Consistency averaged in trials for High Amp Cycles');
%%
Bandpass = [6 12];  % Theta band
threshold_percentile = 75;
selected_cycles_HPC_upper = cell(2, 1);
selected_cycles_PFC_upper = cell(2, 1);
selected_cycles_LFP_theta = cell(2, 1);
selected_cycles_thetaPhase= cell(2, 1);

for i = 1:2
    HPC_LFP = rem_LFPHPC{i};
    PFC_LFP = rem_LFPPFC{i};

    LFPtheta = theta_filter(HPC_LFP, samplingRate);
    thetaPhaseHPC = angle(hilbert(LFPtheta));
%     LFPtheta = theta_filter(signalHPC, samplingRate);
%     LFPtheta = bandPassFSignal(signalHPC, samplingRate, Bandpass);
%     thetaPhase = angle(hilbert(LFPtheta));

    % Identify cycles based on phase wrapping from +pi to -pi
    phase_diff = diff([thetaPhaseHPC; thetaPhaseHPC(1)]);
    cycle_end_indices = find(phase_diff < -pi);
    cycle_start_indices = [cycle_end_indices(1:end-1)+1];
    cycle_end_indices = cycle_end_indices(2:end);
    
    % Extract cycles and their peak values
    peaks_all_cycles = arrayfun(@(start, finish) max(LFPtheta(start:finish)), cycle_start_indices, cycle_end_indices);
    
    % Determine power threshold for upper 75%
    power_threshold = prctile(peaks_all_cycles, threshold_percentile);
    
    % Select cycles with peaks in the upper 75%
    upper_cycles_indices = peaks_all_cycles > power_threshold;
    selected_cycles_HPC_upper{i} = arrayfun(@(idx) HPC_LFP(cycle_start_indices(idx):cycle_end_indices(idx)), find(upper_cycles_indices), 'UniformOutput', false);
    selected_cycles_LFP_theta{i} = arrayfun(@(idx) LFPtheta(cycle_start_indices(idx):cycle_end_indices(idx)), find(upper_cycles_indices), 'UniformOutput', false);
    selected_cycles_thetaPhase{i}=  arrayfun(@(idx) thetaPhaseHPC(cycle_start_indices(idx):cycle_end_indices(idx)), find(upper_cycles_indices), 'UniformOutput', false);
    selected_cycles_PFC_upper{i} = arrayfun(@(idx) PFC_LFP(cycle_start_indices(idx):cycle_end_indices(idx)), find(upper_cycles_indices), 'UniformOutput', false);

end
%%
for i=1:116
    figure
    signalHPC = selected_cycles_HPC_upper{2}{i};
    plot(signalHPC,LineWidth=2)
    pause
    close all
end
%%
num_bins = 20;
bin_edges = linspace(-pi, pi, num_bins + 1);
frequencies = 20:2:180;  
x_labels = (bin_edges(1:end-1) + bin_edges(2:end)) / 2; 

NormalizedPower_log10 = cell(2,1);
NormalizedPower_zphase = cell(2,1);
NormalizedPower_zfreq = cell(2,1);
NormalizedPower_total = cell(2,1);
Power_av = cell(2,1);

for i = 1:2
    Power = []; 
    for cycle_idx = 1:length(selected_cycles_HPC_upper{i})
        signalHPC = selected_cycles_HPC_upper{i}{cycle_idx};
        signalPFC = selected_cycles_PFC_upper{i}{cycle_idx};

        [wt, f] = cwavelet(signalPFC, frequencies, samplingRate);
         % Z-score normalization (Time)
        mean_power_time = mean(wt, 1);
        std_power_time = std(wt, 0, 1);
        NormalizedPower_ztime = (wt - mean_power_time) ./ std_power_time;
        wt = smoothCFSAbdel(NormalizedPower_ztime, flag_SMOOTH, frequency_smoothing_length, time_smoothing_length);

        %LFPtheta = theta_filter(signalHPC, samplingRate);
        %LFPtheta = bandPassFSignal(signalHPC, samplingRate, Bandpass);
        %thetaPhase = angle(hilbert(LFPtheta));
        thetaPhase = selected_cycles_thetaPhase{i}{cycle_idx};
        [~, ~, bin_idx] = histcounts(thetaPhase, bin_edges);

        Power_cycle = zeros(num_bins, length(f));  
        for freq_idx = 1:length(f)
            for bin = 1:num_bins
                idx = bin_idx == bin;
                Power_cycle(bin, freq_idx) = mean(abs(wt(freq_idx, idx)).^2);
            end
        end
        %Power = cat(3, Power, Power_cycle); % Concatenate the power matrix of this cycle to the rest
        Power(:,:,cycle_idx) = Power_cycle;
    end
    
    % Average the power over all cycles for this trial
    Power_avg = nanmean(Power, 3);
    Power_av{i} = Power_avg;
    
    % Log10 Normalization
    NormalizedPower_log10{i} = log10(Power_avg + eps);
    
    % Z-score normalization (Phase)
    mean_power_phase = mean(Power_avg, 1);
    std_power_phase = std(Power_avg, 0, 1);
    NormalizedPower_zphase{i} = (Power_avg - mean_power_phase) ./ std_power_phase;
    
    % Z-score normalization (Frequency)
    mean_power_freq = mean(Power_avg, 2);
    std_power_freq = std(Power_avg, 0, 2);
    NormalizedPower_zfreq{i} = (Power_avg - mean_power_freq) ./ std_power_freq;
end

NormalizedPower_zphase_avg = mean(cat(3, NormalizedPower_zphase{:}), 3);
%%
for cycle_idx = 1:length(selected_cycles_HPC_upper{i})
        signalHPC = selected_cycles_HPC_upper{i}{cycle_idx};
        signalPFC = selected_cycles_PFC_upper{i}{cycle_idx};

        [wt, f] = cwavelet(signalPFC, frequencies, samplingRate);
         % Z-score normalization (Time)
        mean_power_time = mean(wt, 1);
        std_power_time = std(wt, 0, 1);
        NormalizedPower_ztime = (wt - mean_power_time) ./ std_power_time;
        wt = smoothCFSAbdel(NormalizedPower_ztime, flag_SMOOTH, frequency_smoothing_length, time_smoothing_length);

        %LFPtheta = theta_filter(signalHPC, samplingRate);
        LFPtheta = bandPassFSignal(signalHPC, samplingRate, Bandpass);
        thetaPhase = angle(hilbert(LFPtheta));
        [~, ~, bin_idx] = histcounts(thetaPhase, bin_edges);
        bin_indx{cycle_idx}=bin_idx;
end
%%
figure;
%h = pcolor(x_labels, f, squeeze(Power(:,:,1))');
h = pcolor(x_labels, f, NormalizedPower_zphase_avg');
h.EdgeColor = 'none';
shading interp;
colorbar;
xlabel('HPC Theta Phase (radians)');
ylabel('PFC Frequency (Hz)');
title('Z-Scored Power Across Phases in High Amp Cycles');
%%
tic
PPC = cell(2,1);  

for i = 1:2
    PPC{i} = zeros(num_bins, length(frequencies)); 

    all_phases = cell(num_bins, length(frequencies));

    for cycle_idx = 1:length(selected_cycles_HPC_upper{i})
        signalHPC = selected_cycles_HPC_upper{i}{cycle_idx};
        %LFPtheta = theta_filter(signalHPC, samplingRate);
        LFPtheta = bandPassFSignal(signalHPC, samplingRate, Bandpass);
        thetaPhase = angle(hilbert(LFPtheta));

        signalPFC = selected_cycles_PFC_upper{i}{cycle_idx};
        [wt, f] = cwavelet(signalPFC, frequencies, samplingRate);
        % Z-score normalization (Time)
        mean_power_time = mean(wt, 1);
        std_power_time = std(wt, 0, 1);
        NormalizedPower_ztime = (wt - mean_power_time) ./ std_power_time;
        wt = smoothCFSAbdel(NormalizedPower_ztime, flag_SMOOTH, frequency_smoothing_length, time_smoothing_length);

        [~, ~, bin_idx] = histcounts(thetaPhase, bin_edges);

        for freq_idx = 1:length(frequencies)
            for bin = 1:num_bins
                idx = bin_idx == bin;
                phase_angles = angle(wt(freq_idx, idx));

                phase_angles = phase_angles(:);

                if isempty(all_phases{bin, freq_idx})
                    all_phases{bin, freq_idx} = phase_angles;
                else
                    all_phases{bin, freq_idx} = [all_phases{bin, freq_idx}; phase_angles];
                end
            end
        end
    end

    for freq_idx = 1:length(frequencies)
        for bin = 1:num_bins
            PPC{i}(bin, freq_idx) = ppc(all_phases{bin, freq_idx},2);
        end
    end
end
toc
figure
ppc2 = PPC{1}+PPC{2};
ppc2=ppc2/2;
h = pcolor(x_labels, f, ppc2');
h.EdgeColor = 'none';
shading interp;
colorbar;
xlabel('HPC Theta Phase (radians)');
ylabel('PFC Frequency (Hz)');
title('Pairwise Phase Consistency averaged in trials for High Amp Cycles');

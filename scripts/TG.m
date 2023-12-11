clear all; close all; clc
%% 1. reading the data from .mat file
HPClfp = load('./input_data/HPC_100_CH18_0.continuous.mat');
HPClfp = HPClfp.HPC;
PFClfp = load('./input_data/PFC_100_CH22_0.continuous.mat');
PFClfp = PFClfp.PFC;
samplingRate = 1250;
%% 2. Removing the artefacts based on amplitude and time window
[HPClfpCleaned,HPCnoisyInds] = removeArtefacts(HPClfp,samplingRate,[4 8], ...
    [1 0.1]);
[PFClfpCleaned,PFCnoisyInds] = removeArtefacts(PFClfp,samplingRate,[4 8], ...
    [1 0.1]);

%% 3.  
thetaband=[6 12];  %%%%6-12 Hz band passing
% samplingRate=1250; %%%%1000Hz;
LFPtheta = bandPassFSignal(HPClfpCleaned,samplingRate,thetaband);
thetaPhase=angle(hilbert(LFPtheta));

%% 4.

wCoh = wcoherGenzelVersion(s1,s2,scales,wname='amor',varargin);

% fs= 2500;
% %% downsampling the lfp signal
% HPClfp = decimate(HPClfp,2,'fir');
% PFClfp = decimate(PFClfp, 2, 'fir');
% fsNew = 1250;
% timeNewD = linspace(0,length(HPClfp)/fsNew, length(HPClfp));
% % for removing/detecting the artefacts:
% % Default values:
% thresholds = [4 8]; % we need to check for the thresholds (appropriate values)
% aroundArtefact = [1 0.1]; % we need to check for the time windows if we should change them
% % parameters
% threshold1 = thresholds(1); % in sigmas deviating from the mean
% aroundArtefact1 = aroundArtefact(1); % 2, Big and long artefacts
% threshold2 = thresholds(2); % for derivative of z-scored signal
% aroundArtefact2 = aroundArtefact(2); % 0.1 Very fast fluctuations (short time scale)
% %% performing the computations
% t = timeNewD';
% values = HPClfp;
% z = zscore(values);
% d = [diff(z);0];
% bad = false(size(values));
% %% first we detect the large global artefacts
% artefactInterval = t(findIntervalsA(abs(z) > threshold1));
% if numel(artefactInterval)==2
%     artefactInterval=artefactInterval(:)';
% end
% if ~isempty(artefactInterval)
%     artefactInterval = ConsolidateIntervalsFast([artefactInterval(:,1)-aroundArtefact1, artefactInterval(:,2)+aroundArtefact1]);
%     bad = InIntervals(t,artefactInterval);
% else
%     artefactInterval = zeros(0,2);
% end
% %Find noise using the derivative of the zscored signal (2)
% noisyInterval = t(findIntervalsA(abs(d)>threshold2));
% if numel(noisyInterval)==2
%     noisyInterval=noisyInterval(:)';
% end
% if ~isempty(noisyInterval)
%     noisyInterval = ConsolidateIntervalsFast([noisyInterval(:,1)-aroundArtefact2, noisyInterval(:,2)+aroundArtefact2]);
%     bad = bad | InIntervals(t,noisyInterval);
% else
%     noisyInterval = zeros(0,2);
% end
% % Substitute noisy signal with interpolated signal as if artefact did not exist
% values(bad) = interp1(t(~bad),values(~bad),t(bad,1));
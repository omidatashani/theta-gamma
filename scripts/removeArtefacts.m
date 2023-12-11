function [sigValues,artefactInds] = removeArtefacts(lfpSig,downSampleFreq,ampThresh, ...
    timeWinThresh)

%% downsampling the lfp signal
downSig = decimate(lfpSig,2,'fir');
downTime = linspace(0,length(downSig)/downSampleFreq, length(downSig));
% parameters
threshold1 = ampThresh(1); % in sigmas deviating from the mean
aroundArtefact1 = timeWinThresh(1); % 2, Big and long artefacts
threshold2 = ampThresh(2); % for derivative of z-scored signal
aroundArtefact2 = timeWinThresh(2); % 0.1 Very fast fluctuations (short time scale)

%% performing the computations
timeValues = downTime';
sigValues = downSig;
zSig = zscore(sigValues);
diffSig = [diff(zSig);0];
artefactInds = false(size(sigValues));

%% first we detect the large global artefacts
artefactInterval = timeValues(findIntervalsA(abs(zSig) > threshold1));
if numel(artefactInterval)==2
    artefactInterval=artefactInterval(:)';
end
if ~isempty(artefactInterval)
    artefactInterval = ConsolidateIntervalsFast([artefactInterval(:,1)- ...
        aroundArtefact1, artefactInterval(:,2)+aroundArtefact1]);
    artefactInds = InIntervals(timeValues,artefactInterval);
else
    artefactInterval = zeros(0,2);
end

%% Find noise using the derivative of the zscored signal (2)
noisyInterval = timeValues(findIntervalsA(abs(diffSig)>threshold2));
if numel(noisyInterval)==2
    noisyInterval=noisyInterval(:)';
end
if ~isempty(noisyInterval)
    noisyInterval = ConsolidateIntervalsFast([noisyInterval(:,1)- ...
        aroundArtefact2, noisyInterval(:,2)+aroundArtefact2]);
    artefactInds = artefactInds | InIntervals(timeValues,noisyInterval);
else
    noisyInterval = zeros(0,2);
end
% Substitute noisy signal with interpolated signal as if artefact did not exist
sigValues(artefactInds) = interp1(timeValues(~artefactInds), ...
    sigValues(~artefactInds) ,timeValues(artefactInds,1));
%% Uncomment this to run the function and get wcs
% File output_wcoherGenzel.mat contains output
% thetaband = [6 12];
% samplingRate = 1250;
% [filtHPC, IND] = removeArtefacts(HPC, 1250, [2 4], [2 0.1]);
% [filtPFC, pfcIND] = removeArtefacts(PFC, 1250, [2 4], [2 0.1]);
% LFPtheta = bandPassFSignal(filtHPC, samplingRate, thetaband);
% thetaPhase = angle(hilbert(LFPtheta));
% [sample, ThetaTS] = lu_wcsthetaextract(filtHPC, filtPFC, LFPtheta, thetaPhase, samplingRate, [1; 10801]);

%% Uncomment this to plot; load output mat file first
selection = 600;


figure
tiledlayout(4, 10, 'TileSpacing', 'compact', 'Padding', 'compact')
for i=1:10
nexttile(i)
plot(sample(:, selection+i));
title(string(selection+i))

nexttile(i+10)
plot(LFPOutput2{selection+i});
title("LFPOutput2")

nexttile(i+20)
plot(LFPOutput1{selection+i});
title("LFPOutput1")

nexttile(i+30)
reshaped = reshape(sample(:, selection+i), [185 20]);
imagesc('XData', transpose(PhaseBin), 'CData', real(reshaped))
imagesc(real(reshaped))
set(gca, 'YDir', 'normal')
end





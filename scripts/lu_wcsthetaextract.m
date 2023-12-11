function [sample,ThetaTS] = LU_WCSThetaExtract(LFP1,LFP2,LFPtheta,thetaphase,samprate,timerange,varargin)
%   get wavelet cross-spectrum in time&frequency domain representation
%   Extract Wavelet Cross-Spectrum for Each Single Theta Circle and Normalize Spectrum by Theta Phase.
%   Coded by Lu Zhang; math2437@hotmail.com
%   Last updated: (2nd Sept 2020)

%Output:
%sample is the vectorized cross-spectrum normalized in theta phase
%ThetaTS(i) is timeStamps(theta trough) of sample(:,i)


%LFP: raw LFPs
%LFPtheta: filtered LFPs in theta Band
%thetaphase: phase of filtered LFPs, usually calculated by hilbert transform;
%samprate: sampling rate;
%timerange: analysis period.[1 40 120;30 88 160] indicate 1-30s,
%40-88s and 120-160s for three analysis period in total;

% % Default setting
%WaveParam.ntw=100;      smoothing time-window 100
%WaveParam.nsw=3;        smoothing fre-window 3
%WaveParam.smoothW=10;   smoothing Power
%WaveParam.DownSample=5; DownSampling
%WaveParam.wname='morl'; Wavelet Name
%WaveParam.samprateD=samprateD;
%WaveParam.Freq=[1:150]; Frequeny band interested
% % Default setting


if nargin<7
% Default setting
WaveParam.Range=0.2;    %AlignedWindow Size 0.2s
WaveParam.ntw=100;      %smoothing time-window 100
WaveParam.nsw=3;        %smoothing fre-window 3
WaveParam.smoothW=20;   %smoothing Power
WaveParam.DownSample=5; %DownSampling
WaveParam.wname='morl'; %Wavelet Name
WaveParam.samprate=samprate;
WaveParam.Freq=[1:150]; %Frequeny band interested
WaveParam.Zscore=0; %default 0, raw power is calculated; ~=0, power is normalzied for each frequency respectively
%%%%%%%%%%%%%%%%%%%%%across time;
PhaseNum=20;
PhaseStep=2*pi/PhaseNum;
PhaseBin=-pi:PhaseStep:pi;
SavePath=[];
SaveShow=[];

% Default setting
elseif nargin==7
WaveParam=varargin{1};
WaveParam.samprate=samprate;
SavePath=[];
SaveShow=[];

PhaseNum=20;
PhaseStep=2*pi/PhaseNum;
PhaseBin=-pi:PhaseStep:pi;
elseif nargin==8
WaveParam=varargin{1};
WaveParam.samprate=samprate;
Param=varargin{2};
PhaseBin=Param.PhaseBin;
   if isfield(Param,'SavePath')
   SavePath=Param.SavePath;
   SaveShow=Param.SaveShow;
   else
   SavePath=["./"];
   SaveShow=[output];
   end

else
    
end

ThetaWCS={};
ThetaPhase={};
LFPOutput1={};
LFPOutput2={};

ThetaTimeStamps=[];


tic
Duration=round(diff(timerange));
NumThetaPeriod=length(Duration);

Istart=round(timerange(1,:)*samprate)+1;
Iend=round(timerange(2,:)*samprate)+1;
timerangeI=round(timerange*samprate)+1;

%%%%%%%%ii loops      multiple theta period in one recording file
for ii=1:length(Duration)

    %%%%%%%%Inlcuding 1s before the theta start and 1s after theta ends
    %%%%%%%%to avoid edge effects in time of Wavelet
    Istart(ii)=max(Istart(ii)-samprate,1);
    Iend(ii)=min(Iend(ii)+samprate,length(LFPtheta(:,1)));
    tempI=Istart(ii):Iend(ii);
    if isempty(tempI)
        continue
    end
    
    TExclude=[timerangeI(1,ii)-Istart(ii) Iend(ii)-timerangeI(2,ii)]/samprate;
    TimeInclude(:,ii)=[Istart(ii);Iend(ii)]/samprate;
    TInclude=([Istart(ii) Iend(ii)]-Istart(ii))/samprate;
    TInclude=[TInclude(1)+TExclude(1) TInclude(2)-TExclude(2)];
    
    LFPtemp1=LFP1(tempI);
    LFPtemp2=LFP2(tempI);

    LFPthetatemp=LFPtheta(tempI);
    LFPphasetemp=thetaphase(tempI);

    LFPtemp1=zscore(decimate(LFPtemp1,WaveParam.DownSample));
    LFPtemp2=zscore(decimate(LFPtemp2,WaveParam.DownSample));

    
    LFPthetatemp=downsample(LFPthetatemp,WaveParam.DownSample);
    LFPphasetemp=downsample(LFPphasetemp,WaveParam.DownSample);

    samprateD=samprate/WaveParam.DownSample;
    
F=WaveParam.Freq;
wname=WaveParam.wname;
ntw=WaveParam.ntw;
nsw=WaveParam.nsw;
bin_width=1/samprateD;
fc = centfrq(wname);
% scales=sort(fc./F.*samprateD); %current used version
MorletFourierFactor = 4*pi/(6+sqrt(2+6^2));
% % scales=sort(1./(F.*MorletFourierFactor)*samprateD);
% % scales=1./(F*MorletFourierFactor);
scales=1./(F*MorletFourierFactor);
% [F1, F2] = cwtfreqbounds(length(LFPtemp1), samprateD, Wavelet="amor");
% % S1Power = waveletLU(LFPtemp,scales,wname,'ntw',ntw,'nsw',nsw);
% % S1Power = waveletLU(LFPtemp,scales,wname,'ntw',ntw,'nsw',nsw);
%[WCOH,WCS,S1Power,S1Power2]=wcoherGenzelVersion(LFPtemp1,LFPtemp2,scales,'amor','ntw',ntw,'nsw',nsw);
[wCOH, WCS] = wcoherence(LFPtemp1,LFPtemp2,1250);
% [WCOH,WCS,S1Power,S1Power2]=wcoherLU(LFPtemp1,LFPtemp2,scales,wname,'ntw',ntw,'nsw',nsw);
% [WCOH,WCS,S1Power,S1Power2]=wcoher(LFPtemp1,LFPtemp2,scales,wname);%%%%'SAM' refers to sampling rate;
% [WCOH, WCS] = wcoherence(LFPtemp1, LFPtemp2, samprateD, "FrequencyLimits", [1 150], "NumScalesToSmooth", nsw);
% WCOH=abs(WCOH);

% % % if WaveParam.Zscore
% % %    S1Power=zscore(S1Power,0,2); %%%%%%%%%%%%%normalize over time
% % % % S1Power=zscore(S1Power,0,1);  %%%%%%%%%%%%%normalize over frequency
% % % end

Time=([1:length(LFPtemp1)]-1)/samprateD;
FPlot=F(end:-1:1);

temp=LFPphasetemp;
len_t=length(temp);
         thetamax_ts=[];
         for j=1:len_t
             if j==1&&temp(j)<-3.141
                 thetamax_ts=[thetamax_ts,j];
             elseif j==len_t&&temp(j)<-3.141
                 thetamax_ts=[thetamax_ts,j];
             elseif j~=1&&j~=len_t
                 if (temp(j-1)-temp(j))>6.0&&(temp(j+1)-temp(j))>0
                     thetamax_ts=[thetamax_ts,j-1+(pi-temp(j-1))/((pi-temp(j-1))+(pi+temp(j)))];

                 end
             else
                 
             end
         end
thetamax_ts=(thetamax_ts-1)/samprateD;
thetamax_ts(thetamax_ts>TInclude(2))=[];
thetamax_ts(thetamax_ts<TInclude(1))=[];

thetamax_ind=round(thetamax_ts*samprateD);
% thetamax_ts=round(thetamax_ts);
% % figure;
% % % plot(LFPthetatemp(:,2));hold on;plot(thetamax_ind,-pi,'r.')
% % plot(LFPthetatemp(:,1));hold on;plot(thetamax_ind,0,'r.')

% BackW=round(WaveParam.Range*samprateD);
% ForW=BackW;

if ii==1
    tempPhase=LFPphasetemp(:);
else
    tempPhase=[tempPhase;LFPphasetemp(:)];
end

thetamax_ts=thetamax_ts+TimeInclude(1,ii);
thetamax_ts=thetamax_ts(:);
for iii=1:(length(thetamax_ind)-1)
    s=round(max(thetamax_ind(iii),1));
    o=round(min(thetamax_ind(iii+1),len_t));
    ThetaWCS{end+1,1}=WCS(:,s:o)';
    ThetaPhase{end+1,1}=LFPphasetemp(s:o);
    LFPOutput1{end+1,1}=LFPtemp1(s:o);
    LFPOutput2{end+1,1}=LFPtemp2(s:o);

end
ThetaTimeStamps=[ThetaTimeStamps;thetamax_ts(1:(end-1))];


end
%%%%%%%%ii loops      multiple theta period in one recording file

   ThetaWCSN = aveY_discretizeX(ThetaWCS,ThetaPhase,PhaseBin);    
   for i=1:length(ThetaWCSN)
       ThetaWCSNew(:,:,i)=ThetaWCSN{i};
   end
   clear ThetaWCS;
   
   ThetaTS=ThetaTimeStamps; clear ThetaTimeStamps;
   
%     WaveParam.Freq=[20:2:120];
    Fplot=WaveParam.Freq(end:-1:1);
    Div=size(ThetaWCSNew);
    sample=reshape(ThetaWCSNew,Div(1)*Div(2),Div(3));
    Ivalid=find(sum(isnan(sample),1)>0);
    sample(:,Ivalid)=[];
    ThetaTS(Ivalid)=[];
    LFPOutput1(Ivalid)=[];
    LFPOutput2(Ivalid)=[];

    
% % % for i=1:size(sample,2)
% % %     sample(:,i)=sample(:,i)-min(sample(:,i));
% % %     sample(:,i)=sample(:,i)/sum(sample(:,i))+0.00001;
% % %     sample(:,i)=sample(:,i)/sum(sample(:,i));
% % % end

   
% if ~isempty(SavePath)
    %amp = amplitude threshold in removeArtefacts (s is sigma)
    %time = time threshold in removeArtefacts (101 = [1 0.1])
   %save('../output/output_wcoherGenzel_amp2s4s_time201.mat','sample','ThetaTS','LFPOutput1','LFPOutput2','WaveParam','PhaseBin','Fplot','Div','timerange', "-v7.3");
% end
toc



function yy = aveY_discretizeX(y,x,xbin,varargin)

    %%%%%%%y is m*n matrix, x is m*1 vector
    %%%%%%%divided x into xbin
    %%%%%%%averaged the y in each x xbin.
    
    if ~iscell(x)
        bins = discretize(x,xbin);    
        for ii=1:(length(xbin)-1)
            yy(ii,:)=mean(y(bins==ii,:),1);
        end
    elseif iscell(x)
        yyy={};
        for i=1:size(x,1) 
            for j=1:size(x,2)
                yyy{i,j}=aveY_discretizeX(y{i,j},x{i,j},xbin);
            end 
        end
        yy=yyy;
        clear yyy
    elseif isempty(x)
        yy=[];  
    else
        yy=[];  
      
    end

% if narvargin==4
%    yy=CellMat2Mat(cellY);
% end
% 
% function MatY=CellMat2Mat(cellY)
% 
% MatY=[];
% iM=1;
% if ~iscell(cellY)
%    MatY=cat(3,MatY);
% else
%     for i=1:size(cellY,1)
%         for j=1:size(cellY,2)
%             MatY=CellMat2Mat(cellY{i,j});
%         end 
%     end
% end

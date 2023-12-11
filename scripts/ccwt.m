function varargout = ccwt(x,desiredStep,varargin)

if iscell(x) && isnumeric(varargin{1}) && isnumeric(varargin{2}) && ...
        ~istimetable(x)
    scales = varargin{1};
    PSI = varargin{2};
    [varargout{1:nargout}] = wavelet.internal.cwt(x,scales,PSI);
    return;
end


if (numel(varargin)>=2) && isnumeric(varargin{1}) && ...
        (ischar(varargin{2}) || isStringScalar(varargin{2}) || ...
        iscell(varargin{2}) || isstruct(varargin{2}))
    
    if iscell(varargin{2})
        
        [varargout{1:nargout}] = wavelet.internal.cwt(x,varargin{:});
        return;
        
    elseif ischar(varargin{2}) || isStringScalar(varargin{2})
        if isStringScalar(varargin{2})
            varargin{2} = char(varargin{2});
        end
        
        WAVnoNum = deblank(regexprep(varargin{2},'\d',''));
        pat = '-';
        WAVnoPat = deblank(regexprep(WAVnoNum,pat,''));
        pat = '\.';
        tfdot = regexp(WAVnoPat,pat);
        if any(tfdot)
            WAVnoPat(tfdot) = [];
        end
        
        wavinfo = cellstr(wavemngr('tfsn'));
        
        if any(strcmp(WAVnoPat,wavinfo))
            
            [varargout{1:nargout}] = wavelet.internal.cwt(x,varargin{:});
            return;
        end
        
    elseif isstruct(varargin{2})
        
        [varargout{1:nargout}] = wavelet.internal.cwt(x,varargin{:});
        return;
    end
end



%Check nargin and nargout cwt
narginchk(1,14);
nargoutchk(0,5);
TTable = struct();
IsTimeTable = istimetable(x);
% datetimeLabelFlag = false;
if IsTimeTable
    % Get the RowTimes from the timetable and maintain
    OrigTimes = x.Properties.RowTimes;
    % For timetable record if Properties.RowTimes is a datetime array
    if isdatetime(OrigTimes)
        TTable.SampleTimes = x.Properties.RowTimes-x.Properties.RowTimes(1);
    else
        TTable.SampleTimes = x.Properties.RowTimes;
    end
    TTable.Format = TTable.SampleTimes.Format;
    % Convert the RowTimes from a duration or datetime array to
    % time vector.
    [TTable.times,TTable.units,TTable.convertFunc] = ...
        wavelet.internal.getDurationandUnits(TTable.SampleTimes);
    % Check the time vector for uniform sampling
    Tunif = wavelet.internal.isuniform(TTable.times);
    if ~Tunif
        error(message('Wavelet:cwt:NonuniformlySampled'));
    end
    % validate that the times are increasing
    validateattributes(TTable.times,{'double'},{'increasing'},'CWT','RowTimes');
    % Extract valid numeric data from time table
    % Return VariableNames as cell array
    x = wavelet.internal.CheckAndExtractTT(x);
    if size(x,2) ~= 1
        error(message('Wavelet:cwt:ttsinglevecvar','X'));
    end
end
% Allow both real and complex input, double and single precision
validateattributes(x,{'double','single'},{'vector','finite','nonempty'});

Norig = numel(x);
for kk = 1:numel(varargin)
    if isStringScalar(varargin{kk})
        varargin{kk} = char(varargin{kk});
    end
end
[tfFilterBank,fb] = isFilterBank(TTable,varargin{:});

if ~tfFilterBank
    fbcell = parseinputs(Norig,TTable,varargin{:});
    fb = constructFilterBank(fbcell);
end

if any(strcmpi(fb.Wavelet,{'Morse','amor'})) && nargout == 5
    [cfs,freq,coitmp,scalcfs] = fb.wt(x);
elseif nargout < 5
    [cfs,freq,coitmp] = fb.wt(x);
elseif nargout == 5 && strcmpi(fb.Wavelet,'bump')
    error(message('Wavelet:cwt:noscalcfs'));
end

% If sampling frequency is specified, dt = 1/fs
if isempty(fb.SamplingPeriod)
    % The default sampling interval is 1 for normalized frequency
    dt = 1/fb.SamplingFrequency;
    
else
    % Get the dt and Units from the duration object
    [dt,Units,dtfunc] = ...
        wavelet.internal.getDurationandUnits(fb.SamplingPeriod);
    
end

if IsTimeTable && isduration(OrigTimes)
    t = seconds(OrigTimes);
else
    t = 0:dt:Norig*dt-dt;
end

tfnormfreq = isNormalizedFrequency(fb);
[ga,be] = getGammaBeta(fb);

if nargout == 0 && ~isempty(fb.SamplingPeriod)
    plotscalogramperiod(cfs,freq,t,coitmp,Units,dtfunc)
elseif nargout == 0 && (~isempty(fb.SamplingFrequency) || tfnormfreq)
    [FourierFactor,sigmaT] = wavelet.internal.cwt.wavCFandSD(fb.Wavelet,...
        ga,be);
    % If the input is a timetable for nargout == 0, we will plot in hertz
    % using the RowTimes of the timetable
    if IsTimeTable
        plotscalogramfreq(FourierFactor,sigmaT,cfs,freq,t,tfnormfreq,OrigTimes)
    else
        plotscalogramfreq(FourierFactor,sigmaT,cfs,freq,t,tfnormfreq);
    end
    
end
if nargout > 0
    wt = cfs;
    p = freq;
    coi = coitmp;
    
    varargout{1} = wt;
    varargout{2} = p;
    varargout{3} = coi;
    varargout{4} = fb;
    if nargout == 5
        varargout{5} = scalcfs;
    end
end
%-------------------------------------------------------------------------
function fbcell = parseinputs(n,ttable,varargin)
isTimeTable = ~isempty(fieldnames(ttable));
% Set defaults.
params.fs = [];
params.dt = 1;
params.Ts = [];
params.sampinterval = false;
params.engunitflag = true;
params.WAV = 'morse';
params.ga = 3;
params.be = 20;
params.nv = 10;
params.no = [];
params.pad = true;
boundary = 'reflection';
params.normalizedfreq = true;
params.freqlimits= [];
params.periodlimits = [];

if isempty(varargin) && ~isTimeTable
    fbcell = {'SignalLength',n,'Wavelet','Morse','TimeBandwidth',60,...
        'Boundary','reflection'};
    return;
elseif isempty(varargin) && isTimeTable
    times = seconds(ttable.SampleTimes);
    fs = 1/mean(diff(times));
    fbcell = {'SignalLength',n,'Wavelet','Morse','TimeBandwidth',60,...
        'Boundary','reflection','SamplingFrequency',fs};
    return;
end

% Error out if there are any calendar duration objects
tfcalendarDuration = cellfun(@iscalendarduration,varargin);
if any(tfcalendarDuration)
    error(message('Wavelet:FunctionInput:CalendarDurationSupport'));
end

tfsampinterval = cellfun(@(x) isduration(x) && isscalar(x),varargin);

if any(tfsampinterval) && isTimeTable
    error(message('Wavelet:cwt:timetablesampperiod'));
end

if (any(tfsampinterval) && nnz(tfsampinterval) == 1)
    params.sampinterval = true;
    params.Ts = varargin{tfsampinterval>0};
    if (numel(params.Ts) ~= 1 ) || params.Ts <= 0 || isempty(params.Ts)
        error(message('Wavelet:FunctionInput:PositiveScalarDuration'));
    end
    
    params.engunitflag = false;
    params.normalizedfreq = false;
    varargin(tfsampinterval) = [];
end



%Look for Name-Value pairs
numvoices = find(strncmpi('voicesperoctave',varargin,1));

if any(numvoices) && numel(numvoices) == 1
    params.nv = varargin{numvoices+1};
    %validate the value is numeric, even and between 4 and 48
    %validateattributes(params.nv,{'numeric'},{'positive','scalar',...
    %    'even','>=',4,'<=',48},'cwt','VoicesPerOctave');
    varargin(numvoices:numvoices+1) = [];
    
end



morseparams = find(strncmpi('waveletparameters',varargin,1));
timeBandwidth = find(strncmpi('timebandwidth',varargin,1));
if any(morseparams) && any(timeBandwidth)
    error(message('Wavelet:cwt:paramsTB'));
end

if (any(morseparams) && (nnz(morseparams) == 1))
    morseParameter = varargin{morseparams+1};
    %validateattributes(morseParameter,{'numeric'},{'numel',2,...
    %    'positive','nonempty'},'cwt','WaveletParameters');
    
    params.ga = morseParameter(1);
    tb = morseParameter(2);
    % validateattributes(params.ga,{'numeric'},{'scalar',...
    %     'positive','>=',1},'cwt','gamma');
    % validateattributes(tb,{'numeric'},{'scalar',...
    %     '>',params.ga},'cwt','TimeBandwidth');
    % beta must be greater than 1
    params.be = tb/params.ga;
    if params.be>40
        error(message('Wavelet:cwt:TBupperbound'));
    end
    varargin(morseparams:morseparams+1) = [];
    
end

if (any(timeBandwidth) && (nnz(timeBandwidth) == 1))
    params.timebandwidth = varargin{timeBandwidth+1};
    % validateattributes(params.timebandwidth,{'numeric'},{'scalar',...
    %     'positive','>' 3,'<=',120},'cwt','TimeBandwidth');
    params.ga = 3;
    params.be = params.timebandwidth/params.ga;
    varargin(timeBandwidth:timeBandwidth+1) = [];
end

tffreqrange = find(strncmpi('frequencylimits',varargin,1));
tfperiodrange = find(strncmpi('periodlimits',varargin,1));

if any(tffreqrange) && any(tfperiodrange)
    error(message('Wavelet:cwt:freqperiodrange'));
elseif any(tffreqrange) && params.sampinterval
    error(message('Wavelet:cwt:freqrangewithts'));
elseif any(tfperiodrange) && ~params.sampinterval && ~isTimeTable
    error(message('Wavelet:cwt:periodwithsampfreq'));
end


if any(tffreqrange)
    params.freqlimits = varargin{tffreqrange+1};
    validateattributes(params.freqlimits,{'double'},{'numel',2,'finite'});
    varargin(tffreqrange:tffreqrange+1) = [];
end

if any(tfperiodrange)
    params.periodlimits = varargin{tfperiodrange+1};
    validateattributes(params.periodlimits,{'duration'},{'numel',2});
    varargin(tfperiodrange:tfperiodrange+1) = [];
end



extendsignal = find(strncmpi('extendsignal',varargin,1));

if any(extendsignal)
    params.pad = varargin{extendsignal+1};
    
    if ~isequal(params.pad,logical(params.pad))
        error(message('Wavelet:FunctionInput:Logical'));
    end
    if ~params.pad
        boundary = 'periodic';
    end
    varargin(extendsignal:extendsignal+1) = [];
    
end

% NumOctaves name-value pair. Not recommended
tfnumoctaves = find(strncmpi('numoctaves',varargin,1));
if any(tfnumoctaves) && (any(tfperiodrange) || any(tffreqrange))
    error(message('Wavelet:cwt:numoctavesfreqperiod'));
elseif any(tfnumoctaves) && ~(any(tfperiodrange) || any(tffreqrange))
    params.no = varargin{tfnumoctaves+1};
    varargin(tfnumoctaves:tfnumoctaves+1) = [];
    
end


% Only scalar left must be sampling frequency. We will validate that Fs
% is a scalar in the validateattributes call to catch unsupported vector
% inputs.

tfsampfreq = cellfun(@(x) (isvector(x) && isnumeric(x)),varargin);

if any(tfsampfreq) && isTimeTable
    error(message('Wavelet:cwt:timetablesampfreq'));
end

if (any(tfsampfreq) && (nnz(tfsampfreq) == 1) && isempty(params.Ts))
    params.fs = varargin{tfsampfreq};
    params.normalizedfreq = false;
    params.engunits = true;
    varargin(tfsampfreq) = [];
elseif any(tfsampfreq) && ~isempty(params.Ts)
    error(message('Wavelet:FunctionInput:SamplingIntervalOrDuration'));
elseif nnz(tfsampfreq)>1
    error(message('Wavelet:FunctionInput:Invalid_ScalNum'));
end



%Only char variable left must be wavelet
tfwav = cellfun(@ischar,varargin);
if (nnz(tfwav) == 1)
    params.WAV = varargin{tfwav>0};
    params.WAV = ...
        validatestring(params.WAV,{'morse','bump','amor'},'cwt','WAVNAME');
elseif nnz(tfwav)>1
    error(message('Wavelet:FunctionInput:InvalidChar'));
    
end

if any(strcmp(params.WAV,{'bump','amor'})) && (any(morseparams) || any(timeBandwidth))
    error(message('Wavelet:cwt:InvalidParamsWavelet'));
end

if ~isempty(params.no)
    params = validateandconvertno(n,params,desiredStep);
end

if isTimeTable
    % A TimeTable specified without a SamplingFrequency or
    % SamplingPeriod will be treated as hertz
    if ~isempty(params.periodlimits)
        DT = mean(diff(ttable.times));
        params.Ts = ttable.convertFunc(DT);
        params.Ts.Format = ttable.Format;
    else
        times = seconds(ttable.SampleTimes);
        params.fs = 1/mean(diff(times));
    end
end

if ~isempty(params.Ts)
    
    fbcell = {'SignalLength',n,'Wavelet',params.WAV,...
        'VoicesPerOctave',params.nv,'SamplingPeriod',params.Ts,...
        'PeriodLimits',params.periodlimits,'Boundary',boundary};
elseif ~isempty(params.fs)
    
    fbcell = {'SignalLength',n,'Wavelet',params.WAV,...
        'VoicesPerOctave',params.nv,'SamplingFrequency',params.fs,...
        'FrequencyLimits',params.freqlimits,'Boundary',boundary};
else
    fbcell = {'SignalLength',n,'Wavelet',params.WAV,...
        'VoicesPerOctave',params.nv,...
        'FrequencyLimits',params.freqlimits,'Boundary',boundary};
    
end

if any(timeBandwidth)
    fbcell{end+1} = 'TimeBandwidth';
    fbcell{end+1} = params.ga*params.be;
elseif any(morseparams)
    fbcell{end+1} = 'WaveletParameters';
    fbcell{end+1} = [params.ga tb];
end



%-----------------------------------------------------------------------
function plotscalogramperiod(wt,period,t,coitmp,Units,dtfunc)
wt = gather(wt);
coitmp = dtfunc(coitmp);
period = dtfunc(period);

antiAnalytic = (ndims(wt) == 3);

% Use magnitude limits in both the analytic and antianalytic parts
% for the colormap
cmin = min(abs(wt(:)));
cmax = max(abs(wt(:)));
if cmax <= cmin
    cmax = cmin+eps('single');
end

% Plot for CWT of real-valued signal
if ~antiAnalytic
    %%
    
    ylbl = [getString(message('Wavelet:wcoherence:Period')) ' (' Units ') '];
    xlbl = [getString(message('Wavelet:wcoherence:Time'))  ' (' Units ')'];
    titleString = getString(message('Wavelet:cwt:ScalogramTitle'));
    hf = gcf;
    clf;
    AX = axes('parent',hf);
    
    % The following axes hold must occur before any plotting or the default
    % interactivity is restored
    hold(AX,'on');
    hs = image('Parent',AX,...
        'XData',t,'YData',period,...
        'CData',abs(wt), ...
        'CDataMapping','scaled');
    
    
    AX.YLim = [min(period),max(period)];
    AX.XLim = [min(t) max(t)];
    AX.Layer = 'top';
    AX.YDir = 'normal';
    AX.YScale = 'log';
    
    title(AX, titleString);
    ylabel(AX, ylbl)
    xlabel(AX, xlbl)
    
    hcol = colorbar('peer', AX);
    hcol.Label.String = getString(message('Wavelet:cwt:Magnitude'));
    
    
    
    plot(AX,t,coitmp,'w--','linewidth',2);
    baselevel = max(AX.YLim);
    A1 = area(AX,t,coitmp,baselevel);
    A1.EdgeColor = 'none';
    A1.FaceColor = [0.8 0.8 0.8];
    alpha(A1,0.4);
    A1.PickableParts = 'none';
    
    hold(AX,'off');
    hf.NextPlot = 'replace';
    
elseif antiAnalytic
    
    ylbl = [getString(message('Wavelet:wcoherence:Period')) ' (' Units ') '];
    xlbl = [getString(message('Wavelet:wcoherence:Time'))  ' (' Units ')'];
    titleString = {getString(message('Wavelet:cwt:ScalogramTitle'));...
        getString(message('Wavelet:cwt:ScalogramTitlePos'))};
    titleString2 = getString(message('Wavelet:cwt:ScalogramTitleNeg'));
    
    sz = getSizeforPlot;
    hf = gcf;
    origunits = hf.Units;
    hf.Units = 'pixels';
    clf;
    hf.Visible = 'off';
    % Change Bottom and Height of figure to accommodate double title
    hf.Position = [hf.Position(1) hf.Position(2) sz(3)/3 sz(3)/2.5];
    movegui(hf,'center');
    hf.Visible = 'on';
    hf.Units = origunits;
    
    AX(1) = subplot(2,1,1);
    
    % The following axes hold must occur before any plotting or the default
    % interactivity is restored
    hold(AX(1),'on');
    newplot(AX(1));
    hs1 = image('Parent',AX(1),...
        'XData', t, 'YData', period,...
        'CData',abs(wt(:,:,1)), ...
        'CDataMapping','scaled');
    
    AX(1).YLim = [min(period),max(period)];
    AX(1).XLim = [min(t) max(t)];
    AX(1).CLim = [cmin cmax];
    AX(1).Layer = 'top';
    AX(1).YDir = 'normal';
    AX(1).YScale = 'log';
    title(AX(1), titleString);
    ylabel(AX(1), ylbl)
    
    hcol = colorbar('peer', AX(1));
    hcol.Label.String = getString(message('Wavelet:cwt:Magnitude'));
    
    
    
    plot(AX(1),t,coitmp,'w--','linewidth',2);
    baselevel = max(AX(1).YLim);
    A1 = area(AX(1),t,coitmp,baselevel);
    A1.EdgeColor = 'none';
    A1.FaceColor = [0.8 0.8 0.8];
    alpha(A1,0.4);
    A1.PickableParts = 'none';
    
    hold(AX(1),'off');
    
    
    AX(2) = subplot(2,1,2);
    
    % The following axes hold must occur before any plotting or the default
    % interactivity is restored
    hold(AX(2),'on');
    hs2 = image('Parent',AX(2),...
        'XData',t,'YData',period,...
        'CData',abs(wt(:,:,2)), ...
        'CDataMapping','scaled');
    
    
    AX(2).YLim = [min(period),max(period)];
    AX(2).XLim = [min(t) max(t)];
    AX(2).CLim = [cmin cmax];
    AX(2).Layer = 'top';
    AX(2).YDir = 'normal';
    AX(2).YScale = 'log';
    
    title(AX(2), titleString2);
    ylabel(AX(2), ylbl)
    xlabel(AX(2), xlbl)
    
    hcol2 = colorbar('peer', AX(2));
    hcol2.Label.String = getString(message('Wavelet:cwt:Magnitude'));
    
    plot(AX(2),t,coitmp,'w--','linewidth',2);
    baselevel = max(AX(2).YLim);
    A2 = area(AX(2),t,coitmp,baselevel);
    A2.EdgeColor = 'none';
    A2.FaceColor = [0.8 0.8 0.8];
    alpha(A2,0.4);
    A2.PickableParts = 'none';
    
    hold(AX(2),'off');
    hf.NextPlot = 'replace';
    
    
    % Link time axis but not frequency/period
    linkaxes(AX,'x');
    hf.NextPlot = 'replace';
    AX(1).Tag = 'wpos';
    AX(2).Tag = 'wneg';
    hs = [hs1 hs2];
    
end

% Install a custom data cursor

dataCursorBehaviorObj = hgbehaviorfactory('DataCursor');
set(dataCursorBehaviorObj,...
    'UpdateFcn',{@cwtPerCursorUpdateFunction, t, period});
for ii = 1:numel(hs)
    hgaddbehavior(hs(ii),dataCursorBehaviorObj);
    addlistener(hs(ii),'ObjectBeingDestroyed',@resetDataCursorManager);
end

% Override the Y tick labels once, then add callbacks to modify on zoom or
% pan.
for ii = 1:numel(AX)
    AX(ii).YTickLabel = AX(ii).YTick;
    addlistener(AX(ii),'YScale','PostSet', @(es, ed)cbCurrentMode(AX(ii)));
end

hzoom = zoom(hf);
hzoom.ActionPostCallback = @changeYTickLabels;
hpan = pan(hf);
hpan.ActionPreCallback = @prePan;
hpan.ActionPostCallback = @changeYTickLabels;

%-------------------------------------------------------------------------
function plotscalogramfreq(FourierFactor,sigmaT,wt,freq,t,normfreqflag,varargin)
[wt,freq] = gather(wt,freq);
cmin = min(abs(wt(:)));
cmax = max(abs(wt(:)));
if cmax <= cmin
    cmax = cmin+eps('single');
end
antiAnalytic = (ndims(wt) == 3);

if normfreqflag
    frequnitstrs = wavelet.internal.wgetfrequnitstrs;
    ylbl = frequnitstrs{1};
    coifactorfreq = 1;
    
elseif ~normfreqflag
    [freq,eng_exp,uf] = engunits(freq,'unicode');
    coifactorfreq = eng_exp;
    ylbl = wavelet.internal.wgetfreqlbl([uf 'Hz']);
    
end

if normfreqflag
    ut = 'Samples';
    dt = 1;
    coifactortime = 1;
else
    [t,eng_exp,ut] = engunits(t,'unicode','time');
    coifactortime = eng_exp;
    dt = mean(diff(t));
end


N = size(wt,2);

% We have to recompute the cone of influence for whatever scaling
% is done in time and frequency by engunits
% dt = dt*coifactortime;
FourierFactor = FourierFactor/coifactorfreq;
sigmaT = sigmaT*coifactortime;
coiScalar = FourierFactor/sigmaT;
samples = createCoiIndices(N);
coi = coiScalar*dt*samples;
invcoi = 1./coi;
invcoi(invcoi>max(freq)) = max(freq);
if ~isempty(varargin) && isdatetime(varargin{1})
    T = datenum(varargin{1});
    datetimeLabelFlag = true;
else
    T = t;
    datetimeLabelFlag = false;
    
end
if ~antiAnalytic
    %%
    if datetimeLabelFlag
        xlbl = getString(message('Wavelet:getfrequnitstrs:Date'));
    else
        xlbl = [getString(message('Wavelet:getfrequnitstrs:Time')) ' (' ut ')'];
    end
    
    hf = gcf;
    clf;
    AX = axes('parent',hf);
    
    % The following axes hold must occur before any plotting or the default
    % interactivity is restored
    hold(AX,'on');
    hs = image('Parent',AX,...
        'XData',T,'YData',freq,...
        'CData',abs(wt), ...
        'CDataMapping','scaled');
    AX.YLim = [min(freq),max(freq)];
    AX.XLim = [min(T) max(T)];
    AX.Layer = 'top';
    AX.YDir = 'normal';
    AX.YScale = 'log';
    
    title(AX, getString(message('Wavelet:cwt:ScalogramTitle')));
    ylabel(AX, ylbl)
    xlabel(AX, xlbl)
    
    hcol = colorbar('peer', AX);
    hcol.Label.String = getString(message('Wavelet:cwt:Magnitude'));
    
    
    
    plot(AX,T,invcoi,'w--','linewidth',2);
    if datetimeLabelFlag
        datetick('x','KeepLimits');
    end
    
    baselevel = min([min(AX.YLim) min(invcoi)]);
    A1 = area(AX,T,invcoi,baselevel);
    A1.EdgeColor = 'none';
    A1.FaceColor = [0.8 0.8 0.8];
    alpha(A1,0.4);
    A1.PickableParts = 'none';
    
    hold(AX,'off');
    hf.NextPlot = 'replace';
    
elseif antiAnalytic
    
    if datetimeLabelFlag
        xlbl = getString(message('Wavelet:getfrequnitstrs:Date'));
    else
        xlbl = [getString(message('Wavelet:getfrequnitstrs:Time')) ' (' ut ')'];
    end
    titleString = {getString(message('Wavelet:cwt:ScalogramTitle'));...
        getString(message('Wavelet:cwt:ScalogramTitlePos'))};
    titleString2 = getString(message('Wavelet:cwt:ScalogramTitleNeg'));
    
    sz = getSizeforPlot;
    hf = gcf;
    origunits = hf.Units;
    hf.Units = 'pixels';
    clf;
    hf.Visible = 'off';
    % Change Bottom and Height of figure to accommodate double title
    hf.Position = [hf.Position(1) hf.Position(2) sz(3)/3 sz(3)/2.5];
    movegui(hf,'center');
    hf.Visible = 'on';
    hf.Units = origunits;
    
    
    AX(1) = subplot(2,1,1);
    
    % The following axes hold must occur before any plotting or the default
    % interactivity is restored
    hold(AX(1),'on');
    hs1 = image('Parent',AX(1),...
        'XData',T,'YData',freq,...
        'CData',abs(wt(:,:,1)), ...
        'CDataMapping','scaled');
    
    AX(1).YLim = [min(freq),max(freq)];
    AX(1).XLim = [min(T) max(T)];
    AX(1).CLim = [cmin cmax];
    AX(1).Layer = 'top';
    AX(1).YDir = 'normal';
    AX(1).YScale = 'log';
    
    
    title(AX(1), titleString);
    ylabel(AX(1), ylbl)
    
    hcol = colorbar('peer', AX(1));
    hcol.Label.String = getString(message('Wavelet:cwt:Magnitude'));
    
    
    
    plot(AX(1),T,invcoi,'w--','linewidth',2);
    if datetimeLabelFlag
        datetick('x','KeepLimits');
    end
    baselevel = min([min(AX(1).YLim) min(invcoi)]);
    A1 = area(AX(1),T,invcoi,baselevel);
    A1.EdgeColor = 'none';
    A1.FaceColor = [0.8 0.8 0.8];
    alpha(A1,0.4);
    A1.PickableParts = 'none';
    hold(AX(1),'off');
    
    
    AX(2) = subplot(2,1,2);
    
    % The following axes hold must occur before any plotting or the default
    % interactivity is restored
    hold(AX(2),'on');
    hs2 = image('Parent',AX(2),...
        'XData', T,'YData', freq,...
        'CData',abs(wt(:,:,2)), ...
        'CDataMapping','scaled');
    
    AX(2).YLim = [min(freq),max(freq)];
    AX(2).XLim = [min(T) max(T)];
    AX(2).CLim = [cmin cmax];
    AX(2).Layer = 'top';
    AX(2).YDir = 'normal';
    AX(2).YScale = 'log';
    
    title(AX(2), titleString2);
    ylabel(AX(2), ylbl)
    xlabel(AX(2), xlbl)
    
    hcol2 = colorbar('peer', AX(2));
    hcol2.Label.String = getString(message('Wavelet:cwt:Magnitude'));
    
    
    
    plot(AX(2),T,invcoi,'w--','linewidth',2);
    if datetimeLabelFlag
        datetick('x','KeepLimits');
    end
    baselevel = min([min(AX(2).YLim) min(invcoi)]);
    A2 = area(AX(2),T,invcoi,baselevel);
    A2.EdgeColor = 'none';
    A2.FaceColor = [0.8 0.8 0.8];
    alpha(A2,0.4);
    A2.PickableParts = 'none';
    hold(AX(2),'off');
    
    
    % Link time axis but not frequency/period
    linkaxes(AX,'x');
    hf.NextPlot = 'replace';
    AX(1).Tag = 'wpos';
    AX(2).Tag = 'wneg';
    hs = [hs1 hs2];
    
end


% Install a custom data cursor

dataCursorBehaviorObj = hgbehaviorfactory('DataCursor');
if datetimeLabelFlag
    set(dataCursorBehaviorObj,...
        'UpdateFcn',{@cwtFreqCursorUpdateFunction, T, freq, datetimeLabelFlag});
else
    set(dataCursorBehaviorObj,...
        'UpdateFcn',{@cwtFreqCursorUpdateFunction, T, freq});
end
for ii = 1:numel(hs)
    hgaddbehavior(hs(ii),dataCursorBehaviorObj);
end
for ii = 1:numel(hs)
    addlistener(hs(ii),'ObjectBeingDestroyed',@resetDataCursorManager);
end
% Override the Y tick labels once, then add callbacks to modify on zoom or
% pan.

for ii = 1:numel(AX)
    AX(ii).YTickLabel = AX(ii).YTick;
    addlistener(AX(ii),'YScale','PostSet', @(es, ed)cbCurrentMode(AX(ii)));
end
hzoom = zoom(hf);
hzoom.ActionPostCallback = @changeYTickLabels;
hpan = pan(hf);
hpan.ActionPreCallback = @prePan;
hpan.ActionPostCallback = @changeYTickLabels;

%--------------------------------------------------------------------------
function [TF,FB] = isFilterBank(ttable,varargin)
isTimeTable = ~isempty(fieldnames(ttable));
TF = false;
FB = [];

fbstring = strcmpi(varargin,'filterbank');

if nnz(fbstring) == 1 && numel(varargin) == 2
    idx = find(fbstring);
    FB = varargin{idx+1};
    TF = strcmpi(class(FB),'cwtfilterbank');
elseif (nnz(fbstring) == 1 && numel(varargin) ~=2) || nnz(fbstring)>1
    error(message('Wavelet:cwt:filterbankpv'));
end

% For a single-variable timetable check that the sampling frequency or
% sampling period matches the definition in the filter bank
if isTimeTable && TF && ~isempty(FB.SamplingPeriod)
    checkttsamplerate(ttable,FB.SamplingPeriod);
elseif isTimeTable && TF && isempty(FB.SamplingPeriod)
    checkttsamplerate(ttable,FB.SamplingFrequency);
end






%--------------------------------------------------------------------------
function params = validateandconvertno(n,params,desiredStep)
fs = 1;
if ~isempty(params.fs)
    fs = params.fs;
end
% Frequencies
if isempty(params.Ts)
    if ~strcmpi(params.WAV,'morse')
        [minfreq,maxfreq] = cwtfreqbounds(n,fs,'Wavelet',params.WAV);
        maxno = floor(log2(maxfreq/minfreq));
    elseif strcmpi(params.WAV,'morse')
        [minfreq,maxfreq] = cwtfreqbounds(n,fs,...
            'Wavelet',params.WAV,...
            'WaveletParameters',[params.ga params.ga*params.be]);
        maxno = floor(log2(maxfreq/minfreq));
    end
    
    validateattributes(params.no,{'numeric'},{'scalar','finite','nonempty',...
        'integer','<=',maxno},'CWT','NumOctaves');
    minf = 2^(-params.no)*maxfreq;
        
        % New: adjust frequency limits based on desiredStep
        numSteps = floor((maxfreq - minf) / desiredStep);
        maxfreqAdjusted = minf + numSteps * desiredStep;  % Adjust maxfreq to match step size
        
        % FrequencyLimits input to filter bank
        params.freqlimits = [minf maxfreqAdjusted];
    
elseif ~isempty(params.Ts) && ~strcmpi(params.WAV,'morse')
    % Get anonymous function from Sampling Period
    [~,~,convertFunc] = wavelet.internal.getDurationandUnits(params.Ts);
    [maxperiod,minperiod] = cwtfreqbounds(n,params.Ts,'Wavelet',params.WAV);
    maxp = convertFunc(maxperiod);
    minp = convertFunc(minperiod);
    % Octaves are defined in terms of frequencies
    maxno = floor(log2(maxp/minp));
    validateattributes(params.no,{'numeric'},{'scalar','finite','nonempty',...
        'integer','<=',maxno},'CWT','NumOctaves');
    maxp = 2^(params.no)*minp;
    maxp = convertFunc(maxp);
    maxp.Format = params.Ts.Format;
    % PeriodLimits input to filter bank
    params.periodlimits = [minperiod maxp];
elseif ~isempty(params.Ts) && strcmpi(params.WAV,'morse')
    [~,~,convertFunc] = wavelet.internal.getDurationandUnits(params.Ts);
    [maxperiod,minperiod] = cwtfreqbounds(n,params.Ts,...
        'Wavelet',params.WAV,...
        'WaveletParameters',[params.ga params.ga*params.be]);
    maxp = convertFunc(maxperiod);
    minp = convertFunc(minperiod);
    % Octaves are defined in terms of frequencies
    maxno = floor(log2(maxp/minp));
    validateattributes(params.no,{'numeric'},{'scalar','finite','nonempty',...
        'integer','<=',maxno},'CWT','NumOctaves');
    maxp = 2^(params.no)*minp;
    maxp = convertFunc(maxp);
    maxp.Format = params.Ts.Format;
    params.periodlimits = [minperiod maxp];
end

%-------------------------------------------------------------------------
function fb = constructFilterBank(fbcell)
fb = cwtfilterbank(fbcell{:});
%-------------------------------------------------------------------------
function indices = createCoiIndices(N)
if isodd(N)  % is odd
    indices = 1:ceil(N/2);
    indices = [indices, fliplr(indices(1:end-1))];
else % is even
    indices = 1:N/2;
    indices = [indices, fliplr(indices)];
end
%--------------------------------------------------------------------------
function checkttsamplerate(ttable,fs)

if isduration(fs)
    times = ttable.times;
    DT = mean(diff(times));
    [sp,spunits] = wavelet.internal.getDurationandUnits(fs);
    if abs(DT-sp) > 10*eps(max(DT,sp)) || ...
            ~strcmp(spunits,ttable.units) || ...
            ~strcmp(fs.Format,ttable.Format)
        error(message('Wavelet:cwt:SampPeriodAgreement'));
    end
else
    
    % These are converted to seconds
    times = seconds(ttable.SampleTimes);
    sr = 1/mean(diff(times));
    if abs(sr-fs) > 10*eps(max(sr,fs))
        error(message('Wavelet:cwt:SampFreqAgreement'));
    end
    
    
    
end


function output_txt = cwtPerCursorUpdateFunction(~,event_obj, varargin)
%%
pos = event_obj.Position;
dataIndex = event_obj.DataIndex;

thisTime = pos(1);
thisPeriod = pos(2);
output_txt{1} = [getString(message('Wavelet:cwt:Time')) ': ', num2str(thisTime, 4)];
output_txt{2} = [getString(message('Wavelet:wcoherence:Period')) ': ', num2str(thisPeriod, 4)];
cVal = event_obj.Target.CData(dataIndex);
output_txt{3} = [getString(message('Wavelet:cwt:Magnitude')) ': ', num2str(cVal)];



function output_txt = cwtFreqCursorUpdateFunction(~,event_obj, varargin)
%%
Nvarargin = numel(varargin);
pos = event_obj.Position;
dataIndex = event_obj.DataIndex;
if Nvarargin < 3
    thisTime = pos(1);
else
    thisTime = datestr(pos(1));
end
thisFreq = pos(2);
if Nvarargin < 3
    output_txt{1} = [getString(message('Wavelet:cwt:Time')) ': ', num2str(thisTime, 4)];
else
    output_txt{1} = [getString(message('Wavelet:cwt:Time')) ': ', thisTime];
end
output_txt{2} = [getString(message('Wavelet:cwt:Freq')) ': ', num2str(thisFreq, 4)];
cVal = event_obj.Target.CData(dataIndex);
output_txt{3} = [getString(message('Wavelet:cwt:Magnitude')) ': ', num2str(cVal)];



function  changeYTickLabels(varargin)
% After zooming or panning, adjust the tick labels to prefer non-scientific
% format.
haxes = gca;
haxes.YTickLabel = haxes.YTick;

function prePan(varargin)
% During panning, use the matlab graphics labels.
haxes = gca;
haxes.YTickLabelMode = 'auto';


function cbCurrentMode(varargin)
varargin{1}.YScale = 'log';
warning(message('Wavelet:cwt:plotMustBeLogScale'));


function resetDataCursorManager(src,evt)
% Reset when CWT surface is destroyed.

% Respond to ObjectBeingDestroyed event
if strcmpi(evt.EventName,'ObjectBeingDestroyed')
    hf = ancestor(src, 'figure');
    % Reset all axes after first re-enabling interactions
    ax = findall(hf,'type','axes');
    
    arrayfun(@(x)cla(x,'reset'),ax);
    
else
    return;
end

%-------------------------------------------------------------------------
function sz = getSizeforPlot

% Compensates for dual monitor
monitorPositions = get(0,'MonitorPositions');
% Are there dual monitors
isDualMonitor = size(monitorPositions,1) > 1;

if isDualMonitor
    origins = monitorPositions(:,1:2);
    % Identify the primary monitor
    primaryMonitorIndex = find(origins(:,1)==1 & origins(:,2)==1,1);
    
    if isempty(primaryMonitorIndex)
        % pick the first monitor if this doesn't work.
        primaryMonitorIndex = 1;
    else
        primaryMonitorIndex = max(primaryMonitorIndex,1);
    end
    
    sz = monitorPositions(primaryMonitorIndex, :);
else
    sz = get(0, 'ScreenSize');
end
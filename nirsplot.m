function report_table = nirsplot(rawDotNirs,fcut,window,overlap,lambda_mask)
% raw: raw data in Homer format. Ex: raw = load('nirx_sample.nirs','-mat')
% fcut: 1x2 array [fmin fmax] representing the bandpass of the cardiac pulsation (default [0.5 2.5])
% window: length in seconds of the window to partition the signal with (defaut: 5)
% overlap: fraction overlap (0..0.99) between adjacent windows (default: 0, no overlap)
% lambda_mask: binary array mapping the selected two wavelength to correlate
% (default: [1 1 ...], the first two encountered, no matter how many there are)

close all
global nSources nDetectors nChannels secDur sampDur raw 
global qMats qltyThld mergewoiFlag SDMeasListAct

report_table = []; 
qltyThld = 0.9;
mergewoiFlag = true;
SDMeasListAct = [];
save_report_table = false;

if nargin < 4
    overlap = 0;
end

if nargin < 3
    window = 5;
end

if nargin < 2
    fcut = [0.5 2.5];
end

if nargin < 1
    nirsplotLoadFileGUI;
    return;
else
    raw = rawDotNirs;
end

% Creating 's' variable (stimuli matrix) from the information in StimDesign
if ~isfield(raw,'s')
    if isfield(raw,'StimDesign')
        nStim = length(raw.StimDesign);
        sTmp = zeros(size(raw.d,1),nStim);
        fs = 1/mean(diff(raw.t));
        for iStim = 1:nStim
            sTmp(floor(raw.StimDesign(iStim).onset * fs),iStim) = 1;
        end
        raw.s = sTmp;
        clear sTmp;
    else
        error('Stimuli information is not available.');
    end 
end

lambdas = unique(raw.SD.MeasList(:,4));
if nargin < 5
    lambda_mask = zeros(length(lambdas),1);
    lambda_mask(1)=1;
    lambda_mask(2)=1;
end

[sampDur,nChannels] = size(raw.d);
nChannels = nChannels/2;
nSources = size(raw.SD.SrcPos,1);
nDetectors = size(raw.SD.DetPos,1);
secDur = raw.t(end);


% Computation
[qMats] = qualityCompute(raw,fcut,window,overlap,lambda_mask);
% Create GUI
createGUI();
updateQPlots(qMats);

if save_report_table == true
   report_table = saveQuality(qMats);
end
% Wait for calls
end

%-------------------------------------------------------------------------
function createGUI()
global myAxes

% Main figure container
pos.main = [0.125 0.05 0.75 0.8]; % left, bottom, width, height
mainFig = figure('Units','normalized',...
    'Position',pos.main,'Visible','off',...
    'Name','NIRSPlot','NumberTitle','off','MenuBar','none');

% Axes
% SCI
myAxDim.width = 0.9;
myAxDim.height = (1/4)*0.75; % Four axes over the 80% of the figure
myAxDim.xSep = 0.05;
myAxDim.ySep = (1 - myAxDim.height*4) / 5;

pos.inspAx = [myAxDim.xSep,myAxDim.ySep+(0*(myAxDim.height+myAxDim.ySep)),...
    myAxDim.width,myAxDim.height];
myAxes.inspector = axes(mainFig,'Units','normalized',...
    'Position',pos.inspAx,...
    'Title','Inspector');
myAxes.inspector.XLabel.String = 'Time (s)';
myAxes.inspector.YLabel.String = 'Channel #';


pos.comboAx = [myAxDim.xSep,myAxDim.ySep+(1*(myAxDim.height+myAxDim.ySep)),...
    myAxDim.width,myAxDim.height];
myAxes.combo = axes(mainFig,'Units','normalized',...
    'Position',pos.comboAx,...
    'Title','Overall quality');
myAxes.combo.YLabel.String = 'Channel #';

pos.powerAx = [myAxDim.xSep,myAxDim.ySep+(2*(myAxDim.height+myAxDim.ySep)),...
    myAxDim.width,myAxDim.height];
myAxes.power = axes(mainFig,'Units','normalized',...
    'Position',pos.powerAx,...
    'Title','Power peak');
myAxes.power.YLabel.String = 'Channel #';

pos.sciAx = [myAxDim.xSep,myAxDim.ySep+(3*(myAxDim.height+myAxDim.ySep)),...
    myAxDim.width,myAxDim.height];
myAxes.sci = axes(mainFig,'Units','normalized',...
    'Position',pos.sciAx,...
    'Title','SCI');
myAxes.sci.YLabel.String = 'Channel #';

pos.inspectBtn = [myAxDim.xSep, (myAxDim.height+myAxDim.ySep),...
    0.08, myAxDim.ySep*0.8];
inspectBtn = uicontrol(mainFig,'Style', 'pushbutton', 'String', 'Inspect',...
    'FontSize',14,'FontWeight','bold','Units','normalized','Position', pos.inspectBtn,...
    'Callback', @inspectActive);

pos.helpBtn = [(pos.inspectBtn(1)+pos.inspectBtn(3))*1.15,...
    pos.inspectBtn(2),0.05,pos.inspectBtn(4)];
helpBtn = uicontrol(mainFig,'Style','pushbutton','String','?',...
    'FontSize',14,'FontWeight','bold','Units','normalized','Position',...
    pos.helpBtn,'Callback', @showHelp);

pos.bpBtn = [(pos.helpBtn(1)+pos.helpBtn(3))*1.15,...
    pos.inspectBtn(2),0.15,pos.inspectBtn(4)];
bpBtn = uicontrol(mainFig,'Style','pushbutton','String','Channel selection',...
    'FontSize',14,'FontWeight','bold','Units','normalized','Position',...
    pos.bpBtn,'Callback', @selectGoodChannels);

mainFig.Visible = 'on';
end

%-------------------------------------------------------------------------
function showHelp(source,event)
helpFig =  figure('Units','normalized',...
    'Visible','off','Position',[0.3,0.3,0.3,0.2],...
    'Name','NIRSPlot Help','NumberTitle','off','MenuBar','none');
helpStr = sprintf(['Controls\nLeft-click: Select a window and channel \n',...
    'Right-click: Select the complete channel \n',...
    'Up|Down key: move up/down one channel \n',...
    'Left|Right key: move forward/backward one window \n',...
    'ESC key: exit from the Inspector mode']);
pos.helpTxt = [0.1,0,0.9,0.9];
uicontrol(helpFig,'Style','text','String',helpStr,...
    'FontSize',14,'Units','normalized','Position',pos.helpTxt,...
    'HorizontalAlignment','left');
helpFig.Visible='on';
end

%-------------------------------------------------------------------------
function inspectActive(source,event)
global nChannels qMats

button = 0;
flagDispWindow = false;
pastChannel = 0;
pastWindow = 0;
while button ~= 27
    [iWindow,iChannel,button] = my_ginput(1);
    iWindow = round(iWindow);
    iChannel = round(iChannel);
    switch button
        case 1
            flagDispWindow = true;
        case 3
            xLimWindow = [1,(qMats.sampPerWindow*qMats.nWindows)];
            pastChannel = iChannel;
            flagDispWindow = false;
        case 28 % left-arrow key
            if flagDispWindow==true
                iWindow = pastWindow - 1;
            end
            iChannel = pastChannel;
        case 29 % right-arrow key
            if flagDispWindow==true
                iWindow = pastWindow + 1;
            end
            iChannel = pastChannel;
        case 30 % up-arrow key
            if flagDispWindow==true
                iChannel = pastChannel - 1;
                iWindow = pastWindow;
            else
                iChannel = pastChannel - 1;
            end
        case 31 % down-arrow key
            if flagDispWindow==true
                iChannel = pastChannel + 1;
                iWindow = pastWindow;
            else
                iChannel = pastChannel + 1;
            end
    end
    if flagDispWindow == true
        %xLimWindow = [(qMats.sampPerWindow*iWindow)+1,...
        %    qMats.sampPerWindow*(iWindow+1)];
        xLimWindow = [(qMats.sampPerWindow*(iWindow-1))+1,...
            (qMats.sampPerWindow*iWindow)];
    end
    if button ~=27
        if iChannel>0 && iChannel<=nChannels && iWindow>0 && iWindow<=qMats.nWindows
            updateIPlot(iChannel,xLimWindow,iWindow);
            pastWindow = iWindow;
            pastChannel = iChannel;
        end
    end
end
end

%-------------------------------------------------------------------------
function updateQPlots(qMats)
% UpdatePlot updates the quality plots with the 'qualityMats' input arg
% bad channels are ploted according to 'plot_bad' flag
global myAxes nChannels sclAlpha woi 

%Unpacking
sci_array = qMats.sciArray;
power_array = qMats.powerArray;
combo_array = qMats.comboArray;

sclAlpha = 0.2;

% Scalp Contact Index
imSci = imagesc(myAxes.sci,sci_array);
myAxes.sci.CLim = [0,1];
myAxes.sci.YLim =[1, nChannels];
%myAxes.sci.XLim =[1, size(sci_array,2)];
% myAxes.sci.XTick = round(linspace(1,qMats.nWindows,8));
% myAxes.sci.XTickLabel = num2mstr(round(linspace(1,qMats.nWindows*qMats.sampPerWindow,8)));
colormap(myAxes.sci,"gray");
colorbar(myAxes.sci,"eastoutside","Ticks",[0 1]);
myAxes.sci.YLabel.String = 'Channel #';
myAxes.sci.YLabel.FontWeight = 'bold';

% Power peak
imPower = imagesc(myAxes.power,power_array);
myAxes.power.CLim = [0, 0.12];
myAxes.power.YLim =[1, nChannels];
%myAxes.power.XLim =[1, size(power_array,2)];
colormap(myAxes.power,"gray");
colorbar(myAxes.power,"eastoutside","Ticks",[0 0.1 0.12]);
%        h.TickLabels ={'0','0.10*','0.12'};
myAxes.power.YLabel.String = 'Channel #';
myAxes.power.YLabel.FontWeight = 'bold';

% Combo quality
imCombo = imagesc(myAxes.combo,combo_array);
myAxes.combo.CLim = [0, 1];
myAxes.combo.YLim =[1, nChannels];
%myAxes.combo.XLim =[1, size(combo_array,2)];
colormap(myAxes.combo,[0 0 0;1 1 1]);
colorbar(myAxes.combo,"eastoutside","Ticks",[0 1]);
myAxes.combo.YLabel.String = 'Channel #';
myAxes.combo.YLabel.FontWeight = 'bold';

% For visual consistency among axes
myAxes.inspector.YLimMode = 'manual';
myAxes.inspector.YLabel.String = 'Channel #';
myAxes.inspector.XLabel.String = 'Time (s)';
myAxes.inspector.YLabel.FontWeight = 'bold';
myAxes.inspector.XLabel.FontWeight = 'bold';
colorbar(myAxes.inspector,'Visible','off');

woiMatrgb = zeros(nChannels,qMats.nWindows,3);
woiMatrgb(:,:,2) = woi.mat;
alphaMat = woi.mat * sclAlpha;

hold(myAxes.sci,'on');
imwoiMat = imagesc(myAxes.sci,woiMatrgb,'AlphaData',alphaMat);
hold(myAxes.power,'on');
imwoiMat = imagesc(myAxes.power,woiMatrgb,'AlphaData',alphaMat);
hold(myAxes.combo,'on');
imwoiMat = imagesc(myAxes.combo,woiMatrgb,'AlphaData',alphaMat);
%hold(myAxes.inspector,'on');
%imwoiMat = imagesc(myAxes.inspector,woiMatrgb,'AlphaData',woi.mat);


end

%-------------------------------------------------------------------------
function updateIPlot(iChannel,xLimWindow,iWindow)
global myAxes nChannels raw qMats woi sclAlpha
% qMats.cardiacData;  % Lambdas x time x channels
cla(myAxes.inspector);
plot(myAxes.inspector,raw.t(xLimWindow(1):xLimWindow(2)),...
    qMats.cardiacData(1,xLimWindow(1):xLimWindow(2),iChannel),'-b');
hold(myAxes.inspector,'on');
plot(myAxes.inspector,raw.t(xLimWindow(1):xLimWindow(2)),...
    qMats.cardiacData(2,xLimWindow(1):xLimWindow(2),iChannel),'-r');
myAxes.inspector.XLim= [raw.t(xLimWindow(1)),raw.t(xLimWindow(2))];
YLimStd = [min(qMats.cardiacData(:,xLimWindow(1):xLimWindow(2),iChannel),[],'all'),...
   max(qMats.cardiacData(:,xLimWindow(1):xLimWindow(2),iChannel),[],'all')]*1.05;
XLimStd = myAxes.inspector.XLim;

%yLab = myAxes.inspector.YLim(1)+((myAxes.inspector.YLim(2)-myAxes.inspector.YLim(1))*0.71);
%xLab = myAxes.inspector.XLim(1)+((myAxes.inspector.XLim(2)-myAxes.inspector.XLim(1))*0.93);
%text(myAxes.inspector,xLab,yLab,['Channel:',num2str(iChannel)]);
strLgnds = {'Lambda 1','Lambda 2'};

if (xLimWindow(2)-xLimWindow(1)+1) == (qMats.nWindows*qMats.sampPerWindow)
    xRect = 0.5; %Because of the offset at the begining of a window
    yRect = iChannel-0.5;
    wRect = qMats.nWindows;
    hRect = 1;
    poiMatrgb = zeros(nChannels,xLimWindow(2),3);
    poiMatrgb(:,:,2) = repmat(repelem(woi.mat(1,:),qMats.sampPerWindow),nChannels,1);
    alphaMat = poiMatrgb(:,:,2) * sclAlpha;   
    
    impoiMat = imagesc(myAxes.inspector,'XData',...
        [raw.t(xLimWindow(1)),raw.t(xLimWindow(2))],...
        'YData',YLimStd,'CData',poiMatrgb,'AlphaData',alphaMat);
   
    %onsets
    [~,c] = size(raw.s);  
    colorOnsets = [linspace(0.5,1,c)',...
        linspace(0,0.5,c)',linspace(1,0,c)'];
    for j=1:c
        %mapping from 0,1 to 0,25%ofPeakToPeak
        yOnset = (raw.s(xLimWindow(1):xLimWindow(2),j)*(YLimStd(2)-YLimStd(1))*0.25)-abs(YLimStd(1));
        plot(myAxes.inspector,raw.t(xLimWindow(1):xLimWindow(2)),...
            yOnset,'LineWidth',2,...
            'Color',colorOnsets(j,:));
        strLgnds(2+j) = {['Cond ',num2str(j)]};
    end

else
    xRect = iWindow-0.5;
    yRect = iChannel-0.5;
    wRect = 1;
    hRect = 1;
end
updateQPlots(qMats);
myAxes.inspector.YLim = YLimStd;
myAxes.inspector.XLim = XLimStd;
myAxes.inspector.YLabel.String = ['Channel ', num2str(iChannel)];

lgn = legend(myAxes.inspector,strLgnds,'Box','off');

rectangle(myAxes.combo,'Position',[xRect yRect wRect hRect],...
    'EdgeColor','m','FaceColor','none','Linewidth',2);
rectangle(myAxes.power,'Position',[xRect yRect wRect hRect],...
    'EdgeColor','m','FaceColor','none','Linewidth',2);
rectangle(myAxes.sci,'Position',[xRect yRect wRect hRect],...
    'EdgeColor','m','FaceColor','none','Linewidth',2);
end

%-------------------------------------------------------------------------
function [qualityMats] = qualityCompute(raw,fcut,window,overlap,lambda_mask)
global fs n_windows window_samples woi qltyThld

% Set the bandpass filter parameters
fs = 1/mean(diff(raw.t));
fcut_min = fcut(1);
fcut_max = fcut(2);
if fcut_max >= (fs)/2
    fcut_max = (fs)/2 - eps;
    warning(['The highpass cutoff has been reduced from ',...
        num2str(fcut(2)), ' Hz to ', num2str(fcut_max),...
        ' Hz to satisfy the Nyquist sampling criterion']);
end
[B1,A1]=butter(1,[fcut_min*(2/fs) fcut_max*(2/fs)]);

lambdas = unique(raw.SD.MeasList(:,4));
if nargin < 5
    lambda_mask = zeros(length(lambdas),1);
    lambda_mask(1)=1;
    lambda_mask(2)=1;
end
nirs_data = zeros(length(lambdas),size(raw.d,1),size(raw.d,2)/length(lambdas));
cardiac_data = zeros(length(lambdas),size(raw.d,1),size(raw.d,2)/length(lambdas)); % Lambdas x time x channels
for j = 1:length(lambdas)
    % Filter everything but the cardiac component
    idx = find(raw.SD.MeasList(:,4) == lambdas(j));
    nirs_data(j,:,:) = raw.d(:,idx);
    filtered_nirs_data=filtfilt(B1,A1,squeeze(nirs_data(j,:,:)));
    cardiac_data(j,:,:)=filtered_nirs_data./repmat(std(filtered_nirs_data,0,1),size(filtered_nirs_data,1),1); % Normalized heartbeat
end
overlap_samples = floor(window*fs*overlap);
window_samples = floor(window*fs);
n_windows = floor((size(cardiac_data,2)-overlap_samples)/(window_samples-overlap_samples));
cardiac_data = cardiac_data(find(lambda_mask),:,:);
sci_array = zeros(size(cardiac_data,3),n_windows);    % Number of optode is from the user's layout, not the machine
power_array = zeros(size(cardiac_data,3),n_windows);
%fpower_array = zeros(size(cardiac_data,3),n_windows);
cardiac_windows = zeros(length(lambdas),window_samples,size(raw.d,2)/length(lambdas),n_windows);
for j = 1:n_windows
    interval = (j-1)*window_samples-(j-1)*(overlap_samples)+1 : j*window_samples-(j-1)*(overlap_samples);
    cardiac_windows(:,:,:,j) = cardiac_data(:,interval,:);
end
parfor j = 1:n_windows
    cardiac_window = cardiac_windows(:,:,:,j);
    sci_array_channels = zeros(1,size(cardiac_window,3));
    power_array_channels = zeros(1,size(cardiac_window,3));
    fpower_array_channels = zeros(1,size(cardiac_window,3));
    for k = 1:size(cardiac_window,3)
        %cross-correlate the two wavelength signals - both should have cardiac pulsations
        similarity = xcorr(squeeze(cardiac_window(1,:,k)),squeeze(cardiac_window(2,:,k)),'unbiased');
        if any(abs(similarity)>eps)
            % this makes the SCI=1 at lag zero when x1=x2 AND makes the power estimate independent of signal length, amplitude and Fs
            similarity = length(squeeze(cardiac_window(1,:,k)))*similarity./sqrt(sum(abs(squeeze(cardiac_window(1,:,k))).^2)*sum(abs(squeeze(cardiac_window(2,:,k))).^2));
        end
        [pxx,f] = periodogram(similarity,hamming(length(similarity)),length(similarity),fs,'power');
        [pwrest,idx] = max(pxx(f<fcut_max)); % FIX Make it age-dependent
        sci=similarity(length(squeeze(cardiac_window(1,:,k))));
        power=pwrest;
        fpower=f(idx);
        sci_array_channels(k) = sci;
        power_array_channels(k) = power;
        fpower_array_channels(k) = fpower;
    end
    sci_array(:,j) = sci_array_channels;    % Adjust not based on machine
    power_array(:,j) = power_array_channels;
%    fpower_array(:,j) = fpower_array_channels;
end

% Summary analysis
[woi] = getWOI(window_samples,n_windows);
idxPoi = logical(woi.mat(1,:));

mean_sci_link  = mean(sci_array(:,idxPoi),2);
std_sci_link  = std(sci_array(:,idxPoi),0,2);
good_sci_link = sum(sci_array(:,idxPoi)>0.8,2)/size(sci_array(:,idxPoi),2);
mean_sci_window  = mean(sci_array(:,idxPoi),1);
std_sci_window  = std(sci_array(:,idxPoi),0,1);
good_sci_window = sum(sci_array(:,idxPoi)>0.8,1)/size(sci_array(:,idxPoi),1);

mean_power_link  = mean(power_array(:,idxPoi),2);
std_power_link  = std(power_array(:,idxPoi),0,2);
good_power_link = sum(power_array(:,idxPoi)>0.1,2)/size(power_array(:,idxPoi),2);
mean_power_window  = mean(power_array(:,idxPoi),1);
std_power_window  = std(power_array(:,idxPoi),0,1);
good_power_window = sum(power_array(:,idxPoi)>0.1,1)/size(power_array(:,idxPoi),1);

combo_array = (sci_array >= 0.8) & (power_array >= 0.10);
mean_combo_link  = mean(combo_array(:,idxPoi),2);
std_combo_link  = std(combo_array(:,idxPoi),0,2);
%good_combo_link = sum(combo_array(:,idxPoi)>0.1,2)/size(combo_array(:,idxPoi),2);
%good_combo_link = sum(combo_array(:,idxPoi),2)/size(combo_array(:,idxPoi),2);
good_combo_link  = mean(combo_array(:,idxPoi),2);
mean_combo_window  = mean(combo_array(:,idxPoi),1);
std_combo_window  = std(combo_array(:,idxPoi),0,1);
%good_combo_window = sum(combo_array(:,idxPoi)>0.1,1)/size(combo_array(:,idxPoi),1);
%good_combo_window = sum(combo_array(:,idxPoi),1)/size(combo_array(:,idxPoi),1);
good_combo_window = mean(combo_array(:,idxPoi),1);

% Detect artifacts and bad links
bad_links = find(good_combo_link<qltyThld);
bad_windows = find(good_combo_window<qltyThld);

% Packaging sci, peakpower and combo
qualityMats.sciArray    = sci_array;
qualityMats.powerArray  = power_array;
qualityMats.comboArray  = combo_array;
qualityMats.badLinks    = bad_links;
qualityMats.badWindows  = bad_windows;
qualityMats.sampPerWindow = window_samples;
qualityMats.fs = fs;
qualityMats.nWindows = n_windows;
qualityMats.cardiacData = cardiac_data;
qualityMats.gcl = good_combo_link;
qualityMats.gcw = good_combo_window; 

end




%-------------------------------------------------------------------------
function [woi] = getWOI(window_samples,n_windows)
global raw fs nChannels mergewoiFlag
% Assuming no overlapping conditions

% The maximum number of allowed samples is window_samples*n_windows to consider
% an integer number of windows, module(total_samples,n_windows) = 0
allowed_samp = window_samples*n_windows; 
poi = sum(raw.s(1:allowed_samp,:),2);
poi = poi(1:allowed_samp);
% Sometimes 's' variable encodes the stimuli durations by including consecutive
% values of 1. We are interested on the onsets, then we remove consecutive ones.
idxpoi = find(poi);
poi = zeros(size(poi));
poi(idxpoi(diff([0;idxpoi])>1)) = 1;
nOnsets = length(find(poi));
idxStim = find(poi);
interOnsetTimes = raw.t(idxStim(2:end)) - raw.t(idxStim(1:end-1));
medIntTime = median(interOnsetTimes);
iqrIntTime = iqr(interOnsetTimes);
%blckDurTime = (medIntTime/2) + (0.5*iqrIntTime);
blckDurTime = medIntTime + (0.5*iqrIntTime);
blckDurSamp = round(fs*blckDurTime);
blckDurWind = floor(blckDurSamp/window_samples);
woi = struct('mat',zeros(nChannels,n_windows),...
             'start',zeros(1,nOnsets),...
             'end',zeros(1,nOnsets));
woi_array = zeros(1,n_windows);
% Since we are relying on windows, we do not need the POIs variables instead 
% we need WOIs variables information
for i=1:nOnsets
    startPOI = idxStim(i)-blckDurSamp;
    if startPOI < 1
        startPOI = 1;
    end
    startWOI = floor(startPOI/window_samples);
    
    endPOI = idxStim(i)+blckDurSamp;
    if endPOI > allowed_samp
        endPOI = allowed_samp;
    end
    endWOI = ceil(endPOI/window_samples);
    poi(startPOI:endPOI) = 1;
    woi_array(startWOI:endWOI) = 1;
    woi.start(i) = startWOI;
    woi.end(i) = endWOI;
end

% See my comment about the preference of WOIs rather than of POIs, if POI
% information is needed, uncomment next two lines and return POIs variables
% poi = poi';
% poiMat_ = repmat(poi,nChannels,1); 

woiblank = 0;
idxInit = [];
%woitmp = woi_array;
woitmp = woi_array;

% If the gap's duration between two consecutives blocks of interest is less than the
% block's average duration, then those two consecutives blocks will merge. 
% This operation has effect visually (one bigger green block instead of 
% two green blocks with a small gap in between) and for quality results 
% since the windows inside of such a gap are now considered for quality computation.
for i =1:n_windows
    if woitmp(i) == 0
        if isempty(idxInit)
            idxInit = i;
        end
       woiblank = woiblank +1;
    else
        if ~isempty(idxInit)
            if (woiblank <= blckDurWind) 
               woitmp(idxInit:i) = 1;           
            end
            woiblank = 0;
            idxInit = [];
        end
    end
end
if mergewoiFlag == true
    woi_array = woitmp;
end
woi.mat = repmat(woi_array,nChannels,1);
end

%-------------------------------------------------------------------------
function dotNirsOutput = selectGoodChannels(source, events)
global qMats qltyThld raw
SDMeasListAct = ~qMats.badLinks;
bpGoodQuality(qMats, qltyThld, raw);
dotNirsOutput = 0;
end

%-------------------------------------------------------------------------
%!This function is not tested yet!
function report_table = saveQuality(qMats)

report_table = table({qMats.badLinks'}, {qMats.badWindows});
report_table.Properties.VariableNames = {'file_idx','Bad_Links','Bad_Windows'};

for i=1:size(report_table,1)
    a = report_table.Bad_Links{i};
    b = report_table.Bad_Windows{i};
    a1 = num2str(a);
    b1 = num2str(b);
    report_table.Bad_Links{i} = a1;
    report_table.Bad_Windows{i} = b1;
end


 writetable(report_table,'Quality_Report.xls');
end
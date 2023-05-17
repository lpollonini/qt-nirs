function quality_matrices = qtnirs(dotNirsFilePath,varargin)
% QT-NIRS is a Matlab-based tool for the quality assessment of fNIRS data. 
% QT-NIRS can quantify the quality of an fNIRS recording in two different ways, by using a GUI or through a function call.
% Graphically, the QTNIRS GUI allows the user to locate a working folder for processing and quantifying the .nirs files within the working folder. 
% Programmatically, the users also can retrieve a set of quality measures by calling qtnirs from a Matlab script. 
%
% Usage information
% Using QT-NIRS inside of a script allows the users to specify a set of 
% parameters for the quality assessment. The 'dotNirsFilePath' parameter 
% can be the path of a .nirs file or the path to a folder containing 
% several .nirs files or a struct containing:
%              d* 
%              t* 
%      CondNames
%            aux
%              s* 
%             SD* 
%     procResult
%       userdata 
%      procInput
% *required
% In addition to the .nirs path, the user 
% can specify a list of parameters in a pairwise mode including:
% 
% Parameter keyword Description
%  freqCut:   1x2 array [fmin fmax] representing the bandpass of the cardiac pulsation (default [0.5 2.5])
%  window :   length in seconds of the window (defaut: 5)
%  overlap:     fraction overlap (0..0.99) between adjacent windows (default: 0, no overlap)
%  qualityThreshold:   The required quality value (normalized; 0-1) of good-quality windows in every channel (default: 0.75)
%  conditionsMask:   A binary mask or keyword to indicate the conditions
%  for computing the periods of interest. Valid keywords include 'all' and 'resting' to consider all or none of the conditions (default: 'all').
%  lambdaMask:    A binary array mapping the selected two wavelength to compute the SCI (default: [1 1], the first two WLs)
%  dodFlag:     A flag indicating to work from DOD data (default: 0)
%  guiFlag:     A flag indicating whether to start or not the GUI.

%
% 
% An example of QT-NIRS usage is:
% 
% bpFmin = 0.5; bpFmax = 2.5;
% windowSec = 5;
% windowOverlap = 0;
% quality_threshold = 0.9;
% qualityMatrices = qtnirs([pwd,filesep,'tmpDotNirs.nirs'],...
%                 'freqCut',[bpFmin, bpFmax],...
%                 'window',windowSec,...
%                 'overlap',windowOverlap,....
%                 'qualityThreshold',quality_threshold,...j
%                 'conditionsMask','all',...
%                 'dodFlag',0,...
%                 'guiFlag',0);
% 
% The 'qualityMatrices' output variable is an structure that includes the set of fields:
% 
% sci_array: Matrix containing the SCI values (dimension: #timewindows X #channels)
% power_array: Matrix containing the PeakPower values (dimension: #timewindows X #channels)
% combo_array: Matrix containing the QualityMask values (dimension: #timewindows X #channels)
% combo_array_expanded: Matrix containing the QualityMask values for the advanced mode (#timewindows X #channels)
% bad_links: List of the channels (averaged across timewindows) below the 'qualityThreshold' value
% bad_windows: List of the timewindows (averaged across channels) below the 'qualityThreshold' value
% sampPerWindow: Samples within a time window
% fs: Sampling frequency
% n_windows: Number of time windows
% cardiac_data: Filtered version of the fNIRS data (dimension: #WLs X #samples X #channels)
% good_combo_link: List of the channels (averaged across the time windows) above the 'qualityThreshold' value
% good_combo_window: List of the time windows (averaged across the channels) above the 'qualityThreshold' value
% woi: Data structure cointaining the Windows of Interest information (epochs of interest throughout the recording)
% MeasListAct: Array mask of the channels achieving the required level of quality (length: #Channels X #WLs)

% REMOVE THE .git FOLDERS 
% qtpath = mfilename('fullpath');
% qtpath = qtpath(1:strfind(qtpath,[filesep 'qtnirs']));
% pathCell = regexp(path, pathsep, 'split');
% if ispc  % Windows is not case-sensitive
%   onPath = any(strcmpi(qtpath, pathCell));
% else
%   onPath = any(strcmp(qtpath, pathCell));
% end
% %cd(qtpath(1:end-7))
% % Add that folder plus all subfolders to the path.
% %addpath(pwd);
% if ~onPath
%     addpath(genpath(qtpath));
% end
% 
% if nargin < 1
%     qtnirsLoadFileGUI(pwd);
%     return;
% end
if exist('dotNirsFilePath','var')==0
    dotNirsFilePath = pwd;
end
if ischar(dotNirsFilePath)
    if isfile(dotNirsFilePath)
        [filepath,name,ext] = fileparts(dotNirsFilePath);
        switch ext
            case '.nirs'
                rawNirs = load(dotNirsFilePath,'-mat');
            case '.snirf'   
                rawSnirf = SnirfClass(dotNirsFilePath);
                rawNirs.d = rawSnirf.Get_d;
                %rawNirs.s = rawSnirf.Get_s;
                %rawNirs.t = rawSnirf.Get_t;
                rawNirs.t = rawSnirf.data.time;
                rawNirs.s = rawSnirf.GetStims(rawNirs.t);
                rawNirs.SD = rawSnirf.Get_SD;
                rawNirs.aux = rawSnirf.GetAux;
            otherwise
                 error('The input file should be a .nirs file format');
        end
%         if strcmpi(ext,'.nirs')
%             rawNirs = load(dotNirsFilePath,'-mat');
%         else
%             error('The input file should be a .nirs file format');
%             %We should check that the required variables are in the file
%         end
    elseif isfolder(dotNirsFilePath)
        disp(['The input data is a folder. ',...
            'All .nirs files inside ',dotNirsFilePath,' will be evaluated.']);
        qtnirsLoadFileGUI(dotNirsFilePath);
        return;
    else
        error('The input path does not exist.');
    end
elseif isstruct(dotNirsFilePath)
    rawNirs = dotNirsFilePath;
    flagValidStruct = (isfield(rawNirs,'d') && ...
        isfield(rawNirs,'t') && ...
        isfield(rawNirs,'SD') && ...
        (isfield(rawNirs,'StimDesign') || isfield(rawNirs,'s')));

    if flagValidStruct == true
        filepath = pwd;
        name = 'QTNIRSAnalized';
        ext = '.nirs';
    else
        error(['The input data does not have the required fields (t, SD, s).']);
    end
    
elseif isa(dotNirsFilePath,'nirs.core.Data')
    rawNirs = struct;
    rawNirs.SD = nirs.util.probe2sd( dotNirsFilePath.probe );
    rawNirs.d = dotNirsFilePath.data;
    rawNirs.t = dotNirsFilePath.time;
    rawNirs.s = false(size(rawNirs.t));
    for icond = 1:length(dotNirsFilePath.stimulus.keys)
        stim = dotNirsFilePath.stimulus( dotNirsFilePath.stimulus.keys{icond} );
        rawNirs.s = rawNirs.s | stim.getStimVector( rawNirs.t );
    end
    
    rawNirs.s=1*rawNirs.s;
    rawNirs.aux=rawNirs.s;    
    filepath = pwd;
    name = 'QTNIRSAnalized';
    ext = '.nirs';   
end

propertyArgIn = varargin;
while length(propertyArgIn) >= 2
    prop = propertyArgIn{1};
    val = propertyArgIn{2};
    propertyArgIn = propertyArgIn(3:end);
    switch prop
        case 'freqCut'
            if isfloat(val) && length(val)==2
                fcut_ = [min(val) max(val)];
            end
            
        case 'window'
            if length(val)==1
                window_ = ceil(val);
            end
            
        case 'overlap'
            if isfloat(val) && val >= 0 && val <= 1
                overlap_ = val;
            end
            
        case 'qualityThreshold'
            if isfloat(val) && val >= 0 && val <= 1
                q_threshold = val;
            end
        case 'sciThreshold'
            if isfloat(val) && val >= 0 && val <= 1
                sci_threshold = val;
            end
        case 'pspThreshold'
            if isfloat(val) && val >= 0 && val <= 1
                psp_threshold = val;
            end
        case 'conditionsMask'          
            if ischar(val)
                switch val
                    case 'all'
                        cond_mask = ones(1,size(rawNirs.s,2));
                    case 'resting'
                        cond_mask = val;
                    otherwise
                        warning(['Value ',val,' in "cond_mask" parameter not valid.',...
                            'Using the complete scan time.']);
                        cond_mask = 'resting';
                end
            elseif isnumeric(val)
                if ~(any(val>1 | val<0))
                    cond_mask = val;
                else
                    warning(['Value ',num2str(val),' not valid, use only binary values (0/1).',...
                        'Using the complete scan time.']);
                    cond_mask = 'resting';                    
                end
            end
            
        case 'lambdaMask'
            errFlag = false;
            if ~(any(val>1 | val<0))                
                if length(val) ~= length(rawNirs.SD.Lambda)
                    warning(['The number of elements in the "lambdaMask" mask',...
                    ' does not match with the number of WLs in the scan.',...
                    ' Using the first two WLs.']);
                    errFlag = true;
                end
                if sum(val)~=2
                   warning(['Incorrect number of "1" in the "lambdaMask" mask.',...
                       ' Set ONLY two vaules to "1". Using the first two WLs.']);
                   errFlag = true;
                end
                if ~errFlag 
                    lambda_mask_ = logical(val);           
                    lambdas_ = find(lambda_mask_);
                end
            else
                 warning(['Elements in the "lambdaMask" mask',...
                    ' must be "1" and/or "0". Using the first two WLs.']);
            end
            
        case 'dodFlag'
            if val == 1
                if isfield(rawNirs,'procResult') && ~isempty(rawNirs.procResult.dod)
                    dodFlag_ = 1;
                else
                    dodFlag_ = 0;
                    warning('OD data is not available, I will use the raw data.');                
                end
            else
                dodFlag_ = 0;
            end                

    case 'guiFlag'
            if val == 1
                guiFlag_ = 1;
            else
                guiFlag_ = 0;
            end
    end
end
frequency_samp = 1/mean(diff(rawNirs.t));
% Creating 's' variable (stimuli matrix) from the information in StimDesign
if ~isfield(rawNirs,'s')
    if isfield(rawNirs,'StimDesign')
        nStim = length(rawNirs.StimDesign);
        sTmp = zeros(size(rawNirs.d,1),nStim);
        for iStim = 1:nStim
            for iOnset=1:length(rawNirs.StimDesign(iStim).onset)
                onsetTmp = floor(rawNirs.StimDesign(iStim).onset(iOnset) * frequency_samp);
                durTmp = floor(rawNirs.StimDesign(iStim).dur(iOnset)* frequency_samp);
                %sTmp(floor(rawNirs.StimDesign(iStim).onset(iOnset) * frequency_samp),iStim) = 1;
                sTmp(onsetTmp:(onsetTmp+durTmp),iStim) = 1;
            end
        end
        rawNirs.s = sTmp;
        clear sTmp;
    else
        error('Stimuli information is not available.');
    end
end

if ~exist('fcut_','var')
    fcut_ = [0.5 2.5];
end

if ~exist('window_','var')
    window_ = 5;
end
if ~exist('overlap_','var')
    overlap_ = 0;
end
if ~exist('q_threshold','var')
    q_threshold = 0.75;
end
if ~exist('sci_threshold','var')
    sci_threshold = 0.8;
end
if ~exist('psp_threshold','var')
    psp_threshold = 0.1;
end
if ~exist('cond_mask','var') || strcmp(cond_mask,'all')
    cond_mask = ones(1,size(rawNirs.s,2));
end
if ~exist('lambda_mask_','var')
    lambdas_ = unique(rawNirs.SD.MeasList(:,4))';
    lambdas_ = lambdas_(1:2);
    lambda_mask_ = true(1,length(lambdas_));
    if length(lambda_mask_) ~= length(rawNirs.SD.Lambda)
        for ii=3:length(rawNirs.SD.Lambda)
            lambda_mask_(end+1) = 0;
        end
    end
end
if ~exist('dodFlag_','var')
    dodFlag_ = -1;
end
if ~exist('guiFlag_','var')
    guiFlag_ = 1;
end

%------ Sorting for nirstoolbox compatibility ------
varNames = {'source','detector','dummy','type'};
MeasList_table = table(rawNirs.SD.MeasList(:,1),...
    rawNirs.SD.MeasList(:,2),...
    rawNirs.SD.MeasList(:,3),...
    rawNirs.SD.MeasList(:,4),...
    'VariableNames',varNames);

colsToSortBy = {'source', 'detector', 'type'};
[MeasList_table, idxML] = sortrows(MeasList_table, colsToSortBy);
rawNirs.SD.MeasList = table2array(MeasList_table);
rawNirs.d = rawNirs.d(:,idxML);
if dodFlag_ == 1
    rawNirs.procResult.dod = rawNirs.procResult.dod(:,idxML);
end
%---------------------------------------------------

nirsplot_parameters.dotNirsPath = filepath;
nirsplot_parameters.dotNirsFile = name;
nirsplot_parameters.dotNirsExt = ext;
nirsplot_parameters.fcut = fcut_;
nirsplot_parameters.window = window_;
nirsplot_parameters.overlap = overlap_;
nirsplot_parameters.lambda_mask = lambda_mask_;
nirsplot_parameters.lambdas = lambdas_;
nirsplot_parameters.dodFlag = dodFlag_;
nirsplot_parameters.mergewoi_flag = false;
nirsplot_parameters.quality_threshold = q_threshold;
nirsplot_parameters.sci_threshold = sci_threshold;
nirsplot_parameters.psp_threshold = psp_threshold;
% Bug found by Benjamin Zinszer (miscalculation in the number of channels)
nirsplot_parameters.n_channels = size(rawNirs.d,2)/length(rawNirs.SD.Lambda);
nirsplot_parameters.n_sources = size(rawNirs.SD.SrcPos,1);
nirsplot_parameters.n_detectors = size(rawNirs.SD.DetPos,1);
nirsplot_parameters.s = rawNirs.s;
nirsplot_parameters.t = rawNirs.t;


nirsplot_parameters.fs = frequency_samp;
nirsplot_parameters.mergewoiFlag = true;
nirsplot_parameters.cond_mask = cond_mask;
nirsplot_parameters.save_report_table = false;
nirsplot_parameters.sclAlpha = 0.65;
nirsplot_parameters.rectangle_line_width = 1.2;
nirsplot_parameters.guiFlag = guiFlag_;

% Call the GUI for parameter inputs
S=dbstack;
if guiFlag_ == 1
    if length(S)== 1
        qtnirsLoadFileGUI(nirsplot_parameters)
    end
    
    % Create GUI
    prev_window = findobj('type','figure');
    if ~isempty(prev_window)
        close(prev_window);
    end
    [main_fig_axes,main_fig] = createGUI();
else
    
    main_fig = figure('Units','normalized',...
        'Visible','off','Name','QT-NIRS',...
        'NumberTitle','off','MenuBar','none');
    main_fig_axes = [];
    
end
nirsplot_parameters.main_fig_axes = main_fig_axes;
setappdata(main_fig,'nirsplot_parameters',nirsplot_parameters);

switch ext
    case '.nirs'
        setappdata(main_fig,'rawNirs',rawNirs);
    case '.snirf'
        %setappdata(main_fig,'rawNirs',rawSnirf);
        setappdata(main_fig,'rawNirs',rawNirs);
end

% Computation
[quality_matrices] = qualityCompute(main_fig);
nirsplot_parameters.quality_matrices = quality_matrices;

setappdata(main_fig,'nirsplot_parameters',nirsplot_parameters);

if guiFlag_ == 1
    main_fig.Visible = 'on';
    hideMACheckB = findobj('Tag','hideMACheckB');
    updateQPlots(hideMACheckB,[]);
    %rawNT = nirs.io.loadDotNirs(dotNirsFilePath,true);
    %rawNT.draw(1:(nirsplot_parameters.n_channels * 2),[],main_fig_axes.inspector);
else
    close(main_fig);
end

report_table = [];
if nirsplot_parameters.save_report_table == true
    report_table = saveQuality(quality_matrices);
end
% Wait for calls


%% -------------------------------------------------------------------------
    function [main_fig_axes,main_fig] = createGUI()
        
        % Main figure container
        pos.main = [0.20 0.05 0.75 0.85]; % left, bottom, width, height
        main_fig = figure('Units','normalized',...
            'Position',pos.main,'Visible','off',...
            'Name','QT-NIRS','NumberTitle','off','MenuBar','none','Toolbar','figure');
        
        % Axes
        % SCI
        myAxDim.width = 0.9;
        myAxDim.height = (1/4)*0.75; % Four axes over the 80% of the figure
        myAxDim.xSep = 0.04;
        myAxDim.ySep = (1 - myAxDim.height*4) / 5;
        
        pos.inspAx = [myAxDim.xSep,myAxDim.ySep+(0*(myAxDim.height+myAxDim.ySep)),...
            myAxDim.width,myAxDim.height];
        main_fig_axes.inspector = axes(main_fig,'Units','normalized',...
            'Position',pos.inspAx,...
            'Title','Inspector');
        main_fig_axes.inspector.XLabel.String = 'Time (s)';
        main_fig_axes.inspector.YLabel.String = 'Channel #';
        
        pos.comboAx = [myAxDim.xSep,myAxDim.ySep+(1.05*(myAxDim.height+myAxDim.ySep)),...
            myAxDim.width,myAxDim.height];
        main_fig_axes.combo = axes(main_fig,'Units','normalized',...
            'Position',pos.comboAx,...
            'Title','Overall quality');
        main_fig_axes.combo.YLabel.String = 'Channel #';
        main_fig_axes.combo.YLabel.FontWeight = 'bold';
        colorbar(main_fig_axes.combo,'Visible','off','Tag','colorb_combo')
        
        
        pos.powerAx = [myAxDim.xSep,myAxDim.ySep+(2*(myAxDim.height+myAxDim.ySep)),...
            myAxDim.width,myAxDim.height];
        main_fig_axes.power = axes(main_fig,'Units','normalized',...
            'Position',pos.powerAx,...
            'Title','Power peak');
        main_fig_axes.power.YLabel.String = 'Channel #';
        main_fig_axes.power.YLabel.FontWeight = 'bold';
        colorbar(main_fig_axes.power,'Visible','off','Tag','colorb_power')
        
        
        pos.sciAx = [myAxDim.xSep,myAxDim.ySep+(3*(myAxDim.height+myAxDim.ySep)),...
            myAxDim.width,myAxDim.height];
        main_fig_axes.sci = axes(main_fig,'Units','normalized',...
            'Position',pos.sciAx,...
            'Title','SCI');
        main_fig_axes.sci.YLabel.String = 'Channel #';
        main_fig_axes.sci.YLabel.FontWeight = 'bold';
        colorbar(main_fig_axes.sci,'Visible','on','Tag','colorb_sci')
        
        pos.inspectBtn = [myAxDim.xSep, (myAxDim.height+myAxDim.ySep)*1.025,...
            0.08, myAxDim.ySep*0.7];
        uicontrol(main_fig,'Style', 'pushbutton', 'String', 'Inspect',...
            'FontSize',12,'FontWeight','bold','Units','normalized','Position', pos.inspectBtn,...
            'Callback', @inspectActive,'Tag','inspectBtn');
        
        pos.helpBtn = [pos.inspectBtn(1)+pos.inspectBtn(3)+myAxDim.xSep,...
            pos.inspectBtn(2),0.05,pos.inspectBtn(4)];
%         uicontrol(main_fig,'Style','pushbutton','String','?',...
%             'FontSize',12,'FontWeight','bold','Units','normalized','Position',...
%             pos.helpBtn,'Callback', @showHelp,'Tag','helpBtn');
        
        pos.chSelBtn = [pos.inspectBtn(1)+pos.inspectBtn(3)+myAxDim.xSep,...
            pos.inspectBtn(2),0.15,pos.inspectBtn(4)];
        uicontrol(main_fig,'Style','pushbutton','String','Channel selection',...
            'FontSize',12,'FontWeight','bold','Units','normalized','Position',...
            pos.chSelBtn,'Callback', @selectGoodChannels,'Tag','chSelBtn');
        
        % pos.woiSelBtn = [(pos.chSelBtn(1)+pos.chSelBtn(3))*1.15,...
        %     pos.inspectBtn(2),0.15,pos.inspectBtn(4)];
        % woiSelBtn = uicontrol(mainFig,'Style','pushbutton','String','WOI selection',...
        %     'FontSize',14,'FontWeight','bold','Units','normalized','Position',...
        %     pos.woiSelBtn,'Callback', @selectGoodWOI);
        
        
        pos.HideMA = [pos.chSelBtn(1)+pos.chSelBtn(3)+myAxDim.xSep,...
            pos.inspectBtn(2),...
            0.1,...
            pos.inspectBtn(4)];
        uicontrol(main_fig,'Style','checkbox','String','Hide MAs',...
            'FontSize',12,'FontWeight','bold','Units','normalized','Position',...
            pos.HideMA,'Callback', @updateQPlots,'Tag','hideMACheckB');
        
        pos.SaveBtn = [pos.HideMA(1)+pos.HideMA(3)+myAxDim.xSep,...
            pos.inspectBtn(2),0.1,pos.inspectBtn(4)];
        uicontrol(main_fig,'Style','pushbutton','String','Save .nirs',...
            'FontSize',12,'FontWeight','bold','Units','normalized','Position',...
            pos.SaveBtn,'Callback', @save2dotnirs, 'Tag','saveBtn','Enable','off');
        
        pos.AdvView = [pos.SaveBtn(1)+pos.SaveBtn(3)+myAxDim.xSep,...
            pos.inspectBtn(2),...
            0.1,...
            pos.inspectBtn(4)];
        uicontrol(main_fig,'Style','checkbox','String','Advanced mode',...
            'FontSize',12,'FontWeight','bold','Units','normalized','Position',...
            pos.AdvView,'Callback', @updateQPlots,'Tag','advModeCheckB');
        %main_fig.Visible = 'on';
        
    end

%% -------------------------------------------------------------------------
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

%% -------------------------------------------------------------------------
    function inspectActive(source,event)
        
        nirsplot_param = getappdata(source.Parent,'nirsplot_parameters');
        n_channels = nirsplot_param.n_channels;
        qMats = nirsplot_param.quality_matrices;
        overlap_samples = qMats.overlap_samples;
        s = nirsplot_param.s;
        t = nirsplot_param.t;
        allowed_samp = qMats.allowed_samp;
        
        
        button = 0;
        flagDispWindow = false;
        pastChannel = 0;
        pastWindow = 0;
        while button ~= 27
            [iWindow,iChannel,button] = my_ginput(1);
            iWindow = floor(iWindow);
            iChannel = floor(iChannel);
            switch button
                case 1
                    flagDispWindow = true;
                case 3
                    %xLimWindow = [1,(qMats.sampPerWindow*qMats.n_windows)];
                     xLimWindow = [1,allowed_samp];
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
                %xLimWindow = [(qMats.sampPerWindow*(iWindow-1))+1,...
                %    (qMats.sampPerWindow*iWindow)];
                %xLimWindow = [((qMats.sampPerWindow-overlap_samples)*(iWindow-1))+overlap_samples+1,...
                %    ((qMats.sampPerWindow-overlap_samples)*iWindow+overlap_samples)];
                %xLimWindow = [(iWindow-1)*qMats.sampPerWindow-(iWindow-1)*(qMats.overlap_samples)+1,...
                %    iWindow*qMats.sampPerWindow-(iWindow-1)*(qMats.overlap_samples)];
                if nirsplot_param.overlap~=0
                    if mod(iWindow,2)==0
                        jj = floor((iWindow-1)/2)+1;
                        xLimWindow = [(jj-1)*qMats.sampPerWindow+1 , jj*qMats.sampPerWindow];
                        xLimWindow = xLimWindow + qMats.overlap_samples ;
                    else
                        jj = floor(iWindow/2)+1;
                        xLimWindow = [(jj-1)*qMats.sampPerWindow+1 , jj*qMats.sampPerWindow];
                    end
                else
                    xLimWindow = [(iWindow-1)*qMats.sampPerWindow+1 , iWindow*qMats.sampPerWindow];
                end
                
                %disp(['xLimWindow:',num2str(xLimWindow(1)),'-',num2str(xLimWindow(2))]);
            end
            if button ~=27
                if iChannel>0 && iChannel<=n_channels && iWindow>0 && iWindow<=qMats.n_windows
                    updateIPlot(source,iChannel,xLimWindow,iWindow,s,t);
                    pastWindow = iWindow;
                    pastChannel = iChannel;
                end
            end
        end
    end

%% -------------------------------------------------------------------------
    function updateQPlots(source,event)
        % UpdatePlot updates the quality plots with the 'qualityMats' input arg
        % bad channels are ploted according to 'plot_bad' flag
        nirsplot_param = getappdata(source.Parent,'nirsplot_parameters');
        qMats = nirsplot_param.quality_matrices;
        myAxes = nirsplot_param.main_fig_axes;
        n_channels = nirsplot_param.n_channels;
        overlap = nirsplot_param.overlap;
        woi = nirsplot_param.quality_matrices.woi;
        sclAlpha = nirsplot_param.sclAlpha;
        window_time = nirsplot_param.window;
        sciThrld = nirsplot_param.sci_threshold;
        pspThrld = nirsplot_param.psp_threshold;
        raw = getappdata(source.Parent,'rawNirs');
        
        mydivmap = [     0    0.2706    0.1608
                         0    0.4078    0.2157
                    0.1373    0.5176    0.2627
                    0.2549    0.6706    0.3647
                    0.4706    0.7765    0.4745
                    0.5765    0.8235    0.5176
                    0.6784    0.8667    0.5569
                    0.8510    0.9412    0.6392
                    0.9686    0.9882    0.7255
                    1.0000    1.0000    0.8980];
        %Unpacking
        sci_array = qMats.sci_array;
        power_array = qMats.power_array;
        combo_array = qMats.combo_array;
        combo_array_expanded = qMats.combo_array_expanded;
        
        woiMatrgb = zeros(n_channels,qMats.n_windows,3);
        woiMatrgb(:,:,:) = repmat(~woi.mat,1,1,3)*(hex2dec('bf')/255);
        alphaMat = ~woi.mat * sclAlpha;
        if strcmp(source.Tag,'hideMACheckB')
            hideMAVal = source.Value;
            advModeObj = source.Parent.findobj('Tag','advModeCheckB');
            advModeVal = advModeObj.Value;
        else
            advModeVal = source.Value;
            hideMA = source.Parent.findobj('Tag','hideMACheckB');
            hideMAVal = hideMA.Value;
        end
        
        %mygray = [0 0 0; repmat([0.7 0.7 0.7],100,1); 1 1 1];
        %mymap = [0 0 0;repmat([1 0 0],100,1);1 1 1 ];
        overlap_samples_ = nirsplot_param.window*nirsplot_param.fs*nirsplot_param.overlap;
        window_samples_ = nirsplot_param.window*nirsplot_param.fs;
        n_windows_ = (length(raw.t)-overlap_samples_)/(window_samples_-overlap_samples_);
        %colorb_sci = findobj('Tag','colorb_sci');
        
        %-
        % WE USE OPTION 1, BUT ONLY NEED TO "EXTEND" THE COMBO VIEW
        sci_mask = sci_array>=sciThrld;
        power_mask = power_array>=pspThrld;
        
        % Scalp Contact Index
        imagesc(myAxes.sci,sci_mask);
        colormap(myAxes.sci,[0 0 0; 1 1 1]);
        myAxes.sci.CLim = [0,1];
        colorbar(myAxes.sci,'eastoutside',...
            'Tag','colorb_sci',...
            'Ticks',[0.25 0.75],...
            'Limits',[0,1],...
            'TickLabels',{'Bad','Good'});
        % Power peak
        imagesc(myAxes.power,power_mask);
        colormap(myAxes.power,[0 0 0; 1 1 1]);
        myAxes.power.CLim = [0,1];
        colorbar(myAxes.power,'eastoutside',...
            'Tag','colorb_psp',...
            'Ticks',[0.25 0.75],...
            'Limits',[0,1],...
            'TickLabels',{'Bad','Good'});

        
        
        % Combo quality
        if ~hideMAVal
            imagesc(myAxes.combo,combo_array_expanded);
            myAxes.combo.CLim = [0, 4];
            %myAxes.combo.YLim =[1, n_channels];
            %myAxes.combo.XLim =[1, size(combo_array,2)];
            %colormap(myAxes.combo,[0 0 0;1 1 1]);
            % SCI,Power   combo_array_expanded    QualityColor
            %  0,0              0                   [0 0    0]
            %  0,1              1                   [0 0    0]
            %  2,0              2                   [1 0    0]
            %  2,1              3                   [1 1    1]
            %  X,X              4                   [0 0 1] Saturation
            qualityColor = [0 0 0; 0 0 0; 1 0 0; 1 1 1;0 0 1];
            colormap(myAxes.combo,qualityColor);
            colorbar(myAxes.combo,"eastoutside","Ticks",[0.7 1.0 2.05 2.75 3.5],...
                'TickLabels',...
                {[char(hex2dec('2717')),'SCI  ', char(hex2dec('2717')),'Power'],...
                [char(hex2dec('2717')),'SCI  ', char(hex2dec('2713')),'Power'],...
                [char(hex2dec('2713')),'SCI  ', char(hex2dec('2717')),'Power'],...
                [char(hex2dec('2713')),'SCI  ', char(hex2dec('2713')),'Power'],...
                'Saturation'});
        else
            % Combo quality
            imagesc(myAxes.combo,combo_array);
            myAxes.combo.CLim = [0, 1];
            %myAxes.combo.YLim =[1, n_channels];
            myAxes.combo.Colormap = [0 0 0;1 1 1];
            
            tickOffsetWind = 50/nirsplot_param.window;            
            ticksVals = (0:n_windows_)*window_samples_;
            ticksVals = ticksVals(1:tickOffsetWind:length(ticksVals));
            ticksVals = floor(ticksVals./window_samples_);
            ticksLab = 0:50:nirsplot_param.t(end);

            myAxes.combo.XAxis.TickValues=ticksVals(2:end);
            myAxes.combo.XAxis.TickLabels=split(num2str(ticksLab(2:end)));
            colorbar(myAxes.combo,"eastoutside","Ticks",[0.25 0.75],...
                'Limits',[0,1],'TickLabels',{'Bad','Good'});
        end
        
        
        myAxes.combo.YLabel.String = 'Channel #';
        myAxes.combo.YLabel.FontWeight = 'bold';
        % For figures MAs detected/undetected
        myAxes.combo.YAxis.TickValues = 1:n_channels;
        myAxes.combo.YAxis.TickLabels = num2cell([num2str(raw.SD.MeasList(1:length(raw.SD.Lambda):end,1)),repmat('-',n_channels,1),num2str(raw.SD.MeasList(1:length(raw.SD.Lambda):end,2))],2);
        myAxes.combo.YAxis.FontSize = 7;

        hold(myAxes.sci,'on');
        hold(myAxes.power,'on');
        hold(myAxes.combo,'on');
        
        % Drawing green bands
        imagesc(myAxes.sci,woiMatrgb,'AlphaData',alphaMat);
        %hold(myAxes.power,'on');
        imagesc(myAxes.power,woiMatrgb,'AlphaData',alphaMat);
        %hold(myAxes.combo,'on');
        imagesc(myAxes.combo,woiMatrgb,'AlphaData',alphaMat);
        
        %-
        if advModeVal
            % Scalp Contact Index
            imagesc(myAxes.sci,sci_array);
            colormap(myAxes.sci,mydivmap);
            colorbar(myAxes.sci,'eastoutside',...
                'Tag','colorb_sci',...
                'Ticks',[0 sciThrld 1],...
                'Limits',[0,1],...
                'TickLabels',{'0',['SCIThld:',num2str(sciThrld)],'1'});
            % Power peak
            imagesc(myAxes.power,power_array);
            colormap(myAxes.power,mydivmap);
            myAxes.power.CLim = [0,0.5];
            colorbar(myAxes.power,'eastoutside',...
                'Tag','colorb_psp',...
                'Ticks',[0 pspThrld 0.5],...
                'Limits',[0,0.5],...
                'TickLabels',{'0',['PSPThld:',num2str(pspThrld)],...
                '0.5'});             
        end
                
        %the next if-else sentence could be ereased
        
        
        
        % For visual consistency among axes
        myAxes.inspector.YLimMode = 'manual';
        myAxes.inspector.YLabel.String = 'Channel #';
        myAxes.inspector.XLabel.String = 'Time (s)';
        myAxes.inspector.YLabel.FontWeight = 'bold';
        myAxes.inspector.XLabel.FontWeight = 'bold';
        colorbar(myAxes.inspector,'Visible','off');
        
    end

%% -------------------------------------------------------------------------
    function updateIPlot(source,iChannel,xLimWindow,iWindow,s,t)
        disp(['channel:' num2str(iChannel)]);
        raw = getappdata(source.Parent,'rawNirs');
        nirsplot_param = getappdata(source.Parent,'nirsplot_parameters');
        qMats = nirsplot_param.quality_matrices;
        myAxes = nirsplot_param.main_fig_axes;
        n_channels = nirsplot_param.n_channels;
        conditions_mask = nirsplot_param.cond_mask;
        woi = nirsplot_param.quality_matrices.woi;
        fs = nirsplot_param.fs;
        fcut = nirsplot_param.fcut;
        window_time = nirsplot_param.window;
        overlap_samples = qMats.overlap_samples;
        
        rectangle_line_width = nirsplot_param.rectangle_line_width;
        sclAlpha = nirsplot_param.sclAlpha;
        dViewCheckb = findobj('Tag','hideMACheckB');
        
        myAxes.inspector.XLim= [t(xLimWindow(1)),t(xLimWindow(2))];
        YLimStd = [min(qMats.cardiac_data(:,xLimWindow(1):xLimWindow(2),iChannel),[],'all'),...
            max(qMats.cardiac_data(:,xLimWindow(1):xLimWindow(2),iChannel),[],'all')]*1.05;
        XLimStd = myAxes.inspector.XLim;
        
        % Normalization of cardiac_data between [a,b]
        a = -1;
        b =  1;
        cardiac_wl1_norm = a + (((qMats.cardiac_data(1,xLimWindow(1):xLimWindow(2),iChannel)-YLimStd(1)).*(b-a))./ (YLimStd(2)-YLimStd(1)));
        cardiac_wl2_norm = a + (((qMats.cardiac_data(2,xLimWindow(1):xLimWindow(2),iChannel)-YLimStd(1)).*(b-a))./ (YLimStd(2)-YLimStd(1)));
        YLimStd = [a,b].*1.05;
        cla(myAxes.inspector);
        
        plot(myAxes.inspector,t(xLimWindow(1):xLimWindow(2)),...
            cardiac_wl1_norm,'-b');
        hold(myAxes.inspector,'on');
        plot(myAxes.inspector,t(xLimWindow(1):xLimWindow(2)),...
            cardiac_wl2_norm,'-r');
        
        if(isfield(raw.SD,'Lambda'))
             WLs = raw.SD.Lambda(nirsplot_param.lambda_mask);
            strLgnds = {num2str(WLs(1)),num2str(WLs(2))};
        else
            strLgnds = {'\lambda 1','\lambda 2'};
        end
        
        updateQPlots(dViewCheckb,[]);
        
        %if (xLimWindow(2)-xLimWindow(1)+1) == (qMats.n_windows*qMats.sampPerWindow)
        if (xLimWindow(2)-xLimWindow(1)+1) == (qMats.n_windows*(qMats.sampPerWindow-overlap_samples)+overlap_samples)
            xRect = 0.5; %Because of the offset at the begining of a window
            yRect = iChannel-0.5;
            wRect = qMats.n_windows;
            hRect = 1;
            %poiMatrgb = zeros(n_channels,xLimWindow(2)-overlap_samples,3);
            %poiMatrgb(:,:,:) = repmat(repelem(~woi.mat(1,:),qMats.sampPerWindow-overlap_samples),n_channels,1,3).*(hex2dec('bf')/255);
            poiMatrgb=repmat(~woi.poi_mat,1,1,3).*(hex2dec('bf')/255);
            alphaMat = poiMatrgb(:,:,1) * sclAlpha;
            
            impoiMat = imagesc(myAxes.inspector,'XData',...
                [t(xLimWindow(1)),t(xLimWindow(2))],...
                'YData',YLimStd,'CData',poiMatrgb,'AlphaData',alphaMat);
            %ticksVals = linspace(0,XLimStd(2),8);
            ticksVals = (0:50:nirsplot_param.t(end));
            myAxes.inspector.XAxis.TickValues=ticksVals(2:end);
            %ticksLab = round(linspace(0,nirsplot_param.t(end),8));
            ticksLab = 0:50:nirsplot_param.t(end);
            myAxes.inspector.XAxis.TickLabels=split(num2str(ticksLab(2:end)));
            % Drawing onsets
            if ~ischar(conditions_mask)
                c = sum(conditions_mask);
                COI = find(conditions_mask);
                if c<9
                    colorOnsets = colorcube(8);
                else
                    colorOnsets = colorcube(c+1);
                    colorOnsets = colorOnsets(1:end-1,:);
                end
                
                for j=1:c
                    %mapping from 0,1 to 0,25%ofPeakToPeak
                    yOnset = (s(xLimWindow(1):xLimWindow(2),COI(j))*(YLimStd(2)-YLimStd(1))*0.25)-abs(YLimStd(1));
                    plot(myAxes.inspector,t(xLimWindow(1):xLimWindow(2)),...
                        yOnset,'LineWidth',2,...
                        'Color',colorOnsets(j,:));
                    strLgnds(2+j) = {['Cond ',num2str(COI(j))]};
                end
            end
            
        else
            %ticksVals = linspace(myAxes.inspector.XAxis.Limits(1),myAxes.inspector.XAxis.Limits(2),8);
            %ticksVals = linspace(myAxes.inspector.XAxis.Limits(1),myAxes.inspector.XAxis.Limits(2),5);
            ticksVals = myAxes.inspector.XAxis.Limits(1);
            ticksVals = round(ticksVals(1)): (round(ticksVals(1))+window_time);
            myAxes.inspector.XAxis.TickValues=ticksVals(1:end);
            ticksLab = round(ticksVals);
            myAxes.inspector.XAxis.TickLabels=split(num2str(ticksLab(1:end)));
            
            xRect = iWindow-0.5;
            yRect = iChannel-0.5;
            xQlabels = xRect + 1;
            yQlabels = yRect + 1;
            wRect = 1;
            hRect = 1;
            
            %fprintf('SCI:%.3f \t Power:%.3f\n',qMats.sci_array(iChannel,iWindow),qMats.power_array(iChannel,iWindow));
            textHAlign = 'left';
            textVAlign = 'top';
            if xRect > (qMats.n_windows/2)
                textHAlign = 'right';
                 xQlabels = xRect - 0.5;
            end
            if yRect > (n_channels/2)
                textVAlign = 'bottom';
                yQlabels = yRect - 0.5;
            end
            text(myAxes.power,xQlabels,yQlabels,num2str(qMats.power_array(iChannel,iWindow),'%.3f'),...
                'Color','red','FontSize',10,'FontWeight','bold','BackgroundColor','#FFFF00',...
                'Margin',1,'Clipping','on',...
                'HorizontalAlignment',textHAlign,'VerticalAlignment',textVAlign);
            text(myAxes.sci,xQlabels,yQlabels,num2str(qMats.sci_array(iChannel,iWindow),'%.3f'),...
                'Color','red','FontSize',10,'FontWeight','bold','BackgroundColor','#FFFF00',...
                'Margin',1,'Clipping','on',...
                'HorizontalAlignment',textHAlign,'VerticalAlignment',textVAlign);
            %--graphical debug
            % graphicDebug(qMats.cardiac_data(1,xLimWindow(1):xLimWindow(2),iChannel),...
            %     qMats.cardiac_data(2,xLimWindow(1):xLimWindow(2),iChannel),fs,fcut);
            % figure(source.Parent); 
        end
        myAxes.inspector.YLim = YLimStd;
        myAxes.inspector.XLim = XLimStd;
        %myAxes.inspector.YLabel.String = ['Channel ', num2str(iChannel)];
        myAxes.inspector.YLabel.String = ['S', num2str(raw.SD.MeasList(iChannel*length(raw.SD.Lambda),1)),'-D',num2str(raw.SD.MeasList(iChannel*length(raw.SD.Lambda),2)),'(',num2str(iChannel),')'];
        
        lgn = legend(myAxes.inspector,strLgnds,'Box','off','FontSize',10);
        
        rectangle(myAxes.combo,'Position',[xRect yRect wRect hRect],...
            'EdgeColor','m','FaceColor','none','Linewidth',rectangle_line_width);
        rectangle(myAxes.power,'Position',[xRect yRect wRect hRect],...
            'EdgeColor','m','FaceColor','none','Linewidth',rectangle_line_width);
        rectangle(myAxes.sci,'Position',[xRect yRect wRect hRect],...
            'EdgeColor','m','FaceColor','none','Linewidth',rectangle_line_width);
        
    end

%% -------------------------------------------------------------------------
    function [qualityMats] = qualityCompute(main_fig)
        raw = getappdata(main_fig,'rawNirs');
        nirsplot_param = getappdata(main_fig,'nirsplot_parameters');
        fcut = nirsplot_param.fcut;
        window = nirsplot_param.window;
        overlap = nirsplot_param.overlap;
        lambda_mask = nirsplot_param.lambda_mask;
        lambdas = nirsplot_param.lambdas;
        n_channels = nirsplot_param.n_channels;
        qltyThld = nirsplot_param.quality_threshold;
        sciThrld = nirsplot_param.sci_threshold;
        pspThrld = nirsplot_param.psp_threshold;
        scanInfo = nirsplot_param.dotNirsFile;
        dodFlag = nirsplot_param.dodFlag;
        if dodFlag == 1
            dm = mean(abs(raw.d),1);
            raw.d = exp(-raw.procResult.dod).*(ones(size(raw.d,1),1)*dm);
%            raw.d = raw.procResult.dod;
        end
        
        % Set the bandpass filter parameters
        %fs = 1/mean(diff(raw.t));
        fs = nirsplot_param.fs;
        n_samples = size(raw.d,1);
        fcut_min = fcut(1);
        fcut_max = fcut(2);
        if fcut_max >= (fs)/2
            fcut_max = (fs)/2 - eps;
            warning(['The highpass cutoff has been reduced from ',...
                num2str(fcut(2)), ' Hz to ', num2str(fcut_max),...
                ' Hz to satisfy the Nyquist sampling criterion']);
        end
        [B1,A1]=butter(1,[fcut_min*(2/fs) fcut_max*(2/fs)]);
        
        nirs_data = zeros(length(lambdas),n_samples,n_channels);
        cardiac_data = zeros(length(raw.SD.Lambda),n_samples,n_channels); % Lambdas x time x channels
        for j = 1:length(raw.SD.Lambda)
            % Filter everything but the cardiac component
            idx = find(raw.SD.MeasList(:,4) == j);
            nirs_data(j,:,:) = raw.d(:,idx);
            filtered_nirs_data=filtfilt(B1,A1,squeeze(nirs_data(j,:,:)));
            cardiac_data(j,:,:)=filtered_nirs_data./repmat(std(filtered_nirs_data,0,1),size(filtered_nirs_data,1),1); % Normalized heartbeat
        end
        overlap_samples = floor(window*fs*overlap);
        window_samples = floor(window*fs);
        if overlap ==0
            n_windows = floor((n_samples)/(window_samples));
        else % Valid for overlap=50%
            n_windows = 2*floor((n_samples)/(window_samples))-1;
        end
        cardiac_data = cardiac_data(find(lambda_mask),:,:);
        sci_array = zeros(size(cardiac_data,3),n_windows);    % Number of optode is from the user's layout, not the machine
        power_array = zeros(size(cardiac_data,3),n_windows);
        fpower_array = zeros(size(cardiac_data,3),n_windows);
        cardiac_windows = zeros(length(lambdas),window_samples,n_channels,n_windows);
        for j = 1:n_windows
            if overlap~=0
                if mod(j,2)==0
                    jj = floor((j-1)/2)+1;
                    interval = (jj-1)*window_samples+1 : jj*window_samples;
                    interval = interval + overlap_samples ;
                else
                    jj = floor(j/2)+1;
                    interval = (jj-1)*window_samples+1 : jj*window_samples;
                end
            else 
               interval = (j-1)*window_samples+1 : j*window_samples;
            end
            cardiac_windows(:,:,:,j) = cardiac_data(:,interval,:);
%             if j<5 || j>(n_windows-4)
%                 disp(['interval(',num2str(j),'):',num2str(interval(1)),'-',num2str(interval(end))]);
%             end
        end

        for j = 1:n_windows
            cardiac_window = cardiac_windows(:,:,:,j);
            sci_array_channels = zeros(1,size(cardiac_window,3));
            power_array_channels = zeros(1,size(cardiac_window,3));
            fpower_array_channels = zeros(1,size(cardiac_window,3));
            for k = 1:size(cardiac_window,3) % Channels iteration
                %cross-correlate the two wavelength signals - both should have cardiac pulsations
                similarity = xcorr(squeeze(cardiac_window(1,:,k)),squeeze(cardiac_window(2,:,k)),'unbiased');
                if any(abs(similarity)>eps)
                    % this makes the SCI=1 at lag zero when x1=x2 AND makes the power estimate independent of signal length, amplitude and Fs
                    similarity = length(squeeze(cardiac_window(1,:,k)))*similarity./sqrt(sum(abs(squeeze(cardiac_window(1,:,k))).^2)*sum(abs(squeeze(cardiac_window(2,:,k))).^2));
                    similarity(isnan(similarity)) = 0;
                    [pxx,f] = periodogram(similarity,hamming(length(similarity)),length(similarity),fs,'power');
                    [pwrest,idx] = max(pxx(f<fcut_max)); % FIX Make it age-dependent
                    sci=similarity(length(squeeze(cardiac_window(1,:,k))));
                    power=pwrest;
                    fpower=f(idx);
                    sci_array_channels(k) = sci;
                    power_array_channels(k) = power;
                    fpower_array_channels(k) = fpower;
                else
                    warning('Similarity results close to zero');           
                    sci_array_channels(k) = 0;
                    power_array_channels(k) = 0;
                    fpower_array_channels(k) = -1;
                end
                
            end
            sci_array(:,j) = sci_array_channels;    % Adjust not based on machine
            power_array(:,j) = power_array_channels;
            fpower_array(:,j) = fpower_array_channels;
        end
        
        % Summary analysis
        if (1)
            [woi,allowed_samp] = getWOI(window_samples,n_windows,overlap_samples,nirsplot_param);
            idxPoi = logical(woi.mat(1,:));
        else
            [woi,allowed_samp] = getWOI(window_samples,n_windows,overlap_samples,nirsplot_param);
            woi.mat = ones(n_channels,n_windows)
            allowed_samp = ones(1,window_samples*n_windows);
            woi.poi_mat = repmat(allowed_samp,n_channels,1);
            idxPoi = logical(woi.mat(1,:));
        end
        
        mean_sci_link  = mean(sci_array(:,idxPoi),2);
        std_sci_link  = std(sci_array(:,idxPoi),0,2);
        good_sci_link = sum(sci_array(:,idxPoi)>=sciThrld,2)/size(sci_array(:,idxPoi),2);
        mean_sci_window  = mean(sci_array(:,idxPoi),1);
        std_sci_window  = std(sci_array(:,idxPoi),0,1);
        good_sci_window = sum(sci_array(:,idxPoi)>=sciThrld,1)/size(sci_array(:,idxPoi),1);
        
        mean_power_link  = mean(power_array(:,idxPoi),2);
        std_power_link  = std(power_array(:,idxPoi),0,2);
        good_power_link = sum(power_array(:,idxPoi)>=pspThrld,2)/size(power_array(:,idxPoi),2);
        mean_power_window  = mean(power_array(:,idxPoi),1);
        std_power_window  = std(power_array(:,idxPoi),0,1);
        good_power_window = sum(power_array(:,idxPoi)>=pspThrld,1)/size(power_array(:,idxPoi),1);
        
        combo_array = (sci_array >= sciThrld) & (power_array >= pspThrld);
        saturat_mat = fpower_array==-1;
        combo_array_expanded = 2*(sci_array >= sciThrld) + (power_array >= pspThrld) +...
            saturat_mat*4;
        
        mean_combo_link  = mean(combo_array,2);
        std_combo_link  = std(combo_array,0,2);
        
        good_combo_link  = mean(combo_array(:,idxPoi),2);
        mean_combo_window  = mean(combo_array,1);
        std_combo_window  = std(combo_array,0,1);
        
        idx_gcl = good_combo_link>=qltyThld;
        good_combo_window = mean(combo_array(idx_gcl,:),1);
        
        % Detecting experimental blocks
        exp_blocks = zeros(1,length(woi.start));
        for iblock = 1:length(woi.start)
            block_start_w = woi.start(iblock);
            block_end_w = woi.end(iblock);
            exp_blocks(iblock) = mean(good_combo_window(block_start_w:block_end_w));
        end
        
        % Detect artifacts and bad links
        bad_links = find(mean_combo_link<qltyThld);
        bad_windows = find(mean_combo_window<qltyThld);
        
        % Packaging sci, peakpower and combo
        qualityMats.sci_array    = sci_array;
        qualityMats.power_array  = power_array;
        qualityMats.combo_array  = combo_array;
        qualityMats.combo_array_expanded = combo_array_expanded;
        qualityMats.bad_links    = bad_links;
        qualityMats.bad_windows  = bad_windows;
        qualityMats.sampPerWindow = window_samples;
        qualityMats.fs = fs;
        qualityMats.n_windows = n_windows;
        qualityMats.overlap_samples = overlap_samples;
        qualityMats.cardiac_data = cardiac_data; 
        % Bug found by Benjamin Zinszer (miscalculation in the S-D pairs)
        qualityMats.good_combo_link = [raw.SD.MeasList(1:length(raw.SD.Lambda):end,1),...
            raw.SD.MeasList(1:length(raw.SD.Lambda):end,2),good_combo_link];
        qualityMats.good_combo_window = good_combo_window;
        qualityMats.woi = woi;
        qualityMats.allowed_samp = allowed_samp;
        qualityMats.MeasListAct = repelem(idx_gcl,2);%[idx_gcl; idx_gcl];
        qualityMats.MeasList = raw.SD.MeasList;
        qualityMats.thresholds.sci = sciThrld;
        qualityMats.thresholds.peakpower = pspThrld;
        qualityMats.thresholds.quality = qltyThld;
        qualityMats.scanInfo = scanInfo ;
        %
    end




%% -------------------------------------------------------------------------
    function [woi,allowed_samp] = getWOI(window_samples,n_windows,overlap_samples,nirsplot_parameters)
        % Assuming no overlaped trials
        % The maximum number of allowed samples is window_samples*n_windows to consider
        % an integer number of windows, module(total_samples,n_windows) = 0
        
        fs = nirsplot_parameters.fs;
        s = nirsplot_parameters.s;
        t = nirsplot_parameters.t;
        mergewoi_flag = nirsplot_parameters.mergewoi_flag;
        n_channels = nirsplot_parameters.n_channels;        
        %allowed_samp = n_windows*window_samples;
        %allowed_samp = n_windows*(window_samples-overlap_samples)+overlap_samples;
        if nirsplot_parameters.overlap ~= 0 
            allowed_samp = (n_windows*window_samples+window_samples)/2;
        else
            allowed_samp = n_windows*window_samples;
        end
        if strcmp(nirsplot_parameters.cond_mask,'all')
            conditions_mask = ones(1,size(s,2));
        elseif strcmp(nirsplot_parameters.cond_mask,'resting')
            woi = struct;
            woi.mat = ones(n_channels,n_windows);
            woi.poi_mat = ones(n_channels,allowed_samp,1);
            woi.start = 1;
            woi.end = n_windows;
            return
        else
            conditions_mask = logical(nirsplot_parameters.cond_mask);
        end
        
        poi = sum(s(1:allowed_samp,conditions_mask),2);
        poi = poi(1:allowed_samp);
        % Sometimes 's' variable encodes the stimuli durations by including consecutive
        % values of 1. We are interested on the onsets, then we remove consecutive ones.
        idxpoi = find(poi);
        poi = zeros(size(poi));
        poi(idxpoi(diff([0;idxpoi])>1)) = 1;
        nOnsets = length(find(poi));
        idxStim = find(poi);
        interOnsetTimes = t(idxStim(2:end)) - t(idxStim(1:end-1));
        medIntTime = median(interOnsetTimes);
        iqrIntTime = iqr(interOnsetTimes);
        %blckDurTime = (medIntTime/2) + (0.5*iqrIntTime);
        blckDurTime = medIntTime + (0.5*iqrIntTime);
        blckDurSamp = round(fs*blckDurTime);
        blckDurWind = floor(blckDurSamp/((window_samples-overlap_samples)));% floor(blckDurSamp/window_samples);
        woi = struct('mat',zeros(n_channels,n_windows),...
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
            startWOI = ceil(startPOI/((window_samples-overlap_samples)));%floor(startPOI/window_samples);
            if startWOI==0
                startWOI = 1;
            end
            
            endPOI = idxStim(i)+blckDurSamp;
            if endPOI > allowed_samp
                endPOI = allowed_samp;
            end
            endWOI = floor(endPOI/((window_samples-overlap_samples)));%ceil(endPOI/window_samples);
            poi(startPOI:endPOI) = 1;
            woi_array(startWOI:endWOI) = 1;
            woi.start(i) = startWOI;
            woi.end(i) = endWOI;
        end
        
        % See my comment about the preference of WOIs rather than of POIs, if POI
        % information is needed, uncomment next two lines and return POIs variables
         poi = poi';
         poiMat_ = repmat(poi,n_channels,1);
        
        woiblank = 0;
        idxInit = [];
        woitmp = woi_array;
        
        % If the gap's duration between two consecutives blocks of interest is less than the
        % block's average duration, then those two consecutives blocks will merge.
        % This operation has a visual effect (one bigger green block instead of
        % two green blocks with a small gap in between), and for quality
        % results, the windows inside of such a gap are now considered for quality computation.
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
        if mergewoi_flag == true
            woi_array = woitmp;
        end
        woi.mat = repmat(woi_array,n_channels,1);
        woi.poi_mat = poiMat_;
    end

%% -------------------------------------------------------------------------
    function dotNirsOutput = selectGoodChannels(source, events)
        bpGoodQuality(source.Parent);
        uiwait(source.Parent);
        nirsplot_param = getappdata(source.Parent,'nirsplot_parameters');
        if isfield(nirsplot_param.quality_matrices,'active_channels')
            disp(['Threshold was changed to ',num2str(nirsplot_param.quality_threshold)]);
            saveBtn = findobj('Tag','saveBtn');
            saveBtn.Enable = 'on';
        end
        dotNirsOutput = 0;
    end

%% -------------------------------------------------------------------------
%!This function is not tested yet!
    function report_table = saveQuality()
        qMats = getappdata(source.Parent,'qualityMats');
        
        report_table = table({qMats.bad_links'}, {qMats.bad_windows});
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

%% -------------------------------------------------------------------------
    function saving_status = save2dotnirs(source, events)
        nirsplot_param = getappdata(source.Parent,'nirsplot_parameters');
        active_channels = nirsplot_param.quality_matrices.active_channels;
        dotNirsPath = nirsplot_param.dotNirsPath;
        dotNirsFileName = nirsplot_param.dotNirsFile;
        dotNirsExt = nirsplot_param.dotNirsExt;
        outputFolder = 'qtnirs_proc';
        if ~exist(outputFolder,'dir')
            mkdir(outputFolder);            
        end
        switch dotNirsExt
            case '.nirs'
                 raw = getappdata(source.Parent,'rawNirs');                 
                 raw.SD.MeasListAct = [active_channels; active_channels];    
                 raw.tIncMan = ones(length(raw.t),1);
                 if ~isfield(raw, 'aux')
                     raw.aux = [];
                 end
                 save([dotNirsPath,filesep,outputFolder,filesep,dotNirsFileName,'_qt-proc.nirs'],...
                    '-struct','raw','-MAT');
            case '.snirf'
                 raw = getappdata(source.Parent,'rawNirs');                 
                 raw.SD.MeasListAct = [active_channels; active_channels];    
                 raw.tIncMan = ones(length(raw.t),1);               
                 if ~isfield(raw, 'aux')
                     raw.aux = [];
                 end
                 snirf_saved = SnirfClass(raw);
                 snirf_saved.Save([dotNirsPath,filesep,outputFolder,filesep,dotNirsFileName,'_qt-proc.snirf']); 
        end
        
        
        %Notify to the user if the new file was succesfully created
        saving_status = exist([dotNirsPath,filesep,outputFolder,filesep,dotNirsFileName,'_qt-proc',dotNirsExt],'file');
        if saving_status
            msgbox('Operation Completed','Success');
        else
            msgbox('Operation Failed','Error');
        end
        
    end

%%
    function graphicDebug(window1,window2,fs,fcut)
        %cross-correlate the two wavelength signals - both should have cardiac pulsations
        [similarity,lags] = xcorr(window1,window2,'unbiased');
        if any(abs(similarity)>eps)
            % this makes the SCI=1 at lag zero when x1=x2 AND makes the power estimate independent of signal length, amplitude and Fs
            similarity = length(window1)*similarity./sqrt(sum(abs(window1).^2)*sum(abs(window2).^2));
        else
            warning('Similarity results close to zero');
        end
        [pxx,f] = periodogram(similarity,hamming(length(similarity)),length(similarity),fs,'power');
        f3=figure(3);
        clf(f3);
        subplot(1,3,1);
        plot(window1-2,'k-');
        hold on;
        plot(window2+2,'k--');
        ylabel('Raw intensity','FontSize',15);% x_{\lambda_1}, x_{\lambda_2}');
        legend('$x_{\lambda_1}$', '$x_{\lambda_2}$','Interpreter','latex','FontSize',16);
        ylim([-7 7]);
        subplot(1,3,2);
        plot(f3.Children(1),lags,similarity,'k-');
        yline(0.8,'r--');
        ylim([-1 1]);
        ylabel('$\bar{x}_{\lambda_1} \otimes \bar{x}_{\lambda_2}$','Interpreter','latex','FontSize',16);
        subplot(1,3,3);
        plot(f3.Children(1),f,pxx,'k-');
        ylabel('$F(x_{\lambda_1} \otimes x_{\lambda_2})$','Interpreter','latex','FontSize',16);
        ylim([0 0.5]);
        yline(0.1,'r--');
    end


%% -------------------------------------------------------------------------
    function [qualityMats] = qualityCompute2(main_fig)
        raw = getappdata(main_fig,'rawNirs');
        nirsplot_param = getappdata(main_fig,'nirsplot_parameters');
        fcut = nirsplot_param.fcut;
        window = nirsplot_param.window;
        overlap = nirsplot_param.overlap;
        lambda_mask = nirsplot_param.lambda_mask;
        lambdas = nirsplot_param.lambdas;
        n_channels = nirsplot_param.n_channels;
        qltyThld = nirsplot_param.quality_threshold;
        sciThrld = nirsplot_param.sci_threshold;
        pspThrld = nirsplot_param.psp_threshold;
        
        % Set the bandpass filter parameters
        %fs = 1/mean(diff(raw.t));
        fs = nirsplot_param.fs;
        fcut_min = fcut(1);
        fcut_max = fcut(2);
        if fcut_max >= (fs)/2
            fcut_max = (fs)/2 - eps;
            warning(['The highpass cutoff has been reduced from ',...
                num2str(fcut(2)), ' Hz to ', num2str(fcut_max),...
                ' Hz to satisfy the Nyquist sampling criterion']);
        end
        [B1,A1]=butter(1,[fcut_min*(2/fs) fcut_max*(2/fs)]);
        overlap_samples = floor(window*fs*overlap);
        window_samples = floor(window*fs);
        n_windows = floor((size(raw.d,1)-overlap_samples)/(window_samples-overlap_samples));

        cardiac_windows = zeros(length(lambdas),window_samples,n_channels,n_windows);
        nirs_data = zeros(length(lambdas),window_samples,n_channels);
        cardiac_data = zeros(length(lambdas),window_samples,n_channels); % #Lambdas x #Windows x #Channels
        for lam = 1:length(lambdas)
            for j = 1:n_windows
                interval = (j-1)*window_samples-(j-1)*(overlap_samples)+1 : j*window_samples-(j-1)*(overlap_samples);
                idx = raw.SD.MeasList(:,4) == lambdas(lam);
                nirs_data(lam,:,:) = raw.d(interval,idx);
                filtered_nirs_data=filtfilt(B1,A1,squeeze(nirs_data(lam,:,:)));
                cardiac_data(lam,:,:)=filtered_nirs_data./repmat(std(filtered_nirs_data,0,1),size(filtered_nirs_data,1),1); % Normalized heartbeat
                cardiac_windows(lam,:,:,j) = cardiac_data(lam,:,:);
            end
        end
        overlap_samples = floor(window*fs*overlap);
        window_samples = floor(window*fs);
        n_windows = floor((size(raw.d,1)-overlap_samples)/(window_samples-overlap_samples));
        cardiac_data = cardiac_data(lambda_mask,:,:);
        sci_array = zeros(size(cardiac_data,3),n_windows);    % Number of optode is from the user's layout, not the machine
        power_array = zeros(size(cardiac_data,3),n_windows);
        parfor j = 1:n_windows
            cardiac_window = cardiac_windows(:,:,:,j);
            sci_array_channels = zeros(1,size(cardiac_window,3));
            power_array_channels = zeros(1,size(cardiac_window,3));
            fpower_array_channels = zeros(1,size(cardiac_window,3));
            for k = 1:size(cardiac_window,3) % Channels iteration
                %cross-correlate the two wavelength signals - both should have cardiac pulsations
                similarity = xcorr(squeeze(cardiac_window(1,:,k)),squeeze(cardiac_window(2,:,k)),'unbiased');
                if any(abs(similarity)>eps)
                    % this makes the SCI=1 at lag zero when x1=x2 AND makes the power estimate independent of signal length, amplitude and Fs
                    similarity = length(squeeze(cardiac_window(1,:,k)))*similarity./sqrt(sum(abs(squeeze(cardiac_window(1,:,k))).^2)*sum(abs(squeeze(cardiac_window(2,:,k))).^2));
                else
                    warning('Similarity results close to zero');
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
        [woi] = getWOI(window_samples,n_windows,nirsplot_param);
        idxPoi = logical(woi.mat(1,:));
        
        mean_sci_link  = mean(sci_array(:,idxPoi),2);
        std_sci_link  = std(sci_array(:,idxPoi),0,2);
        good_sci_link = sum(sci_array(:,idxPoi)>=sciThrld,2)/size(sci_array(:,idxPoi),2);
        mean_sci_window  = mean(sci_array(:,idxPoi),1);
        std_sci_window  = std(sci_array(:,idxPoi),0,1);
        good_sci_window = sum(sci_array(:,idxPoi)>=sciThrld,1)/size(sci_array(:,idxPoi),1);
        
        mean_power_link  = mean(power_array(:,idxPoi),2);
        std_power_link  = std(power_array(:,idxPoi),0,2);
        good_power_link = sum(power_array(:,idxPoi)>=pspThrld,2)/size(power_array(:,idxPoi),2);
        mean_power_window  = mean(power_array(:,idxPoi),1);
        std_power_window  = std(power_array(:,idxPoi),0,1);
        good_power_window = sum(power_array(:,idxPoi)>=pspThrld,1)/size(power_array(:,idxPoi),1);
        
        combo_array = (sci_array >= sciThrld) & (power_array >= pspThrld);
        combo_array_expanded = 2*(sci_array >= sciThrld) + (power_array >= pspThrld);
                
        mean_combo_link  = mean(combo_array,2);
        std_combo_link  = std(combo_array,0,2);
        
        good_combo_link  = mean(combo_array(:,idxPoi),2);
        mean_combo_window  = mean(combo_array,1);
        std_combo_window  = std(combo_array,0,1);
        
        idx_gcl = good_combo_link>=qltyThld;
        good_combo_window = mean(combo_array(idx_gcl,:),1);
        
        % Detecting experimental blocks
        exp_blocks = zeros(1,length(woi.start));
        for iblock = 1:length(woi.start)
            block_start_w = woi.start(iblock);
            block_end_w = woi.end(iblock);
            exp_blocks(iblock) = mean(good_combo_window(block_start_w:block_end_w));
        end
        
        % Detect artifacts and bad links
        bad_links = find(mean_combo_link<qltyThld);
        bad_windows = find(mean_combo_window<qltyThld);
        
       % Packaging sci, peakpower and combo
        qualityMats.sci_array    = sci_array;
        qualityMats.power_array  = power_array;
        qualityMats.combo_array  = combo_array;
        qualityMats.combo_array_expanded = combo_array_expanded;
        qualityMats.bad_links    = bad_links;
        qualityMats.bad_windows  = bad_windows;
        qualityMats.sampPerWindow = window_samples;
        qualityMats.fs = fs;
        qualityMats.n_windows = n_windows;
        qualityMats.overlap_samples = overlap_samples;
        qualityMats.cardiac_data = cardiac_data; 
        % Bug found by Benjamin Zinszer (miscalculation in the S-D pairs)
        qualityMats.good_combo_link = [raw.SD.MeasList(1:length(raw.SD.Lambda):end,1),...
            raw.SD.MeasList(1:length(raw.SD.Lambda):end,2),good_combo_link];
        qualityMats.good_combo_window = good_combo_window;
        qualityMats.woi = woi;
        qualityMats.allowed_samp = allowed_samp;
        qualityMats.MeasListAct = repelem(idx_gcl,2);%[idx_gcl; idx_gcl];
        qualityMats.MeasList = raw.SD.MeasList;
        qualityMats.thresholds.sci = sciThrld;
        qualityMats.thresholds.peakpower = pspThrld;
        qualityMats.thresholds.quality = qltyThld;
        %
    end
end %end of nirsplot function definition
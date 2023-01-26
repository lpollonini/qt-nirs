classdef qtnirsLoadFileGUI < matlab.apps.AppBase
    
    % Properties that correspond to app components
    properties (Access = public)
        NIRSPlotGUIUIFigure     matlab.ui.Figure
        NIRSPlotLabel           matlab.ui.control.Label
        LoadButton              matlab.ui.control.Button
        SourcesLabel            matlab.ui.control.Label
        DetectorsLabel          matlab.ui.control.Label
        ChannelsLabel           matlab.ui.control.Label
        DurLabel                matlab.ui.control.Label
        nSrcLab                 matlab.ui.control.Label
        nDetecLab               matlab.ui.control.Label
        nChanLab                matlab.ui.control.Label
        secDurLab               matlab.ui.control.Label
        nirsfileEditFieldLabel  matlab.ui.control.Label
        nirsfileEditField       matlab.ui.control.EditField
        FreqMinEditFieldLabel   matlab.ui.control.Label
        FreqMinEditField        matlab.ui.control.NumericEditField
        FreqMaxEditFieldLabel   matlab.ui.control.Label
        FreqMaxEditField        matlab.ui.control.NumericEditField
        LengthsecSpinnerLabel   matlab.ui.control.Label
        LengthsecSpinner        matlab.ui.control.Spinner
        WindowLabel             matlab.ui.control.Label
        CutoffFrequencyLabel    matlab.ui.control.Label
        PlotButton              matlab.ui.control.Button
        OverlappingLabel        matlab.ui.control.Label
        WindowsOverlapCtrl      matlab.ui.control.CheckBox
        QualityThresholdLabel   matlab.ui.control.Label
        QualityThresholdField   matlab.ui.control.NumericEditField
        condCheckBoxes          matlab.ui.control.CheckBox
        condCheckBLabel         matlab.ui.control.Label
        treeDotNirs             matlab.ui.container.Tree
        treeNodesDotNirs        matlab.ui.container.TreeNode
        dodCheckBox             matlab.ui.control.CheckBox
        wlLabel                 matlab.ui.control.Label
        wlCheckBoxes            matlab.ui.control.CheckBox
        scipspLabel             matlab.ui.control.Label
        sciThresholdLabel       matlab.ui.control.Label
        pspThresholdLabel       matlab.ui.control.Label
        sciThresholdField       matlab.ui.control.NumericEditField
        pspThresholdField       matlab.ui.control.NumericEditField
        bg                      matlab.ui.container.ButtonGroup
        rb1                     matlab.ui.control.RadioButton
        rb2                     matlab.ui.control.RadioButton  
    end
    
    
    properties (Access = public)
        rawDotNirs  % This is the content from .nirs files loaded as '-mat' format.
        nSources    % Number of sources
        nDetectors  % Number of detectors
        nChannels   % Number of channels
        nCond       % Number of conditions
        sampDur     % Number of data samples
        secDur      % Time in seconds
    end
    
    properties (Access = private)
        bpFmin % [bpFmin bpFmax] representing the bandpass of the cardiac pulsation (default [0.5 2.5])
        bpFmax % [bpFmin bpFmax] representing the bandpass of the cardiac pulsation (default [0.5 2.5])
        windowSec % length in seconds of the window to partition the signal with (defaut: 5)
        windowOverlap % fraction overlap (0..0.99) between adjacent windows (default: 0, no overlap)
        reportTable % Quality output table report
        quality_threshold % Minimum required quality for channels
        sci_threshold % SCI threshold (default:0.8)
        psp_threshold % PSP threshold (default:0.1)
        condMask % Binary mask for selecting the conditions of interest
        lambdaMask % Binary mask for selection the WLs to be used for computing quality metrics.
        dotNirsPath % Path of the .nirs file to analyze
        dotNirsFile % Path of the .nirs file to analyze
        nlambda     % Number of WLs in the scan file
    end
    
    
    % Callbacks that handle component events
    methods (Access = private)
        
        function startupFcn(app, input_param)
            if ischar(input_param)
                app.dotNirsPath = input_param;
                 app.dotNirsFile = [];
                 
            elseif isstruct(input_param)
                if isfield(input_param,'dotNirsFile')
                    app.dotNirsPath = input_param.dotNirsPath;
                    app.dotNirsFile = [input_param.dotNirsFile '.nirs'];

                    
                end
                if isfield(input_param,'n_sources')
                    app.nSources        = input_param.n_sources;
                    app.nSrcLab.Text    = num2str(input_param.n_sources);
                end
                if isfield(input_param,'n_detectors')
                    app.nDetectors      = input_param.n_detectors;
                    app.nDetecLab.Text  = num2str(input_param.n_detectors);
                end
                if isfield(input_param,'n_channels')
                    app.nChannels       = input_param.n_channels;
                    app.nChanLab.Text   = num2str(input_param.n_channels);
                end
                if isfield(input_param,'t')
                    app.secDur          = input_param.t(end);
                    app.secDurLab.Text  = num2str(app.secDur,'%.2f');
                end
                if isfield(input_param,'dodFlag')
                    if input_param.dodFlag ~= -1
                        app.dodCheckBox.Enable = 'on';
                        app.dodCheckBox.Value = input_param.dodFlag;
                    else
                        app.dodCheckBox.Enable = 'off';
                        app.dodCheckBox.Value = 0;
                    end
                else
                    app.dodCheckBox.Enable = 'off';
                end
                if isfield(input_param,'fcut')
                    app.bpFmin = min(input_param.fcut);
                    app.bpFmax = max(input_param.fcut);
                    app.FreqMinEditField.Value = app.bpFmin;
                    app.FreqMaxEditField.Value = app.bpFmax;
                end
                if isfield(input_param,'window')
                    app.windowSec = input_param.window;
                    app.LengthsecSpinner.Value = app.windowSec;
                end
                if isfield(input_param,'overlap')
                    app.windowOverlap = input_param.overlap;
                    %app.WindowsoverlapSlider.Value = app.windowOverlap*100;
                    app.WindowsOverlapCtrl.Value = app.windowOverlap>0;
                end
                if isfield(input_param,'quality_threshold')
                    app.quality_threshold = input_param.quality_threshold;
                    app.QualityThresholdField.Value = app.quality_threshold*100;
                end
                if isfield(input_param,'sci_threshold')
                    app.sci_threshold = input_param.sci_threshold;
                    app.sciThresholdField.Value = app.sci_threshold;
                end
                if isfield(input_param,'psp_threshold')
                    app.psp_threshold = input_param.psp_threshold;
                    app.pspThresholdField.Value = app.psp_threshold;
                end
                if isfield(input_param,'cond_mask')
                    % Create and fill StimuliSelectionBoxes
                    xOffset = 16;
                    yOffset = 145;
                    app.condMask = input_param.cond_mask;
                    if isnumeric(app.condMask)
                        app.nCond = length(input_param.cond_mask);
                        app.condCheckBoxes.delete;
                        for iCB=1:app.nCond
                            app.condCheckBoxes(iCB) = uicheckbox(app.NIRSPlotGUIUIFigure);
                            app.condCheckBoxes(iCB).Position =[xOffset,yOffset,...
                                35, 15];
                            app.condCheckBoxes(iCB).Value = input_param.cond_mask(iCB);
                            app.condCheckBoxes(iCB).Text = num2str(iCB);
                            mod_ = mod(iCB,5);
                            if mod_ == 0
                                yOffset = yOffset - 25;
                                xOffset = 16;
                            else
                                xOffset = xOffset+40;
                            end
                        end
                    else
                        app.nCond = size(input_param.s,2);
                    end
                end
            end
            LoadDotNirsDir(app, app.dotNirsPath);
        end
        
        % Tree selectionChanged function: treeDotNirs
        function treeSelection(app, event)
            if isempty(event.SelectedNodes.Children) %Is it a leaf (file)?
                app.PlotButton.Enable = 1;
                dotNirsFileSel = event.SelectedNodes.Text;
                app.dotNirsFile = dotNirsFileSel;
                disp(['Loading: ',app.dotNirsPath,filesep,dotNirsFileSel]);
                [filepath,name,ext] = fileparts([app.dotNirsPath,filesep,dotNirsFileSel]);
                switch ext
                    case '.nirs'
                        app.rawDotNirs = load([app.dotNirsPath,filesep,app.dotNirsFile],'-mat');
                    case '.snirf'
                        rawsnirf = SnirfClass([app.dotNirsPath,filesep,app.dotNirsFile]);
                        %app.rawDotNirs.d = rawsnirf.Get_d;
                        %app.rawDotNirs.s = rawsnirf.Get_s;
                        %app.rawDotNirs.t = rawsnirf.Get_t;
                        app.rawDotNirs.t = rawsnirf.data.time;
                        app.rawDotNirs.s = rawsnirf.GetStims(app.rawDotNirs.t);
                        app.rawDotNirs.SD = rawsnirf.Get_SD;
                    otherwise
                        error('The input file should be a .nirs file format');
                end
                app.nlambda         = length(app.rawDotNirs.SD.Lambda);
                app.nChannels       = size(app.rawDotNirs.SD.MeasList,1)/app.nlambda;
                app.nSources        = size(app.rawDotNirs.SD.SrcPos,1);
                app.nDetectors      = size(app.rawDotNirs.SD.DetPos,1);
                app.secDur          = app.rawDotNirs.t(end);
                app.nSrcLab.Text    = num2str(app.nSources);
                app.nDetecLab.Text  = num2str(app.nDetectors);
                app.nChanLab.Text   = num2str(app.nChannels);
                app.secDurLab.Text  = num2str(app.secDur,'%.2f');
                
                % Checking for OD data
                if isfield(app.rawDotNirs,'procResult')
                    if ~isempty(app.rawDotNirs.procResult.dod)
                        app.dodCheckBox.Enable = 'on';
                    end
                else
                    app.dodCheckBox.Value = 0;
                    app.dodCheckBox.Enable = 'off';
                end
                
                if ~isfield(app.rawDotNirs,'s')
                    frequency_samp = 1/mean(diff(app.rawDotNirs.t));
                    if isfield(app.rawDotNirs,'StimDesign')
                        nStim = length(app.rawDotNirs.StimDesign);
                        sTmp = zeros(size(app.rawDotNirs.d,1),nStim);
                        for iStim = 1:nStim
                            sTmp(floor(app.rawDotNirs.StimDesign(iStim).onset * frequency_samp),iStim) = 1;
                        end
                        app.rawDotNirs.s = sTmp;
                        clear sTmp;
                    else
                        error('Stimuli information is not available.');
                    end
                end
                
                % Create StimuliSelectionBoxes
                xOffset = 14;
                yOffset = 145;
                app.nCond = size(app.rawDotNirs.s,2);
                app.condCheckBoxes.delete;
                % 'Resting' option
                app.condCheckBoxes(1) = uicheckbox(app.NIRSPlotGUIUIFigure);
                app.condCheckBoxes(1).Position =[xOffset,yOffset,...
                    60, 15];
                app.condCheckBoxes(1).Value = 1;
                app.condCheckBoxes(1).Text = 'Resting';
                app.condCheckBoxes(1).ValueChangedFcn = createCallbackFcn(app,@cleanCheckboxes,true);
                xOffset = xOffset + 80;
                
                for iCB=2:(app.nCond +1)
                    app.condCheckBoxes(iCB) = uicheckbox(app.NIRSPlotGUIUIFigure);
                    app.condCheckBoxes(iCB).Position =[xOffset,yOffset,...
                        35, 15];
                    app.condCheckBoxes(iCB).Value = 0;
                    app.condCheckBoxes(iCB).Text = num2str(iCB-1);
                    app.condCheckBoxes(iCB).ValueChangedFcn = createCallbackFcn(app,@cleanCheckboxes,true);

                    mod_ = mod(iCB,4);
                    if mod_ == 0
                        yOffset = yOffset - 25;
                        xOffset = 14;
                    else
                        xOffset = xOffset+40;
                    end
                end
                
                % Checking the number of WLs, if more than two, allow the
                % user to select two of them              
                if app.nlambda > 2
                    app.wlLabel.Visible = 'on';
                    xOffset = 54;
                    yOffset = 185;
                    for iCB=1:app.nlambda
                        app.wlCheckBoxes(iCB) = uicheckbox(app.NIRSPlotGUIUIFigure);
                        app.wlCheckBoxes(iCB).Position =[xOffset,yOffset,...
                            45, 15];
                        app.wlCheckBoxes(iCB).Value = 0;
                        app.wlCheckBoxes(iCB).Text = num2str(app.rawDotNirs.SD.Lambda(iCB));
                        mod_ = mod(iCB,4);
                        if mod_ == 0
                            yOffset = yOffset - 25;
                            xOffset = 14;
                        else
                            xOffset = xOffset+45;
                        end
                    end
                else
                     app.wlLabel.Visible = 'off';
                     app.wlCheckBoxes.delete;
                end
            else
                app.PlotButton.Enable = 0;
            end
        end
        
        % Button pushed function: LoadButton
        function LoadDotNirsDirPushed(app, event)
            dotNirsPath_tmp = uigetdir(app.dotNirsPath,'Select .nirs files folder');
            drawnow;
            figure(app.NIRSPlotGUIUIFigure);
            if dotNirsPath_tmp ~=0    
                app.dotNirsPath = dotNirsPath_tmp;
                LoadDotNirsDir(app, app.dotNirsPath);
            end
        end
        
        % Button pushed function: LoadButton
        function LoadDotNirsDir(app, dotNirsPath)
            app.dotNirsPath = dotNirsPath;
            %Search for .nirs files in the folder
            dotNirsFound = dir([app.dotNirsPath,filesep,'*.*nir*']);
            cd(app.dotNirsPath);
            % Add nodes to the tree with those .nirs files found
            if ~isempty(dotNirsFound)
                folder = split(dotNirsFound(1,1).folder,filesep);
                folder = folder{end};
                app.treeDotNirs.Children(1).Children.delete;
                app.treeDotNirs.Children(1).Text = folder;
                
                for i=1:length(dotNirsFound)
                    tn = uitreenode(app.treeDotNirs.Children(1),'Text',dotNirsFound(i).name);
                    if strcmp(app.dotNirsFile,dotNirsFound(i).name)==1
                        app.treeDotNirs.SelectedNodes = tn;
                    end
                end
                expand(app.treeDotNirs);
            else
                opts = struct('WindowStyle','modal',...
                    'Interpreter','none');
                errordlg('No .nirs files were found, try other folders.','Invalid dir',opts);
            end
            
        end
        
        
        % Callback function
        function goQuality(app, event)
            if isfield(app,'reportTable')
                app = rmfield(app,'reportTable');
            else
                app.reportTable = struct([]);
            end
            if app.rb1.Value ==1
                PlotButtonIndiv(app, event);
                generateQReport(app.reportTable);
            else
                PlotButtonGroup(app, event);
                generateQReport(app.reportTable);
            end
        end
        
        function PlotButtonIndiv(app, event)
            app.bpFmin = app.FreqMinEditField.Value;
            app.bpFmax = app.FreqMaxEditField.Value;
            app.windowSec = app.LengthsecSpinner.Value;
            app.windowOverlap = app.WindowsOverlapCtrl.Value*0.5;
            
            app.quality_threshold = app.QualityThresholdField.Value/100;
            app.sci_threshold = app.sciThresholdField.Value;
            app.psp_threshold = app.pspThresholdField.Value;
            app.condMask = zeros(1,app.nCond);
            if app.condCheckBoxes(1).Value == 1
                app.condMask = 'resting';
            else
                for iCB=1:app.nCond
                    app.condMask(iCB) = logical(app.condCheckBoxes(iCB+1).Value);
                end
            end
            
            if isvalid(app.wlCheckBoxes)
               for iCB=1:app.nlambda
                  app.lambdaMask(iCB) = app.wlCheckBoxes(iCB).Value; 
               end
            else
                app.lambdaMask = [1 1];
            end
            
            app.reportTable = qtnirs([app.dotNirsPath,filesep,app.dotNirsFile],...
                'freqCut',[app.bpFmin, app.bpFmax],...
                'window',app.windowSec,...
                'overlap',app.windowOverlap,....
                'qualityThreshold',app.quality_threshold,...
                'sciThreshold',app.sci_threshold,...
                'pspThreshold',app.psp_threshold,...
                'conditionsMask',app.condMask,...
                'dodFlag',app.dodCheckBox.Value,...
                'guiFlag',1,...
                'lambdaMask',app.lambdaMask);
        end
        
        
        function PlotButtonGroup(app, event)
            %Search for .nirs files in the folder
            dotNirsFound = dir([app.dotNirsPath,filesep,'*.*nir*']);
            
            app.bpFmin = app.FreqMinEditField.Value;
            app.bpFmax = app.FreqMaxEditField.Value;
            app.windowSec = app.LengthsecSpinner.Value;
            app.windowOverlap = app.WindowsOverlapCtrl.Value*0.5;
            
            app.quality_threshold = app.QualityThresholdField.Value/100;
            app.sci_threshold = app.sciThresholdField.Value;
            app.psp_threshold = app.pspThresholdField.Value;
            app.condMask = zeros(1,app.nCond);
            if app.condCheckBoxes(1).Value == 1
                app.condMask = 'resting';
            else
                for iCB=1:app.nCond
                    app.condMask(iCB) = logical(app.condCheckBoxes(iCB+1).Value);
                end
            end
            
            if isvalid(app.wlCheckBoxes)
               for iCB=1:app.nlambda
                  app.lambdaMask(iCB) = app.wlCheckBoxes(iCB).Value; 
               end
            else
                app.lambdaMask = [1 1];
            end
            
            for i = 1:numel(dotNirsFound)
                if((strcmp(dotNirsFound(i).name(end-4:end),'.nirs')==1) || ...
                        (strcmp(dotNirsFound(i).name(end-5:end),'.snirf')==1))
                    app.reportTable = [app.reportTable;qtnirs([app.dotNirsPath,filesep,dotNirsFound(i).name],...
                        'freqCut',[app.bpFmin, app.bpFmax],...
                        'window',app.windowSec,...
                        'overlap',app.windowOverlap,....
                        'qualityThreshold',app.quality_threshold,...
                        'sciThreshold',app.sci_threshold,...
                        'pspThreshold',app.psp_threshold,...
                        'conditionsMask',app.condMask,...
                        'dodFlag',app.dodCheckBox.Value,...
                        'guiFlag',0,...
                        'lambdaMask',app.lambdaMask)];
                    %app.reportTable(i) = qMats;
                    %clear qMats;
                    fprintf('Scan %i of %i processed.\n',i,numel(dotNirsFound));
                else
                    error('unsupported type');
                end                
            end            
        end
        
        function cleanCheckboxes(app,event)
            if strcmp(event.Source.Text,'Resting')==1 && event.Value==1
                for i =2:numel(app.condCheckBoxes)
                    app.condCheckBoxes(i).Value = 0;
                end
            elseif strcmp(event.Source.Text,'Resting')==0 && event.Value==1
                app.condCheckBoxes(1).Value = 0;
            end
        end
    end
    
    % Component initialization
    methods (Access = private)
        
        % Create UIFigure and components
        function createComponents(app)
            w.height = 700;
            w.width = 234;
            w.x = 150;
            w.y = 50;
            
            uicomp = zeros(0,4);
            % Create NIRSPlotGUIUIFigure and hide until all components are created
            app.NIRSPlotGUIUIFigure = uifigure('Visible', 'off','Tag','nirsplotGUI');
            app.NIRSPlotGUIUIFigure.Position = [w.x w.y w.width w.height];
            app.NIRSPlotGUIUIFigure.Name = 'QT-NIRS GUI';
            
            % Create NIRSPlotLabel
            uicomp(1,:) = [80 w.height-25 83 22];
            app.NIRSPlotLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.NIRSPlotLabel.FontSize = 18;
            app.NIRSPlotLabel.FontWeight = 'bold';
            app.NIRSPlotLabel.FontAngle = 'italic';
            app.NIRSPlotLabel.FontColor = [0 0 1];
            app.NIRSPlotLabel.Position = uicomp(1,:);
            app.NIRSPlotLabel.HorizontalAlignment = 'center';
            app.NIRSPlotLabel.Text = 'QT-NIRS';
            
            % Create nirsfileEditFieldLabel
            uicomp(2,:) = [24 530 60 uicomp(1,4)];
            app.nirsfileEditFieldLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.nirsfileEditFieldLabel.HorizontalAlignment = 'right';
            app.nirsfileEditFieldLabel.Position = uicomp(2,:);
            app.nirsfileEditFieldLabel.Text = '.nirs .snirf';
            
            % Create tree
            uicomp(2,:) = [24 550 200 120];
            app.treeDotNirs = uitree(app.NIRSPlotGUIUIFigure);
            app.treeDotNirs.Position = uicomp(2,:);
            app.treeDotNirs.SelectionChangedFcn = createCallbackFcn(app,@treeSelection,true);
            uitreenode(app.treeDotNirs,'Text','');           
            
            % Create LoadButton
            app.LoadButton = uibutton(app.NIRSPlotGUIUIFigure, 'push');            
            app.LoadButton.ButtonPushedFcn = createCallbackFcn(app,@LoadDotNirsDirPushed,true);
            app.LoadButton.Position = [120 515 79 22];
            app.LoadButton.Text = 'Sel. Folder';
            
            % Create dodCheckBox
            app.dodCheckBox = uicheckbox(app.NIRSPlotGUIUIFigure);
            app.dodCheckBox.Value = 0;
            app.dodCheckBox.Text = 'Load OD';
            app.dodCheckBox.Position = [24 510 70 22];
            app.dodCheckBox.Enable = 'off';
            
            % Create SourcesLabel
            app.SourcesLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.SourcesLabel.Position = [24 485 50 22];
            app.SourcesLabel.Text = 'Sources';
            
            % Create DetectorsLabel
            app.DetectorsLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.DetectorsLabel.Position = [24 465 57 22];
            app.DetectorsLabel.Text = 'Detectors';
            
            % Create ChannelsLabel
            app.ChannelsLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.ChannelsLabel.Position = [125 485 56 22];
            app.ChannelsLabel.Text = 'Channels';
            
            % Create DurLabel
            app.DurLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.DurLabel.Position = [125 465 42 22];
            app.DurLabel.Text = 'Time(s)';
            
            % Create WLlabel
            app.wlLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.wlLabel.Position = [24 185 42 22];
            app.wlLabel.Text = 'WL';
            app.wlLabel.Visible = 'off';
            
            % Create nSrcLab
            app.nSrcLab = uilabel(app.NIRSPlotGUIUIFigure);
            app.nSrcLab.Position = [94 485 37 22];
            app.nSrcLab.Text = '0';
            
            % Create nDetecLab
            app.nDetecLab = uilabel(app.NIRSPlotGUIUIFigure);
            app.nDetecLab.Position = [94 465 40 22];
            app.nDetecLab.Text = '0';
            
            % Create nChanLab
            app.nChanLab = uilabel(app.NIRSPlotGUIUIFigure);
            app.nChanLab.Position = [195 485 40 22];
            app.nChanLab.Text = '0';
            
            % Create secDurLab
            app.secDurLab = uilabel(app.NIRSPlotGUIUIFigure);
            app.secDurLab.Position = [185 465 80 22];
            app.secDurLab.Text = '0';
            
            %----
            % from 465 to 370
            % central label
            % Create scipspLabel
            app.scipspLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.scipspLabel.Position = [69 435 97 22];
            app.scipspLabel.Text = 'Metric thresholds';
            % Sci label and field
            % Create sciThresholdLabel
            app.sciThresholdLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.sciThresholdLabel.HorizontalAlignment = 'left';
            app.sciThresholdLabel.Position = [24 405 30 22];
            app.sciThresholdLabel.Text = 'SCI';
            
            % Create sciThresholdField
            app.sciThresholdField = uieditfield(app.NIRSPlotGUIUIFigure, 'numeric');
            app.sciThresholdField.Position = [54 405 40 22];
            app.sciThresholdField.Value = 0.8;
            % PSP label and field
            % Create pspThresholdLabel
            app.pspThresholdLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.pspThresholdLabel.HorizontalAlignment = 'left';
            app.pspThresholdLabel.Position = [124 405 30 22];
            app.pspThresholdLabel.Text = 'PSP';
            
            % Create pspThresholdField
            app.pspThresholdField = uieditfield(app.NIRSPlotGUIUIFigure, 'numeric');
            app.pspThresholdField.Position = [155 405 40 22];
            app.pspThresholdField.Value = 0.1;
            
                      
            % Create CutoffFrequencyLabel
            app.CutoffFrequencyLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.CutoffFrequencyLabel.Position = [69 370 97 22];
            app.CutoffFrequencyLabel.Text = 'Cutoff Frequency';
            
            % Create FreqMinEditFieldLabel
            app.FreqMinEditFieldLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.FreqMinEditFieldLabel.HorizontalAlignment = 'left';
            app.FreqMinEditFieldLabel.Position = [24 340 30 22];
            app.FreqMinEditFieldLabel.Text = 'Min.';
            
            % Create FreqMinEditField
            app.FreqMinEditField = uieditfield(app.NIRSPlotGUIUIFigure, 'numeric');
            app.FreqMinEditField.Position = [54 340 40 22];
            app.FreqMinEditField.Value = 0.5;
            
            % Create FreqMaxEditFieldLabel
            app.FreqMaxEditFieldLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.FreqMaxEditFieldLabel.HorizontalAlignment = 'left';
            app.FreqMaxEditFieldLabel.Position = [124 340 30 22];
            app.FreqMaxEditFieldLabel.Text = 'Max.';
            
            % Create FreqMaxEditField
            app.FreqMaxEditField = uieditfield(app.NIRSPlotGUIUIFigure, 'numeric');
            app.FreqMaxEditField.Position = [155 340 40 22];
            app.FreqMaxEditField.Value = 2.5;
            
            % Create WindowLabel
            app.WindowLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.WindowLabel.Position = [94 310 48 22];
            app.WindowLabel.Text = 'Window';
            
            % Create LengthsecSpinnerLabel
            app.LengthsecSpinnerLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.LengthsecSpinnerLabel.Position = [24 285 72 22];
            app.LengthsecSpinnerLabel.Text = 'Length (sec)';
            
            % Create LengthsecSpinner
            app.LengthsecSpinner = uispinner(app.NIRSPlotGUIUIFigure);
            app.LengthsecSpinner.Position = [107 285 59 22];
            app.LengthsecSpinner.Value = 5;
            
            % Create OverlappingLabel
            app.OverlappingLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.OverlappingLabel.Position = [24 250 70 28];
            app.OverlappingLabel.Text = {'Overlapping'};
            
 
            % Create dodWindowsOverlap
            app.WindowsOverlapCtrl = uicheckbox(app.NIRSPlotGUIUIFigure);
            app.WindowsOverlapCtrl.Value = 0;
            app.WindowsOverlapCtrl.Text = '';
            app.WindowsOverlapCtrl.Position = [99 250 70 22];
            app.WindowsOverlapCtrl.Enable = 'on';
            
            % Create QualityThresholdLabel
            app.QualityThresholdLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.QualityThresholdLabel.HorizontalAlignment = 'left';
            app.QualityThresholdLabel.Position = [24 215 80 28];
            app.QualityThresholdLabel.Text = {'Quality'; '(0-100)'};
            
            % Create QualityThresholdField
            app.QualityThresholdField = uieditfield(app.NIRSPlotGUIUIFigure, 'numeric');
            app.QualityThresholdField.Position = [99 215 50 22];
            app.QualityThresholdField.Value = 90;
            
            % Create stimCheckBLabel
            app.condCheckBLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.condCheckBLabel.HorizontalAlignment = 'left';
            app.condCheckBLabel.Position = [90, 155, 120, 28];
            app.condCheckBLabel.Text = 'Conditions';
            
            % Create PlotButton
            app.PlotButton = uibutton(app.NIRSPlotGUIUIFigure, 'push');
            app.PlotButton.ButtonPushedFcn = createCallbackFcn(app, @goQuality, true);
            app.PlotButton.Position = [120 35 70 22];
            app.PlotButton.Enable = 0;
            app.PlotButton.Text = 'Go';
            
            % Create radiobuttons for Subject or Group level analysis
            app.bg = uibuttongroup(app.NIRSPlotGUIUIFigure,'Position',[15 32 90 45]);   
            app.bg.BorderType = 'none';
            app.rb1 = uiradiobutton(app.bg,'Position',[10 25 91 15]);
            app.rb2 = uiradiobutton(app.bg,'Position',[10 7 91 15]);
            app.rb1.Text = 'Scan';
            app.rb2.Text = 'Group';
            app.rb1.Value = true;
                       
            % Show the figure after all components are created
            app.NIRSPlotGUIUIFigure.Visible = 'on';
        end
    end
    
    % App creation and deletion
    methods (Access = public)
        
        % Construct app
        function app = qtnirsLoadFileGUI(varargin)
            
            % Create UIFigure and components
            createComponents(app)
            
            % Register the app with App Designer
            registerApp(app, app.NIRSPlotGUIUIFigure)
            app.dotNirsPath = pwd;
            % Execute the startup function
            if ~isempty(varargin)
                runStartupFcn(app, @(app)startupFcn(app, varargin{:}))
            end
            
            if nargout == 0
                clear app
            end
        end
        
        % Code that executes before app deletion
        function delete(app)
            
            % Delete UIFigure when app is deleted
            delete(app.NIRSPlotGUIUIFigure)
        end
    end
end
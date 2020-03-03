classdef nirsplotLoadFileGUI < matlab.apps.AppBase
    
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
        WindowsoverlapSlider    matlab.ui.control.Slider
        QualityThresholdLabel   matlab.ui.control.Label
        QualityThresholdField   matlab.ui.control.NumericEditField
        condCheckBoxes          matlab.ui.control.CheckBox
        condCheckBLabel         matlab.ui.control.Label
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
        condMask % Binary mask for selecting the conditions of interest
        dotNirsFilePath % Path of the .nirs file to analyze
    end
    
    
    % Callbacks that handle component events
    methods (Access = private)
        
        function startupFcn(app, input_param)
            %unpacking inputs
            if isfield(input_param,'dotNirsPath')
                app.dotNirsFilePath = input_param.dotNirsPath;
                [PathFolder,FileName,ext] = fileparts(app.dotNirsFilePath);
                cd(PathFolder);
                app.nirsfileEditField.Value = FileName;       
            
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
                app.WindowsoverlapSlider.Value = app.windowOverlap*100;
            end
            if isfield(input_param,'quality_threshold')
                app.quality_threshold = input_param.quality_threshold;
                app.QualityThresholdField.Value = app.quality_threshold;
            end
            if isfield(input_param,'cond_mask')
                % Create and populate StimuliSelectionBoxes
                xOffset = 16;
                yOffset = 145;
                app.nCond = length(input_param.cond_mask);
                app.condMask = input_param.cond_mask;
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
            end
        end
        
        % Button pushed function: LoadButton
        function LoadButtonPushed(app, event)
            [FileName,PathFolder,~] = uigetfile('*.nirs','Please select the .nirs file to import');
            if ~isequal(FileName,0)
                app.dotNirsFilePath = [PathFolder FileName];
                cd(PathFolder);
                app.rawDotNirs = load(app.dotNirsFilePath,'-mat');
                app.nirsfileEditField.Value = FileName;
                
                [app.sampDur,app.nChannels] = size(app.rawDotNirs.d);
                app.nChannels = app.nChannels/2;
                app.nSources        = size(app.rawDotNirs.SD.SrcPos,1);
                app.nDetectors      = size(app.rawDotNirs.SD.DetPos,1);
                app.secDur          = app.rawDotNirs.t(end);
                app.nSrcLab.Text    = num2str(app.nSources);
                app.nDetecLab.Text  = num2str(app.nDetectors);
                app.nChanLab.Text   = num2str(app.nChannels);
                app.secDurLab.Text  = num2str(app.secDur,'%.2f');
                
                % Create StimuliSelectionBoxes
                xOffset = 16;
                yOffset = 145;
                app.nCond = size(app.rawDotNirs.s,2);
                app.condCheckBoxes.delete;
                for iCB=1:app.nCond
                    app.condCheckBoxes(iCB) = uicheckbox(app.NIRSPlotGUIUIFigure);
                    app.condCheckBoxes(iCB).Position =[xOffset,yOffset,...
                        35, 15];
                    app.condCheckBoxes(iCB).Value = 1;
                    app.condCheckBoxes(iCB).Text = num2str(iCB);
                    mod_ = mod(iCB,5);
                    if mod_ == 0
                        yOffset = yOffset - 25;
                        xOffset = 16;
                    else
                        xOffset = xOffset+40;
                    end
                end
                
            end
        end
        
        % Callback function
        function PlotButtonPushed2(app, event)
            app.bpFmin = app.FreqMinEditField.Value;
            app.bpFmax = app.FreqMaxEditField.Value;
            app.windowSec = app.LengthsecSpinner.Value;
            app.windowOverlap = (app.WindowsoverlapSlider.Value)/100;
            app.quality_threshold = app.QualityThresholdField.Value;
            app.condMask = zeros(1,app.nCond);
            for iCB=1:app.nCond
                app.condMask(iCB) = logical(app.condCheckBoxes(iCB).Value);
            end
            
            setappdata(app.NIRSPlotGUIUIFigure,'dotNirs_filename',app.nirsfileEditField.Value);
            
            %app.reportTable = nirsplot(app.rawDotNirs, [app.bpFmin, app.bpFmax],...
            %    app.windowSec,app.windowOverlap,app.quality_threshold,...
            %    app.condMask);
            app.reportTable = nirsplot(app.dotNirsFilePath,...
                'freqCut',[app.bpFmin, app.bpFmax],...
                'window',app.windowSec,...
                'overlap',app.windowOverlap,....
                'qualityThreshold',app.quality_threshold,...
                'conditionsMask',app.condMask);
        end
    end
    
    % Component initialization
    methods (Access = private)
        
        % Create UIFigure and components
        function createComponents(app)
            
            % Create NIRSPlotGUIUIFigure and hide until all components are created
            app.NIRSPlotGUIUIFigure = uifigure('Visible', 'off','Tag','nirsplotGUI');
            app.NIRSPlotGUIUIFigure.Position = [100 100 234 655];
            app.NIRSPlotGUIUIFigure.Name = 'NIRSPlot GUI';
            
            % Create NIRSPlotLabel
            app.NIRSPlotLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.NIRSPlotLabel.FontSize = 18;
            app.NIRSPlotLabel.FontWeight = 'bold';
            app.NIRSPlotLabel.FontAngle = 'italic';
            app.NIRSPlotLabel.FontColor = [0 0 1];
            app.NIRSPlotLabel.Position = [85 598 83 22];
            app.NIRSPlotLabel.Text = 'NIRSPlot';
            
            
            % Create nirsfileEditFieldLabel
            app.nirsfileEditFieldLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.nirsfileEditFieldLabel.HorizontalAlignment = 'right';
            app.nirsfileEditFieldLabel.Position = [24 545 47 22];
            app.nirsfileEditFieldLabel.Text = '.nirs file';
            
            % Create nirsfileEditField
            app.nirsfileEditField = uieditfield(app.NIRSPlotGUIUIFigure, 'text');
            app.nirsfileEditField.Position = [91 545 100 22];
            
            % Create LoadButton
            app.LoadButton = uibutton(app.NIRSPlotGUIUIFigure, 'push');
            app.LoadButton.ButtonPushedFcn = createCallbackFcn(app, @LoadButtonPushed, true);
            app.LoadButton.Position = [120 515 79 22];
            app.LoadButton.Text = 'Load';
            
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
            app.ChannelsLabel.Position = [24 445 56 22];
            app.ChannelsLabel.Text = 'Channels';
            
            % Create DurLabel
            app.DurLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.DurLabel.Position = [24 425 32 22];
            app.DurLabel.Text = 'Time';
            
            % Create nSrcLab
            app.nSrcLab = uilabel(app.NIRSPlotGUIUIFigure);
            app.nSrcLab.Position = [109 485 37 22];
            app.nSrcLab.Text = '0';
            
            % Create nDetecLab
            app.nDetecLab = uilabel(app.NIRSPlotGUIUIFigure);
            app.nDetecLab.Position = [109 465 40 22];
            app.nDetecLab.Text = '0';
            
            % Create nChanLab
            app.nChanLab = uilabel(app.NIRSPlotGUIUIFigure);
            app.nChanLab.Position = [109 445 40 22];
            app.nChanLab.Text = '0';
            
            % Create secDurLab
            app.secDurLab = uilabel(app.NIRSPlotGUIUIFigure);
            app.secDurLab.Position = [109 425 80 22];
            app.secDurLab.Text = '0';
            
            % Create CutoffFrequencyLabel
            app.CutoffFrequencyLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.CutoffFrequencyLabel.Position = [69 395 97 22];
            app.CutoffFrequencyLabel.Text = 'Cutoff Frequency';
            
            % Create FreqMinEditFieldLabel
            app.FreqMinEditFieldLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.FreqMinEditFieldLabel.HorizontalAlignment = 'right';
            app.FreqMinEditFieldLabel.Position = [24 365 60 22];
            app.FreqMinEditFieldLabel.Text = 'Freq. Min.';
            
            % Create FreqMinEditField
            app.FreqMinEditField = uieditfield(app.NIRSPlotGUIUIFigure, 'numeric');
            app.FreqMinEditField.Position = [99 365 50 22];
            app.FreqMinEditField.Value = 0.5;
            
            % Create FreqMaxEditFieldLabel
            app.FreqMaxEditFieldLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.FreqMaxEditFieldLabel.HorizontalAlignment = 'right';
            app.FreqMaxEditFieldLabel.Position = [24 340 63 22];
            app.FreqMaxEditFieldLabel.Text = 'Freq. Max.';
            
            % Create FreqMaxEditField
            app.FreqMaxEditField = uieditfield(app.NIRSPlotGUIUIFigure, 'numeric');
            app.FreqMaxEditField.Position = [99 340 50 22];
            app.FreqMaxEditField.Value = 2.5;
            
            % Create WindowLabel
            app.WindowLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.WindowLabel.Position = [94 305 48 22];
            app.WindowLabel.Text = 'Window';
            
            % Create LengthsecSpinnerLabel
            app.LengthsecSpinnerLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.LengthsecSpinnerLabel.Position = [24 280 72 22];
            app.LengthsecSpinnerLabel.Text = 'Length (sec)';
            
            % Create LengthsecSpinner
            app.LengthsecSpinner = uispinner(app.NIRSPlotGUIUIFigure);
            app.LengthsecSpinner.Position = [107 280 59 22];
            app.LengthsecSpinner.Value = 5;
            
            % Create OverlappingLabel
            app.OverlappingLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.OverlappingLabel.Position = [18 240 70 28];
            app.OverlappingLabel.Text = {'Overlapping'; '(%)'};
            
            % Create WindowsoverlapSlider
            app.WindowsoverlapSlider = uislider(app.NIRSPlotGUIUIFigure);
            app.WindowsoverlapSlider.MajorTicks = [0 25 50 75 100];
            app.WindowsoverlapSlider.MajorTickLabels = {'0', '25', '50', '75', '100'};
            app.WindowsoverlapSlider.Position = [99 255 94 3];
            
            % Create QualityThresholdLabel
            app.QualityThresholdLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.QualityThresholdLabel.HorizontalAlignment = 'left';
            app.QualityThresholdLabel.Position = [16 195 80 28];
            app.QualityThresholdLabel.Text = {'Quality'; '(0.1-1.0)'};
            
            % Create QualityThresholdField
            app.QualityThresholdField = uieditfield(app.NIRSPlotGUIUIFigure, 'numeric');
            app.QualityThresholdField.Position = [99 195 50 22];
            app.QualityThresholdField.Value = 0.9;
            
            % Create stimCheckBLabel
            app.condCheckBLabel = uilabel(app.NIRSPlotGUIUIFigure);
            app.condCheckBLabel.HorizontalAlignment = 'left';
            app.condCheckBLabel.Position = [90, 160, 120, 28];
            app.condCheckBLabel.Text = 'Conditions';
            
            % Create PlotButton
            app.PlotButton = uibutton(app.NIRSPlotGUIUIFigure, 'push');
            app.PlotButton.ButtonPushedFcn = createCallbackFcn(app, @PlotButtonPushed2, true);
            app.PlotButton.Position = [120 35 70 22];
            app.PlotButton.Text = 'Plot';
            
            % Show the figure after all components are created
            app.NIRSPlotGUIUIFigure.Visible = 'on';
        end
    end
    
    % App creation and deletion
    methods (Access = public)
        
        % Construct app
        function app = nirsplotLoadFileGUI(varargin)
            
            % Create UIFigure and components
            createComponents(app)
            
            % Register the app with App Designer
            registerApp(app, app.NIRSPlotGUIUIFigure)
            
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
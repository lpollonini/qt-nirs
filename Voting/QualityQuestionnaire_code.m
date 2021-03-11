classdef QualityQuestionnaire_code < matlab.apps.AppBase
    
    % Properties that correspond to app components
    properties (Access = public)
        UIFigure            matlab.ui.Figure
        VisualizationPanel  matlab.ui.container.Panel
        GridLayout          matlab.ui.container.GridLayout
        IntAxes             matlab.ui.control.UIAxes
        HbAxes              matlab.ui.control.UIAxes
        AssessmentPanel     matlab.ui.container.Panel
        GridLayout2         matlab.ui.container.GridLayout
        LowButton           matlab.ui.control.Button
        MediumButton        matlab.ui.control.Button
        HighButton          matlab.ui.control.Button
        Panel_2             matlab.ui.container.Panel
        treeDotNirs         matlab.ui.container.Tree
        ScansNode           matlab.ui.container.TreeNode
        Node2               matlab.ui.container.TreeNode
        Node3               matlab.ui.container.TreeNode
        Node4               matlab.ui.container.TreeNode
        LoadscanButton      matlab.ui.control.Button
        Panel               matlab.ui.container.Panel
        Button_2            matlab.ui.control.Button
        Button_3            matlab.ui.control.Button
        scanIDLab           matlab.ui.control.Label
    end
    
    properties (Access = public)
        rawDotNirs  % This is the content from .nirs files loaded as '-mat' format.
        nSources    % Number of sources
        nDetectors  % Number of detectors
        nChannels   % Number of channels
        nCond       % Number of conditions
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
        flagdod     % flag to detect OD data
        flagdc      % flag to detect Hb data
        ptrOffsetSamp     % pointer to current offset sample
        ptrOnsetSamp      % pointer to current onset sample
        ptrChannel  	  %pointer to current channel
        ptrWindow 
        ptrSubject        % pointer to the current subject
        Fs          	  % Sampling frequency
        overlap_samples
        window_samples
        n_windows
        n_samples
        scoresFile
        scoresMat
        nChannelsFx
        nSubjFx
        flagScFl
    end
    % Component initialization
    methods (Access = private)
        
        % Create UIFigure and components
        function createComponents(app)
            
            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 730 697];
            app.UIFigure.Name = 'UI Figure';
            
            % Create VisualizationPanel
            app.VisualizationPanel = uipanel(app.UIFigure);
            app.VisualizationPanel.Position = [17 265 700 411];
            
            % Create GridLayout
            app.GridLayout = uigridlayout(app.VisualizationPanel);
            app.GridLayout.ColumnWidth = {'1x'};
            
            % Create IntAxes
            app.IntAxes = uiaxes(app.GridLayout);
            title(app.IntAxes, 'Raw');
            xlabel(app.IntAxes, 'Time (s)');
            ylabel(app.IntAxes, 'Intensity');
            app.IntAxes.Layout.Row = 1;
            app.IntAxes.Layout.Column = 1;  
            
            % Create HbAxes
            app.HbAxes = uiaxes(app.GridLayout);
            title(app.HbAxes, 'Hb');
            xlabel(app.HbAxes, 'Time (s)');
            ylabel(app.HbAxes, '\DeltaHb (\muMol)');           
            app.HbAxes.Layout.Row = 2;
            app.HbAxes.Layout.Column = 1;
            
            % Create AssessmentPanel
            app.AssessmentPanel = uipanel(app.UIFigure);
            app.AssessmentPanel.Position = [18 16 698 236];
            
            % Create GridLayout2
            app.GridLayout2 = uigridlayout(app.AssessmentPanel);
            app.GridLayout2.ColumnWidth = {'1x', '3x'};
            app.GridLayout2.RowHeight = {'1x'};
            
            % Create Panel
            app.Panel = uipanel(app.GridLayout2);
            app.Panel.Layout.Row = 1;
            app.Panel.Layout.Column = 2;
            
            % Create Button_2
            app.Button_2 = uibutton(app.Panel, 'push');
            app.Button_2.Position = [15 13 74 55];
            app.Button_2.Text = '<<';
            app.Button_2.Enable = 'off';
            app.Button_2.ButtonPushedFcn = createCallbackFcn(app,@go2window,true);
            
            % Create Button_3
            app.Button_3 = uibutton(app.Panel, 'push');
            app.Button_3.Position = [400 13 74 55];
            app.Button_3.Text = '>>';
            app.Button_3.Enable = 'off';
            app.Button_3.ButtonPushedFcn = createCallbackFcn(app,@go2window,true);
            
            % Create LowButton
            app.LowButton = uibutton(app.Panel, 'push');
            app.LowButton.FontSize = 16;
            %app.LowButton.FontWeight = 'bold';
            app.LowButton.FontColor = [1 0 0];
            app.LowButton.Position = [77 102 78 70];
            app.LowButton.Text = 'Bad';
            app.LowButton.Enable = 'off';
            app.LowButton.ButtonPushedFcn = createCallbackFcn(app,@go2window,true);
            
%             % Create MediumButton
%             app.MediumButton = uibutton(app.Panel, 'push');
%             app.MediumButton.FontSize = 16;
%             app.MediumButton.FontWeight = 'bold';
%             app.MediumButton.FontColor = [0.9294 0.6941 0.1255];
%             app.MediumButton.Position = [206 102 78 70];
%             app.MediumButton.Text = 'Medium';
            
            % Create HighButton
            app.HighButton = uibutton(app.Panel, 'push');
            app.HighButton.FontSize = 16;
            %app.HighButton.FontWeight = 'bold';
            app.HighButton.FontColor = [0.3922 0.8314 0.0745];
            app.HighButton.Position = [334 102 78 70];
            app.HighButton.Text = 'Good';
            app.HighButton.Enable = 'off';
            app.HighButton.ButtonPushedFcn = createCallbackFcn(app,@go2window,true);
            
            % Create Panel_2
            app.Panel_2 = uipanel(app.GridLayout2);
            app.Panel_2.Layout.Row = 1;
            app.Panel_2.Layout.Column = 1;
            
            % Create Tree
            app.treeDotNirs = uitree(app.Panel_2);
            app.treeDotNirs.Position = [1 60 166 135];
            uitreenode(app.treeDotNirs,'Text','');
            
            % Create LoadscanButton
            app.LoadscanButton = uibutton(app.Panel_2, 'push');
            app.LoadscanButton.Position = [59 30 100 22];
            app.LoadscanButton.Text = 'Load scan';
            app.LoadscanButton.ButtonPushedFcn = createCallbackFcn(app,@go2scan,true);
            
            % Create scanID label
            app.scanIDLab = uilabel(app.Panel_2);
            app.scanIDLab.Position = [10 5 150 22];
            app.scanIDLab.Text = 'Scan:';
            app.scanIDLab.FontSize = 14;

            
            
            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
        
        % Button pushed function: LoadButton
        function LoadDotNirsDirPushed(app)
%             app.dotNirsPath = uigetdir(app.dotNirsPath,'Select .nirs files folder');
%             drawnow;
%             figure(app.UIFigure);
            %app.dotNirsPath = 'C:\Users\smonter2\OneDrive - University Of Houston\Nirsplot-OneDrive-SAMH\Data\CPGM\dotNirs\modified';
            if app.dotNirsPath~=0
                LoadDotNirsDir(app);
            end
        end
        
        
        function LoadDotNirsDir(app)
            %Search for .nirs files in the folder
            dotNirsFound = dir([app.dotNirsPath,filesep,'*.*nir*']);            
            % Add nodes to the tree with those .nirs files found
            if ~isempty(dotNirsFound)
                folder = split(dotNirsFound(1,1).folder,filesep);
                folder = folder{end};
                app.treeDotNirs.Children(1).Children.delete;
                app.treeDotNirs.Children(1).Text = folder;
                app.nSubjFx = length(dotNirsFound);
                for i=1:app.nSubjFx
                    tn = uitreenode(app.treeDotNirs.Children(1),'Text',dotNirsFound(i).name);
                    %if strcmp(app.dotNirsFile,dotNirsFound(i).name)==1
                    %    app.treeDotNirs.SelectedNodes = tn;
                    % end
                    if app.flagScFl==0
                        raw = load([app.dotNirsPath,filesep,dotNirsFound(i).name],'-mat');
                        nChannels = size(raw.d,2)/2;
                        Fs = 1/mean(diff(raw.t));
                        n_samples = length(raw.t);
                        % overlap_samples = ceil(app.windowSec*Fs*app.windowOverlap);
                        window_samples = floor(app.windowSec*Fs);
                        if app.windowOverlap ==0
                            nWindows = floor((n_samples)/(window_samples));
                        else % Valid for overlap=50%
                            nWindows = 2*floor((n_samples)/(window_samples))-1;
                        end
                        app.scoresMat{i}=zeros(nChannels,nWindows);
                    end
                end
                expand(app.treeDotNirs);
            else
                opts = struct('WindowStyle','modal',...
                    'Interpreter','none');
                errordlg('No .nirs files were found, try other folders.','Invalid dir',opts);
            end
            
        end
        
        
        function go2scan(app,event)
            if ~isempty(event)
                app.ptrWindow = 1;
                app.ptrChannel = 1;
            end
            if isempty(app.treeDotNirs.SelectedNodes.Children) %Is it a leaf (file)?
                app.scanIDLab.Text = 'Scan:';
                dotNirsFileSel = app.treeDotNirs.SelectedNodes.Text;
                app.dotNirsFile = dotNirsFileSel;
                disp(['Loading: ',app.dotNirsPath,filesep,dotNirsFileSel]);
                [filepath,name,ext] = fileparts([app.dotNirsPath,filesep,dotNirsFileSel]);
                app.scanIDLab.Text = ['Scan: ',name];
                name = strsplit(name,'_');
                name = name{1}(2:end);
                app.ptrSubject = str2num(name);
                switch ext
                    case '.nirs'
                        app.rawDotNirs = load([app.dotNirsPath,filesep,app.dotNirsFile],'-mat');
                    case '.snirf'
                        rawsnirf = SnirfClass([app.dotNirsPath,filesep,app.dotNirsFile]);
                        %app.rawDotNirs.d = rawsnirf.Get_d;
                        app.rawDotNirs.s = rawsnirf.Get_s;
                        app.rawDotNirs.t = rawsnirf.Get_t;
                        app.rawDotNirs.SD = rawsnirf.Get_SD;
                    otherwise
                        error('The input file should be a .nirs file format');
                end
                app.nlambda         = length(app.rawDotNirs.SD.Lambda);
                app.nChannels       = size(app.rawDotNirs.SD.MeasList,1)/app.nlambda;
                app.nSources        = size(app.rawDotNirs.SD.SrcPos,1);
                app.nDetectors      = size(app.rawDotNirs.SD.DetPos,1);
                app.secDur          = app.rawDotNirs.t(end);
                app.Fs  = 1/mean(diff(app.rawDotNirs.t));
                app.n_samples = length(app.rawDotNirs.t);
                %                 app.nSrcLab.Text    = num2str(app.nSources);
                %                 app.nDetecLab.Text  = num2str(app.nDetectors);
                %                 app.nChanLab.Text   = num2str(app.nChannels);
                %                 app.secDurLab.Text  = num2str(app.secDur,'%.2f');
                
                % Checking for OD data
                if isfield(app.rawDotNirs,'procResult')
                    if ~isempty(app.rawDotNirs.procResult.dod)
                        %app.dodCheckBox.Enable = 'on';
                        app.flagdod = true;
                    end
                    if ~isempty(app.rawDotNirs.procResult.dc)
                        %app.dodCheckBox.Enable = 'on';
                        app.flagdc = true;
                    end
                else
                    app.flagdod = false;
                    app.flagdc =false;
                end
                

                %Sort channels wl1-wl2-wl1-...
                [~, idx]=sortrows(app.rawDotNirs.SD.MeasList,[3 1 2]);
                app.rawDotNirs.SD.MeasList = app.rawDotNirs.SD.MeasList(idx,:);
                app.rawDotNirs.d = app.rawDotNirs.d(:,idx);
                app.rawDotNirs.SD.MeasListAct = app.rawDotNirs.SD.MeasListAct(idx);
                app.rawDotNirs.SD.MeasListVis = app.rawDotNirs.SD.MeasListVis(idx);
                app.overlap_samples = ceil(app.windowSec*app.Fs*app.windowOverlap);
                app.window_samples = floor(app.windowSec*app.Fs);
                if app.windowOverlap ==0
                    app.n_windows = floor((app.n_samples)/(app.window_samples));
                else % Valid for overlap=50%
                    app.n_windows = 2*floor((app.n_samples)/(app.window_samples))-1;
                end
                app.Button_2.Enable = 'on';
                app.Button_3.Enable = 'on';
                app.LowButton.Enable = 'on';
                app.HighButton.Enable = 'on';
                %app.ptrOffsetSamp = 0;
%                 if app.flagScFl~=1
%                     app.ptrOnsetSamp  = 1;
%                     app.ptrOffsetSamp = app.ptrOnsetSamp + app.window_samples;
%                 end
                app.go2window(app);
            end
        end
        
        function updateAxes(app,iChannel,xLimWindow)
            if  ~app.flagdod
                error('No Hb data in the selected .nirs file.');
            end
            t = app.rawDotNirs.t;
            %disp(['channel:' num2str(iChannel)]);
            raw = app.rawDotNirs;            
            disp(xLimWindow)
            % Intensity
            signal1 = raw.d(xLimWindow(1):xLimWindow(2),iChannel);
            signal2 = raw.d(xLimWindow(1):xLimWindow(2),iChannel+1);
            
            YLimStd = [min([signal1;signal2],[],'all'),...
                max([signal1; signal2],[],'all')].*[0.95,1.05];
            %XLimStd = myAxes.IntAxes.XLim;
            
            cla(app.IntAxes);
            
            plot(app.IntAxes,t(xLimWindow(1):xLimWindow(2)),...
                signal1,'-k');
            hold(app.IntAxes,'on');
            plot(app.IntAxes,t(xLimWindow(1):xLimWindow(2)),...
                signal2,'--k');
            %xlim(myAxes.IntAxes,XLimStd);
            app.IntAxes.XLim= [t(xLimWindow(1)),t(xLimWindow(2))];
            ylim(app.IntAxes,YLimStd);
            legend(app.IntAxes,'WL1','WL2');
            title(app.IntAxes, ['Raw (' num2str(app.ptrWindow) '/' num2str(app.n_windows) ')']);
            ylabel(app.IntAxes, ['Intensity (' num2str(app.ptrChannel) '/' num2str(app.nChannels) ')']);


                      
            % Hb
            signal1 = raw.procResult.dc(xLimWindow(1):xLimWindow(2),iChannel);
            signal2 = raw.procResult.dc(xLimWindow(1):xLimWindow(2),iChannel+1);
            YLimStd = [min([signal1;signal2],[],'all'),...
                max([signal1; signal2],[],'all')].*[0.95,1.05];
            XLimStd = app.HbAxes.XLim;
            cla(app.HbAxes);
            
            plot(app.HbAxes,t(xLimWindow(1):xLimWindow(2)),...
                signal1,'-r');
            hold(app.HbAxes,'on');
            plot(app.HbAxes,t(xLimWindow(1):xLimWindow(2)),...
                signal2,'-b');
            xlim(app.HbAxes,XLimStd);
            %ylim(myAxes.HbAxes,YLimStd);
            app.HbAxes.XLim= [t(xLimWindow(1)),t(xLimWindow(2))];
            legend(app.HbAxes,'HbO_2','HbR');
            title(app.HbAxes, ['Hb (' num2str(app.ptrWindow) '/' num2str(app.n_windows) ')']);
            ylabel(app.HbAxes, ['\DeltaHb (' num2str(app.ptrChannel) '/' num2str(app.nChannels) ')']);
        end
             
        
        function go2window(app,event)
            %             if ~isa(event,'matlab.ui.eventdata.ButtonPushedData')
            %                 app.ptrOnsetSamp  = 1;
            %                 app.ptrOffsetSamp = app.ptrOnsetSamp + app.window_samples;
            %             else
            if isa(event,'matlab.ui.eventdata.ButtonPushedData')
                if strcmp(event.Source.Text,'Bad')
                    disp('bad');
                    %app.ptrOnsetSamp = app.ptrOffsetSamp + 1 - app.overlap_samples;
                    %app.ptrOffsetSamp = app.ptrOnsetSamp + app.window_samples;
                    app.scoresMat{app.ptrSubject}(app.ptrChannel,app.ptrWindow) = -1;
                    app.ptrWindow = app.ptrWindow +1;
                end
                if strcmp(event.Source.Text,'Good')
                    disp('good');
                    %app.ptrOnsetSamp = app.ptrOffsetSamp + 1 - app.overlap_samples;
                    %app.ptrOffsetSamp = app.ptrOnsetSamp + app.window_samples;
                    app.scoresMat{app.ptrSubject}(app.ptrChannel,app.ptrWindow) = 1;
                    app.ptrWindow = app.ptrWindow +1;
                end
                if strcmp(event.Source.Text,'>>')
                    disp('>>');
                    %app.ptrOnsetSamp = app.ptrOffsetSamp + 1 - app.overlap_samples;
                    %app.ptrOffsetSamp = app.ptrOnsetSamp + app.window_samples;
                    app.ptrWindow = app.ptrWindow +1;
                end
                if strcmp(event.Source.Text,'<<')
                    disp('<<');
                    %app.ptrOffsetSamp = app.ptrOnsetSamp-1+app.overlap_samples;
                    %app.ptrOnsetSamp = app.ptrOnsetSamp - app.window_samples +app.overlap_samples - 1;
                    app.ptrWindow = app.ptrWindow -1;
                end
            end
            app.ptrOnsetSamp = (app.ptrWindow-1)*app.window_samples + 1 -app.overlap_samples;
            app.ptrOffsetSamp = app.ptrOnsetSamp + app.window_samples;
           
            %            end
            if app.ptrOnsetSamp<1
                app.ptrOnsetSamp  = 1;
                app.ptrOffsetSamp = app.ptrOnsetSamp + app.window_samples+app.overlap_samples;
                app.ptrWindow = 1;
                %return;
            end
            if app.ptrOffsetSamp > app.n_samples
                %app.ptrOffsetSamp = app.ptrOnsetSamp-1;
                %app.ptrOnsetSamp = app.ptrOnsetSamp - (app.window_samples+app.overlap_samples) -1;
                app.ptrChannel = app.ptrChannel+1;
                app.ptrWindow = 1;
                app.ptrOnsetSamp  = 1;
                app.ptrOffsetSamp = app.ptrOnsetSamp + app.window_samples+app.overlap_samples;
                %return;
            end
            
            if (app.ptrChannel <= app.nChannels)
                updateAxes(app,app.ptrChannel,[app.ptrOnsetSamp,app.ptrOffsetSamp]);
                cast = app.scoresMat{app.ptrSubject}(app.ptrChannel,app.ptrWindow);
                %fprintf('Window: %i  vote:%i\n',app.ptrWindow,cast);
                switch cast
                    case -1
                        app.LowButton.FontWeight = 'bold';
                        app.HighButton.FontWeight = 'normal';
                        app.LowButton.FontSize = 22;
                        app.HighButton.FontSize = 16;
                    case 1
                        app.HighButton.FontWeight = 'bold';
                        app.LowButton.FontWeight = 'normal';
                        app.LowButton.FontSize = 16;
                        app.HighButton.FontSize = 22;
                    case 0
                        app.LowButton.FontWeight = 'normal';
                        app.HighButton.FontWeight = 'normal';
                        app.LowButton.FontSize = 16;
                        app.HighButton.FontSize = 16;
                end
            else
                idxsm = checkVotes(app);
                if isempty(idxsm)
                    disp('Scan finished, proceed to the next one.');
                else
                    disp('Some windows have not been evaluated yet.');
                end
            end
            
        end
        
        function idxsm = checkVotes(app)
           iSubj = app.ptrSubject;
           iScoresMat = app.scoresMat{iSubj};
           idxsm = find(iScoresMat);
            
        end
        function initApp(app)
            app.ptrChannel = 1;
            app.ptrWindow = 0;
            %app.ptrOnsetSamp = 1;
            app.ptrSubject = 1;
            app.windowSec = 5;
            app.windowOverlap = 0.0;
            app.scoresFile = 'scores.mat';
            app.nChannelsFx = 16;
            app.nSubjFx = 18;
            app.dotNirsPath = uigetdir(app.dotNirsPath,'Select .nirs files folder');
            if app.dotNirsPath==0
                 app.dotNirsPath = pwd;
            else
                cd(app.dotNirsPath);
                %drawnow;
                %figure(app.UIFigure);
                if exist('scores.mat','file')==2
                    disp('Scores.mat found, starting from the last visited window.');
                    load(app.scoresFile);
                    app.scoresMat = scoresMat;
                    app.flagScFl = 1;
                else
                    disp('Scores.mat was not found, and it will be created.');
                    app.scoresMat = cell(app.nSubjFx,1);
                    scoresMat = app.scoresMat;
                    save('scores.mat','scoresMat');
                    disp('Starting from the first window/channel/subject.');
                    app.flagScFl = 0;
                end
            end
        end
    end
    
    % App creation and deletion
    methods (Access = public)
        
        % Construct app
        function app = QualityQuestionnaire_code
            
            initApp(app);
            % Create UIFigure and components
            createComponents(app);
            LoadDotNirsDir(app);
            if app.flagScFl == 1
                %app.ptrOnsetSamp = 1;
                app.ptrOffsetSamp = 1;
                for i=1:length(app.scoresMat)
                   if ~isempty(find(app.scoresMat{i}~=0,1,'first'))
                       [app.ptrWindow,app.ptrChannel]=find(app.scoresMat{i}(:,:)'==0,1,'first');
                       app.ptrSubject = i;
                   end
                end
                curNode=app.treeDotNirs.Children.Children(app.ptrSubject);
                app.treeDotNirs.SelectedNodes = curNode;
                go2scan(app,[]);
            end
            % Register the app with App Designer
            registerApp(app, app.UIFigure);
            
            if nargout == 0
                clear app
            end
        end
        
        % Code that executes before app deletion
        function delete(app)
            scoresMat = app.scoresMat;
            save('scores.mat','scoresMat');
            disp('Scores saved in scores.mat file.');
            % Delete UIFigure when app is deleted
            delete(app.UIFigure);
        end
    end
end
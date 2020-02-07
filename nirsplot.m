function report_table = nirsplot(rawDotNirs,fcut_,window_,overlap_,q_threshold,cond_mask,lambda_mask_)
% raw: raw data in Homer format. Ex: raw = load('nirx_sample.nirs','-mat')
% fcut: 1x2 array [fmin fmax] representing the bandpass of the cardiac pulsation (default [0.5 2.5])
% window: length in seconds of the window to partition the signal with (defaut: 5)
% overlap: fraction overlap (0..0.99) between adjacent windows (default: 0, no overlap)
% lambda_mask: binary array mapping the selected two wavelength to correlate
% (default: [1 1 ...], the first two encountered, no matter how many there are)

close all
%global nSources nDetectors nChannels secDur sampDur raw
%global qMats qltyThld mergewoiFlag SDMeasListAct

report_table = [];

if nargin<5
    q_threshold = 0.75;
end

if nargin < 4
    overlap_ = 0;
end

if nargin < 3
    window_ = 5;
end

if nargin < 2
    fcut_ = [0.5 2.5];
end

if nargin < 1
    nirsplotLoadFileGUI;
    return;
end

% Creating 's' variable (stimuli matrix) from the information in StimDesign
frequency_samp = 1/mean(diff(rawDotNirs.t));
if ~isfield(rawDotNirs,'s')
    if isfield(rawDotNirs,'StimDesign')
        nStim = length(rawDotNirs.StimDesign);
        sTmp = zeros(size(rawDotNirs.d,1),nStim);
        for iStim = 1:nStim
            sTmp(floor(rawDotNirs.StimDesign(iStim).onset * frequency_samp),iStim) = 1;
        end
        rawDotNirs.s = sTmp;
        clear sTmp;
    else
        error('Stimuli information is not available.');
    end
end

lambdas_ = unique(rawDotNirs.SD.MeasList(:,4));
if nargin < 7
    lambda_mask_ = zeros(length(lambdas_),1);
    lambda_mask_(1)=1;
    lambda_mask_(2)=1;
end

if nargin<6
  cond_mask = ones(1,size(rawDotNirs.s,2));
end

% Create GUI
[main_fig_axes,main_fig] = createGUI();
nirsplot_parameters.fcut = fcut_;
nirsplot_parameters.window = window_;
nirsplot_parameters.overlap = overlap_;
nirsplot_parameters.lambda_mask = lambda_mask_;
nirsplot_parameters.lambdas = lambdas_;
nirsplot_parameters.mergewoi_flag = true;
nirsplot_parameters.quality_threshold = q_threshold;
nirsplot_parameters.n_channels = size(rawDotNirs.d,2)/2;
nirsplot_parameters.n_sources = size(rawDotNirs.SD.SrcPos,1);
nirsplot_parameters.n_detectors = size(rawDotNirs.SD.DetPos,1);
nirsplot_parameters.s = rawDotNirs.s;
nirsplot_parameters.t = rawDotNirs.t;
nirsplot_parameters.fs = frequency_samp;
nirsplot_parameters.mergewoiFlag = true;
nirsplot_parameters.cond_mask = cond_mask;
nirsplot_parameters.save_report_table = false;
nirsplot_parameters.sclAlpha = 0.2;
nirsplot_parameters.main_fig_axes = main_fig_axes;
setappdata(main_fig,'nirsplot_parameters',nirsplot_parameters);
setappdata(main_fig,'rawDotNirs',rawDotNirs);

% Computation
%[quality_matrices] = qualityCompute(raw,fcut,window,overlap,lambda_mask,lambdas,mergewoi_flag);
[quality_matrices] = qualityCompute(main_fig);
nirsplot_parameters.quality_matrices = quality_matrices;
setappdata(main_fig,'nirsplot_parameters',nirsplot_parameters);

main_fig.Visible = 'on';

updateQPlots(main_fig);

if nirsplot_parameters.save_report_table == true
    report_table = saveQuality(quality_matrices);
end
% Wait for calls


%% -------------------------------------------------------------------------
    function [main_fig_axes,main_fig] = createGUI()
        
        % Main figure container
        pos.main = [0.125 0.05 0.75 0.8]; % left, bottom, width, height
        main_fig = figure('Units','normalized',...
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
        main_fig_axes.inspector = axes(main_fig,'Units','normalized',...
            'Position',pos.inspAx,...
            'Title','Inspector');
        main_fig_axes.inspector.XLabel.String = 'Time (s)';
        main_fig_axes.inspector.YLabel.String = 'Channel #';
        
        
        pos.comboAx = [myAxDim.xSep,myAxDim.ySep+(1*(myAxDim.height+myAxDim.ySep)),...
            myAxDim.width,myAxDim.height];
        main_fig_axes.combo = axes(main_fig,'Units','normalized',...
            'Position',pos.comboAx,...
            'Title','Overall quality');
        main_fig_axes.combo.YLabel.String = 'Channel #';
        
        pos.powerAx = [myAxDim.xSep,myAxDim.ySep+(2*(myAxDim.height+myAxDim.ySep)),...
            myAxDim.width,myAxDim.height];
        main_fig_axes.power = axes(main_fig,'Units','normalized',...
            'Position',pos.powerAx,...
            'Title','Power peak');
        main_fig_axes.power.YLabel.String = 'Channel #';
        
        pos.sciAx = [myAxDim.xSep,myAxDim.ySep+(3*(myAxDim.height+myAxDim.ySep)),...
            myAxDim.width,myAxDim.height];
        main_fig_axes.sci = axes(main_fig,'Units','normalized',...
            'Position',pos.sciAx,...
            'Title','SCI');
        main_fig_axes.sci.YLabel.String = 'Channel #';
        
        pos.inspectBtn = [myAxDim.xSep, (myAxDim.height+myAxDim.ySep),...
            0.08, myAxDim.ySep*0.8];
        inspectBtn = uicontrol(main_fig,'Style', 'pushbutton', 'String', 'Inspect',...
            'FontSize',14,'FontWeight','bold','Units','normalized','Position', pos.inspectBtn,...
            'Callback', @inspectActive);
        
        pos.helpBtn = [(pos.inspectBtn(1)+pos.inspectBtn(3))*1.15,...
            pos.inspectBtn(2),0.05,pos.inspectBtn(4)];
        helpBtn = uicontrol(main_fig,'Style','pushbutton','String','?',...
            'FontSize',14,'FontWeight','bold','Units','normalized','Position',...
            pos.helpBtn,'Callback', @showHelp);
        
        pos.chSelBtn = [(pos.helpBtn(1)+pos.helpBtn(3))*1.15,...
            pos.inspectBtn(2),0.15,pos.inspectBtn(4)];
        chSelBtn = uicontrol(main_fig,'Style','pushbutton','String','Channel selection',...
            'FontSize',14,'FontWeight','bold','Units','normalized','Position',...
            pos.chSelBtn,'Callback', @selectGoodChannels);
        
        % pos.woiSelBtn = [(pos.chSelBtn(1)+pos.chSelBtn(3))*1.15,...
        %     pos.inspectBtn(2),0.15,pos.inspectBtn(4)];
        % woiSelBtn = uicontrol(mainFig,'Style','pushbutton','String','WOI selection',...
        %     'FontSize',14,'FontWeight','bold','Units','normalized','Position',...
        %     pos.woiSelBtn,'Callback', @selectGoodWOI);
        
        pos.SaveBtn = [2.5*(pos.chSelBtn(3)+pos.chSelBtn(3)),...
            pos.inspectBtn(2),0.1,pos.inspectBtn(4)];
        SaveBtn = uicontrol(main_fig,'Style','pushbutton','String','Save .nirs',...
            'FontSize',14,'FontWeight','bold','Units','normalized','Position',...
            pos.SaveBtn,'Callback', @save2dotnirs);
        
        pos.AdvView = [myAxDim.width,...
            pos.inspectBtn(2)-0.05,...
            0.1,...
            pos.inspectBtn(4)];        
        viewMCheckB = uicontrol(main_fig,'Style','checkbox','String','Advanced view',...
            'FontSize',10,'FontWeight','bold','Units','normalized','Position',...
            pos.AdvView,'Callback', @viewMode,'Tag','viewMCheckB');
        
        debFig=figure(3);
        main_fig_axes.debAx =  axes(debFig);
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
        s = nirsplot_param.s;
        t = nirsplot_param.t;
        
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
                    xLimWindow = [1,(qMats.sampPerWindow*qMats.n_windows)];
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
                if iChannel>0 && iChannel<=n_channels && iWindow>0 && iWindow<=qMats.n_windows
                    updateIPlot(source,iChannel,xLimWindow,iWindow,s,t);
                    pastWindow = iWindow;
                    pastChannel = iChannel;
                end
            end
        end
    end

%% -------------------------------------------------------------------------
    function updateQPlots(main_fig)
        % UpdatePlot updates the quality plots with the 'qualityMats' input arg
        % bad channels are ploted according to 'plot_bad' flag
        nirsplot_param = getappdata(main_fig,'nirsplot_parameters');
        qMats = nirsplot_param.quality_matrices;
        myAxes = nirsplot_param.main_fig_axes;
        n_channels = nirsplot_param.n_channels;
        woi = nirsplot_param.quality_matrices.woi;
        viewMCheckB = findobj('Tag','viewMCheckB');
        
        %Unpacking
        sci_array = qMats.sci_array;
        power_array = qMats.power_array;
        combo_array = qMats.combo_array;
        combo_array_expanded = qMats.combo_array_expanded;
        
        sclAlpha = 0.2;
        
        viewMode(viewMCheckB,[]);
%         % Scalp Contact Index
%         imSci = imagesc(myAxes.sci,sci_array);
%         myAxes.sci.CLim = [0,1];
%         myAxes.sci.YLim =[1, n_channels];
%         %myAxes.sci.XLim =[1, size(sci_array,2)];
%         % myAxes.sci.XTick = round(linspace(1,qMats.n_windows,8));
%         % myAxes.sci.XTickLabel = num2mstr(round(linspace(1,qMats.n_windows*qMats.sampPerWindow,8)));
%         colormap(myAxes.sci,"gray");
%         colorbar(myAxes.sci,"eastoutside","Ticks",[0 0.8 1]);
%         myAxes.sci.YLabel.String = 'Channel #';
%         myAxes.sci.YLabel.FontWeight = 'bold';
        
%         % Power peak
%         imPower = imagesc(myAxes.power,power_array);
%         myAxes.power.CLim = [0, 0.12];
%         myAxes.power.YLim =[1, n_channels];
%         %myAxes.power.XLim =[1, size(power_array,2)];
%         colormap(myAxes.power,"gray");
%         colorbar(myAxes.power,"eastoutside","Ticks",[0 0.1 0.12]);
%         %        h.TickLabels ={'0','0.10*','0.12'};
%         myAxes.power.YLabel.String = 'Channel #';
%         myAxes.power.YLabel.FontWeight = 'bold';
        
%         % Combo quality
%         imCombo = imagesc(myAxes.combo,combo_array_expanded);
%         myAxes.combo.CLim = [0, 3];
%         myAxes.combo.YLim =[1, n_channels];
%         %myAxes.combo.XLim =[1, size(combo_array,2)];
%         %colormap(myAxes.combo,[0 0 0;1 1 1]);
%         % SCI,Power   combo_array_expanded    QualityColor
%         %  0,0              0                   [0 0    0]
%         %  0,1              1                   [1 0    0]
%         %  2,0              2                   [1 0.64 0]
%         %  2,1              3                   [0 1    0]
%         colormap(myAxes.combo,[0 0 0; 1 0 0; 1 0.64 0; 0 1 0]);
%         colorbar(myAxes.combo,"eastoutside","Ticks",[0 1 2 3],...
%             'TickLabels',...
%             {[char(hex2dec('2717')),'SCI  ', char(hex2dec('2717')),'Power'],...
%              [char(hex2dec('2717')),'SCI  ', char(hex2dec('2713')),'Power'],...
%              [char(hex2dec('2713')),'SCI  ', char(hex2dec('2717')),'Power'],...
%              [char(hex2dec('2713')),'SCI  ', char(hex2dec('2713')),'Power']});
%         myAxes.combo.YLabel.String = 'Channel #';
%         myAxes.combo.YLabel.FontWeight = 'bold';
        
        % For visual consistency among axes
        myAxes.inspector.YLimMode = 'manual';
        myAxes.inspector.YLabel.String = 'Channel #';
        myAxes.inspector.XLabel.String = 'Time (s)';
        myAxes.inspector.YLabel.FontWeight = 'bold';
        myAxes.inspector.XLabel.FontWeight = 'bold';
        colorbar(myAxes.inspector,'Visible','off');
        
        woiMatrgb = zeros(n_channels,qMats.n_windows,3);
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

%% -------------------------------------------------------------------------
    function updateIPlot(source,iChannel,xLimWindow,iWindow,s,t)
        
        nirsplot_param = getappdata(source.Parent,'nirsplot_parameters');
        qMats = nirsplot_param.quality_matrices;
        myAxes = nirsplot_param.main_fig_axes;
        n_channels = nirsplot_param.n_channels;
        conditions_mask = nirsplot_param.cond_mask;
        woi = nirsplot_param.quality_matrices.woi;
        fs = nirsplot_param.fs;
        
        sclAlpha = 0.2;
        myAxes.inspector.XLim= [t(xLimWindow(1)),t(xLimWindow(2))];
%        YLimStd = [min(qMats.cardiac_data(:,xLimWindow(1):xLimWindow(2),iChannel),[],'all'),...
%            max(qMats.cardiac_data(:,xLimWindow(1):xLimWindow(2),iChannel),[],'all')]*1.05;
%        YLimStd = [min(qMats.cardiac_data(:,:,iChannel),[],'all'),...
%            max(qMats.cardiac_data(:,:,iChannel),[],'all')]*1.05;        
        YLimStd = [min(qMats.cardiac_data(:,xLimWindow(1):xLimWindow(2),iChannel),[],'all'),...
            max(qMats.cardiac_data(:,xLimWindow(1):xLimWindow(2),iChannel),[],'all')]*1.05;
        XLimStd = myAxes.inspector.XLim;
        
        % qMats.cardiac_data;  % #Lambda x time x #channels
        wl1Norm = (qMats.cardiac_data(1,xLimWindow(1):xLimWindow(2),iChannel)-YLimStd(1))./ (YLimStd(2)-YLimStd(1));
        wl2Norm = (qMats.cardiac_data(2,xLimWindow(1):xLimWindow(2),iChannel)-YLimStd(1))./ (YLimStd(2)-YLimStd(1));
        YLimStd = [0,1.05];
        cla(myAxes.inspector);
%        plot(myAxes.inspector,t(xLimWindow(1):xLimWindow(2)),...
%            qMats.cardiac_data(1,xLimWindow(1):xLimWindow(2),iChannel),,'-b');
        plot(myAxes.inspector,t(xLimWindow(1):xLimWindow(2)),...            
            wl1Norm,'-b');
        hold(myAxes.inspector,'on');
%        plot(myAxes.inspector,t(xLimWindow(1):xLimWindow(2)),...
%            qMats.cardiac_data(2,xLimWindow(1):xLimWindow(2),iChannel),'-r');
        plot(myAxes.inspector,t(xLimWindow(1):xLimWindow(2)),...            
            wl2Norm,'-r');
        
        %yLab = myAxes.inspector.YLim(1)+((myAxes.inspector.YLim(2)-myAxes.inspector.YLim(1))*0.71);
        %xLab = myAxes.inspector.XLim(1)+((myAxes.inspector.XLim(2)-myAxes.inspector.XLim(1))*0.93);
        %text(myAxes.inspector,xLab,yLab,['Channel:',num2str(iChannel)]);
        strLgnds = {'Lambda 1','Lambda 2'};
        updateQPlots(source.Parent);

        if (xLimWindow(2)-xLimWindow(1)+1) == (qMats.n_windows*qMats.sampPerWindow)
            xRect = 0.5; %Because of the offset at the begining of a window
            yRect = iChannel-0.5;
            wRect = qMats.n_windows;
            hRect = 1;
            poiMatrgb = zeros(n_channels,xLimWindow(2),3);
            poiMatrgb(:,:,2) = repmat(repelem(woi.mat(1,:),qMats.sampPerWindow),n_channels,1);
            alphaMat = poiMatrgb(:,:,2) * sclAlpha;
            
            impoiMat = imagesc(myAxes.inspector,'XData',...
                [t(xLimWindow(1)),t(xLimWindow(2))],...
                'YData',YLimStd,'CData',poiMatrgb,'AlphaData',alphaMat);
            
            %onsets
            c = sum(conditions_mask);
            COI = find(conditions_mask);
            colorOnsets = [linspace(0.5,1,c)',...
                linspace(0,0.5,c)',linspace(1,0,c)'];
            for j=1:c
                %mapping from 0,1 to 0,25%ofPeakToPeak
                yOnset = (s(xLimWindow(1):xLimWindow(2),COI(j))*(YLimStd(2)-YLimStd(1))*0.25)-abs(YLimStd(1));
                plot(myAxes.inspector,t(xLimWindow(1):xLimWindow(2)),...
                    yOnset,'LineWidth',2,...
                    'Color',colorOnsets(j,:));
                strLgnds(2+j) = {['Cond ',num2str(COI(j))]};
            end
            
        else
            xRect = iWindow-0.5;
            yRect = iChannel-0.5;
            wRect = 1;
            hRect = 1;
            
            fprintf('SCI:%.2f \t Power:%.2f\n',qMats.sci_array(iChannel,iWindow),qMats.power_array(iChannel,iWindow));
            textHAlign = 'left';
            textVAlign = 'top';
            if xRect > (qMats.n_windows/2)
               textHAlign = 'right'; 
            end
            if yRect > (n_channels/2)
               textVAlign = 'bottom';
            end
            text(myAxes.power,xRect,yRect,num2str(qMats.power_array(iChannel,iWindow),2),...
                'Color','red','FontSize',10,'FontWeight','bold','BackgroundColor','#FFFF00',...
                'Margin',1,'Clipping','on',...
                'HorizontalAlignment',textHAlign,'VerticalAlignment',textVAlign);
            text(myAxes.sci,xRect,yRect,num2str(qMats.sci_array(iChannel,iWindow),2),...
                'Color','red','FontSize',10,'FontWeight','bold','BackgroundColor','#FFFF00',...
                'Margin',1,'Clipping','on',...
                'HorizontalAlignment',textHAlign,'VerticalAlignment',textVAlign);
            %graphical debug
            graphicDebug(qMats.cardiac_data(1,xLimWindow(1):xLimWindow(2),iChannel),...
                qMats.cardiac_data(2,xLimWindow(1):xLimWindow(2),iChannel),fs);
        end
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

%% -------------------------------------------------------------------------
    function [qualityMats] = qualityCompute(main_fig)
        raw = getappdata(main_fig,'rawDotNirs');
        nirsplot_param = getappdata(main_fig,'nirsplot_parameters');
        fcut = nirsplot_param.fcut;
        window = nirsplot_param.window;
        overlap = nirsplot_param.overlap;
        lambda_mask = nirsplot_param.lambda_mask;
        lambdas = nirsplot_param.lambdas;
        n_channels = nirsplot_param.n_channels;
        qltyThld = nirsplot_param.quality_threshold;
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
        
        nirs_data = zeros(length(lambdas),size(raw.d,1),n_channels);
        cardiac_data = zeros(length(lambdas),size(raw.d,1),n_channels); % Lambdas x time x channels
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
        cardiac_windows = zeros(length(lambdas),window_samples,n_channels,n_windows);
        for j = 1:n_windows
            interval = (j-1)*window_samples-(j-1)*(overlap_samples)+1 : j*window_samples-(j-1)*(overlap_samples);
            cardiac_windows(:,:,:,j) = cardiac_data(:,interval,:);
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
        combo_array_expanded = 2*(sci_array >= 0.8) + (power_array >= 0.10);
        
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
        
%         f2=figure(2);
%         af2=axes(f2);
%         idxActiveWindows = idxPoi;
%         idxUnactiveWindows = ~idxPoi;
%         bp=bar(af2,good_combo_window);
%         bp.FaceColor = 'flat';
%         ylim([0,1]);
%         bp.CData(idxUnactiveWindows,:) = repmat([0.6, 0.6, 0.6],sum(idxUnactiveWindows),1);
%         bp.CData(idxActiveWindows,:) = repmat([0 1 0],sum(idxActiveWindows),1);
      
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
        qualityMats.cardiac_data = cardiac_data;
        qualityMats.good_combo_link = good_combo_link;
        qualityMats.good_combo_window = good_combo_window;
        qualityMats.woi = woi;
        %
    end




%% -------------------------------------------------------------------------
    function [woi] = getWOI(window_samples,n_windows,nirsplot_parameters)
        % Assuming no overlaped trials
        % The maximum number of allowed samples is window_samples*n_windows to consider
        % an integer number of windows, module(total_samples,n_windows) = 0
        
        fs = nirsplot_parameters.fs;
        n_channels = nirsplot_parameters.n_channels;
        s = nirsplot_parameters.s;
        t = nirsplot_parameters.t;
        mergewoi_flag = nirsplot_parameters.mergewoi_flag;
        n_channels = nirsplot_parameters.n_channels;
        conditions_mask = logical(nirsplot_parameters.cond_mask);        

        allowed_samp = window_samples*n_windows;
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
        blckDurWind = floor(blckDurSamp/window_samples);
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
            startWOI = floor(startPOI/window_samples);
            if startWOI==0
                startWOI = 1;
            end
            
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
        % poiMat_ = repmat(poi,n_channels,1);
        
        woiblank = 0;
        idxInit = [];
        %woitmp = woi_array;
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
    end

%% -------------------------------------------------------------------------
    function dotNirsOutput = selectGoodChannels(source, events)
        bpGoodQuality(source.Parent);
        uiwait(source.Parent);
        nirsplot_param = getappdata(source.Parent,'nirsplot_parameters');
        disp(['Threshold was changed to ',num2str(nirsplot_param.quality_threshold)]);
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
         raw = getappdata(source.Parent,'rawDotNirs');
        active_channels = nirsplot_param.quality_matrices.active_channels;
        % saving the indices of good-quality channels     
            SD = raw.SD;
            SD.MeasListAct = [active_channels; ones(size(SD.MeasList,1)/2,1)];
            t = raw.t;
            d = raw.d;            
            s = raw.s;            
            aux = raw.aux;
            tIncMan = ones(length(t),1);
            save('dotNirs_nirsplot.nirs','SD','t','d','s','aux','tIncMan');
            %Notify to the user if the new file was succesfully created     
            saving_status = exist('dotNirs_nirsplot.nirs','file');
            if saving_status
                msgbox('Operation Completed','Success');
            else
                msgbox('Operation Failed','Error');
            end

    end


    function viewMode(source,event)
        nirsplot_param = getappdata(source.Parent,'nirsplot_parameters');
        qMats = nirsplot_param.quality_matrices;
        myAxes = nirsplot_param.main_fig_axes;
        n_channels = nirsplot_param.n_channels;
        
        sci_array = qMats.sci_array;
        power_array = qMats.power_array;
        combo_array = qMats.combo_array;
        combo_array_expanded = qMats.combo_array_expanded;
      
        checked = source.Value;
        
        if checked
            % Scalp Contact Index
            imagesc(myAxes.sci,sci_array);
            myAxes.sci.CLim = [0,1];
            myAxes.sci.YLim =[1, n_channels];
            %myAxes.sci.XLim =[1, size(sci_array,2)];
            % myAxes.sci.XTick = round(linspace(1,qMats.n_windows,8));
            % myAxes.sci.XTickLabel = num2mstr(round(linspace(1,qMats.n_windows*qMats.sampPerWindow,8)));
            colormap(myAxes.sci,"hot");
            colorbar(myAxes.sci,"eastoutside","Ticks",[0 0.8 1]);
            myAxes.sci.YLabel.String = 'Channel #';
            myAxes.sci.YLabel.FontWeight = 'bold';

            % Power peak
            imPower = imagesc(myAxes.power,power_array);
            myAxes.power.CLim = [0, 0.12];
            myAxes.power.YLim =[1, n_channels];
            %myAxes.power.XLim =[1, size(power_array,2)];
            colormap(myAxes.power,"hot");
            colorbar(myAxes.power,"eastoutside","Ticks",[0 0.1 0.12]);
            %        h.TickLabels ={'0','0.10*','0.12'};
            myAxes.power.YLabel.String = 'Channel #';
            myAxes.power.YLabel.FontWeight = 'bold';
            
            
            % Combo quality
            imagesc(myAxes.combo,combo_array_expanded);
            myAxes.combo.CLim = [0, 3];
            myAxes.combo.YLim =[1, n_channels];
            %myAxes.combo.XLim =[1, size(combo_array,2)];
            %colormap(myAxes.combo,[0 0 0;1 1 1]);
            % SCI,Power   combo_array_expanded    QualityColor
            %  0,0              0                   [0 0    0]
            %  0,1              1                   [1 0    0]
            %  2,0              2                   [1 0.64 0]
            %  2,1              3                   [0 1    0]
            colormap(myAxes.combo,[0 0 0; 1 0 0; 1 0.64 0; 0 1 0]);
            colorbar(myAxes.combo,"eastoutside","Ticks",[0 1 2 3],...
                'TickLabels',...
                {[char(hex2dec('2717')),'SCI  ', char(hex2dec('2717')),'Power'],...
                 [char(hex2dec('2717')),'SCI  ', char(hex2dec('2713')),'Power'],...
                 [char(hex2dec('2713')),'SCI  ', char(hex2dec('2717')),'Power'],...
                 [char(hex2dec('2713')),'SCI  ', char(hex2dec('2713')),'Power']});
            myAxes.combo.YLabel.String = 'Channel #';
            myAxes.combo.YLabel.FontWeight = 'bold';
        else
            % Scalp Contact Index
            imagesc(myAxes.sci,sci_array);
            myAxes.sci.CLim = [0,1];
            myAxes.sci.YLim =[1, n_channels];
            %myAxes.sci.XLim =[1, size(sci_array,2)];
            % myAxes.sci.XTick = round(linspace(1,qMats.n_windows,8));
            % myAxes.sci.XTickLabel = num2mstr(round(linspace(1,qMats.n_windows*qMats.sampPerWindow,8)));
            colormap(myAxes.sci,"gray");
            colorbar(myAxes.sci,"eastoutside","Ticks",[0 0.8 1]);
            myAxes.sci.YLabel.String = 'Channel #';
            myAxes.sci.YLabel.FontWeight = 'bold';

            % Power peak
            imagesc(myAxes.power,power_array);
            myAxes.power.CLim = [0, 0.12];
            myAxes.power.YLim =[1, n_channels];
            %myAxes.power.XLim =[1, size(power_array,2)];
            colormap(myAxes.power,"gray");
            colorbar(myAxes.power,"eastoutside","Ticks",[0 0.1 0.12]);
            %        h.TickLabels ={'0','0.10*','0.12'};
            myAxes.power.YLabel.String = 'Channel #';
            myAxes.power.YLabel.FontWeight = 'bold';

            % Combo quality
            imagesc(myAxes.combo,combo_array);
            myAxes.combo.CLim = [0, 1];
            myAxes.combo.YLim =[1, n_channels];
            %myAxes.combo.XLim =[1, size(combo_array,2)];
            colormap(myAxes.combo,[0 0 0;1 1 1]);
            colorbar(myAxes.combo,"eastoutside","Ticks",[0 1]);
            myAxes.combo.YLabel.String = 'Channel #';
            myAxes.combo.YLabel.FontWeight = 'bold';
        end
    end

    function graphicDebug(window1,window2,fs)
        %cross-correlate the two wavelength signals - both should have cardiac pulsations
        similarity = xcorr(window1,window2,'unbiased');
        if any(abs(similarity)>eps)
            % this makes the SCI=1 at lag zero when x1=x2 AND makes the power estimate independent of signal length, amplitude and Fs
            similarity = length(window1)*similarity./sqrt(sum(abs(window1).^2)*sum(abs(window2).^2));
        else
            warning('Similarity results close to zero');
        end
        [pxx,f] = periodogram(similarity,hamming(length(similarity)),length(similarity),fs,'power');
        f3=figure(3);
        clf(f3);
        subplot(2,1,1);
        plot(f3.Children(1),similarity);
        subplot(2,1,2);
        plot(f3.Children(1),pxx);
    end
    
end %end of nirsplot function definition
function generateQReport(reportTable)
%generateQReport generates a graphical quality report from QTNIRS results
% Arguments
%
% Input
% reportTable   : matrix with Hb data [#samples x (#channels)] from WL1
%
% Output
% none
%
% Notes:
% 
% Version 0.1   Initial state
%
%

nscans = numel(reportTable);
nchannels = size(reportTable(1).MeasList,1)/2;

% Group-level quality report
if nscans>1
    report_mat=zeros(1,nchannels);
    allSubj = zeros(nscans,nchannels);
    scanInfonames = cell(1,nscans);
    master_threshold = reportTable(1).thresholds.quality;
    for i =1:nscans
        report_mat = report_mat + (reportTable(i).good_combo_link(:,3)>=master_threshold)';
        allSubj(i,:) = reportTable(i).good_combo_link(:,3);
        scanInfonames{i} = reportTable(i).scanInfo; 
    end
    %Bar plot: display the number of scans that achieved a quality above
    %the threshold for every channels.
    f1 = figure('Name','Data quality report: channels above quality threshold','NumberTitle','off');
    bar(report_mat);
    title(sprintf('%s (SCI=%.2f,PSP=%.2f)','Group-level high-quality channels',...
        reportTable(1).thresholds.sci,reportTable(1).thresholds.peakpower),'Interpreter','none');
    xticks(1:size(reportTable(1).MeasList,1)/2);
    xticklabels(string([num2str(reportTable(1).good_combo_link(:,1)),...
        repmat('-',size(reportTable(1).good_combo_link,1),1),...
        num2str(reportTable(1).good_combo_link(:,2),2)]));
    xtickangle(35);
    xlabel('Source-Detector pair');
    ylabel('#Scans');
    ylim([0, nscans]);
    legend(['HQ above ',num2str(master_threshold*100),'% of time']);
    f1.WindowState = 'maximized';
    saveas(f1,'QC_Group_good-channels_scanwise.png');
    close(f1);

    %Heat map: display the achieved quality for every channel and every
    %subject
    f2 = figure('Name','Data quality report: channel quality','NumberTitle','off');
    imagesc(allSubj.*100);
    colormap(bone(20));
    title('Group-level quality evaluation');
    xticks(1:nchannels);
    xticklabels(string([num2str(reportTable(1).good_combo_link(:,1)),...
        repmat('-',size(reportTable(1).good_combo_link,1),1),...
        num2str(reportTable(1).good_combo_link(:,2),2)]));
    xtickangle(35);
    xlabel('Source-Detector pair');
    yticks(1:nscans);
    yticklabels(scanInfonames);
    ylabel('Scans');
    c = colorbar();
    c.Label.String = 'Recording time with good quality (%)';
    a = gca;
    a.TickLabelInterpreter = 'none';
    f2.WindowState = 'maximized';
    saveas(f2,'QC_Group_quality_channelwise.png');
    close(f2);

    % QT style
    for i =1:nscans
        approx_scan_duration = reportTable(i).sampPerWindow*reportTable(i).n_windows / reportTable(i).fs;
        f3 = figure('Name','Data quality report: SCI, PSP, (SCI and PSP)','NumberTitle','off');
        
        if approx_scan_duration > 2000
            shift_in_sec = 200;
            interval = shift_in_sec*reportTable(i).fs/reportTable(i).sampPerWindow;
        elseif approx_scan_duration > 1000
            shift_in_sec = 100;
            interval = shift_in_sec*reportTable(i).fs/reportTable(i).sampPerWindow;
        else
            shift_in_sec = 50;
            interval = shift_in_sec*reportTable(i).fs/reportTable(i).sampPerWindow;
        end
        

        window_seconds = reportTable(i).sampPerWindow/reportTable(i).fs;
        ticksVals = interval:interval:reportTable(i).n_windows;
        ticksLab = round(ticksVals*window_seconds);      

        subplot(3,1,1);
        imagesc(reportTable(i).sci_array>=reportTable(i).thresholds.sci);
        colormap([0 0 0; 1 1 1]);
        asci = gca;
        asci.CLim = [0,1];
        colorbar('eastoutside',...
            'Tag','colorb_sci',...
            'Ticks',[0.25 0.75],...
            'Limits',[0,1],...
            'TickLabels',{'Bad','Good'});
        asci.XAxis.TickValues=ticksVals;
        asci.XAxis.TickLabels=split(num2str(ticksLab));
        asci.YLabel.String = 'Channel #';
        asci.YLabel.FontWeight = 'bold';
        title(sprintf('%s - SCI >= %.2f',reportTable(i).scanInfo,reportTable(i).thresholds.sci),"Interpreter","none");

        % Power peak
        subplot(3,1,2);
        imagesc(reportTable(i).power_array>=reportTable(i).thresholds.peakpower);
        colormap([0 0 0; 1 1 1]);
        apsp = gca;
        apsp.CLim = [0,1];
        colorbar('eastoutside',...
            'Tag','colorb_psp',...
            'Ticks',[0.25 0.75],...
            'Limits',[0,1],...
            'TickLabels',{'Bad','Good'});
        apsp.XAxis.TickValues=ticksVals;
        apsp.XAxis.TickLabels=split(num2str(ticksLab));
        apsp.YLabel.String = 'Channel #';
        apsp.YLabel.FontWeight = 'bold';
        title(sprintf('PSP >= %.2f',reportTable(i).thresholds.peakpower));

        % Combo panel
        subplot(3,1,3);
        imagesc(reportTable(i).combo_array_expanded);
        acombo = gca;
        acombo.CLim = [0, 4];        
        qualityColor = [0 0 0; 0 0 0; 1 0 0; 1 1 1;0 0 1];
        colormap(acombo,qualityColor);
        colorbar(acombo,"eastoutside","Ticks",[0.7 1.0 2.05 2.75 3.5],...
                'TickLabels',...
                {[char(hex2dec('2717')),'SCI  ', char(hex2dec('2717')),'Power'],...
                [char(hex2dec('2717')),'SCI  ', char(hex2dec('2713')),'Power'],...
                [char(hex2dec('2713')),'SCI  ', char(hex2dec('2717')),'Power'],...
                [char(hex2dec('2713')),'SCI  ', char(hex2dec('2713')),'Power'],...
                'Saturation'});
        acombo.XAxis.TickValues=ticksVals;
        acombo.XAxis.TickLabels=split(num2str(ticksLab));
        acombo.YLabel.String = 'Channel #';
        acombo.YLabel.FontWeight = 'bold';
        acombo.XLabel.String = 'Time(s)';
        title(sprintf('SCI >= %.2f and PSP >= %.2f',...
            reportTable(i).thresholds.sci,...
            reportTable(i).thresholds.peakpower));
        
        f3.WindowState = 'maximized';
        saveas(f3,sprintf('QC_%s_SCI-PSP-Combo.png',reportTable(i).scanInfo));
        close(f3);
    end
    
else
    report_mat = reportTable.good_combo_link(:,3);

    f1=figure('Name','Channel-level Report','NumberTitle','off');
    bar(report_mat);
    title(sprintf('%s (SCI=%.2f,PSP=%.2f)',reportTable.scanInfo,...
        reportTable.thresholds.sci,reportTable.thresholds.peakpower),'Interpreter','none');
    xticks(1:nchannels);
    xticklabels(string([num2str(reportTable(1).good_combo_link(:,1)),...
        repmat('-',size(reportTable(1).good_combo_link,1),1),...
        num2str(reportTable(1).good_combo_link(:,2),2)]));
    xtickangle(35);
    xlabel('Source-Detector pair');
    ylabel('%Quality');
    
    % SCI and PSP product
    f2=figure('Name','Channel-level Report','NumberTitle','off');
    imagesc((reportTable.sci_array.^2).*reportTable.power_array);
    title(reportTable.scanInfo,'Interpreter','none');
    yticks(1:nchannels);
    yticklabels(string([num2str(reportTable(1).good_combo_link(:,1)),...
        repmat('-',size(reportTable(1).good_combo_link,1),1),...
        num2str(reportTable(1).good_combo_link(:,2),2)]));
    colormap(bone(20));
    c = colorbar();
    c.Limits = [0 0.5];
    c.Label.String = 'SCI^2 \times PSP';
    xlabel('Time (n-sec windows)');
    ylabel('Source-Detector pair');
end

end
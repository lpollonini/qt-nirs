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

else
    report_mat = reportTable.good_combo_link(:,3);
    
    %report_mat = report_mat./max(report_mat);
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
    imagesc(reportTable.sci_array.*reportTable.power_array);
    title(reportTable.scanInfo,'Interpreter','none');
    yticks(1:nchannels);
    yticklabels(string([num2str(reportTable(1).good_combo_link(:,1)),...
        repmat('-',size(reportTable(1).good_combo_link,1),1),...
        num2str(reportTable(1).good_combo_link(:,2),2)]));
    colormap(bone(20));
    c = colorbar();
    c.Limits = [0 0.5];
    c.Label.String = 'SCI \times PSP';
    xlabel('Time (n-sec windows)');
    ylabel('Source-Detector pair');
end

end
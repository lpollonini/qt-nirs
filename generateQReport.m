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

% Group-level quality report
if numel(reportTable)>1
    report_mat=zeros(1,size(reportTable(1).MeasList,1)/2);
    master_threshold = reportTable(1).thresholds.quality;
    for i =1:numel(reportTable)
        report_mat = report_mat + (reportTable(i).good_combo_link(:,3)>=master_threshold)';
    end
    %report_mat = report_mat./max(report_mat);
    f=figure('Name','Channel-level Report','NumberTitle','off');
    bar(report_mat);
    title('High-quality Channels');
    xticks(1:size(reportTable(1).MeasList,1)/2);
    xticklabels(string([num2str(reportTable(1).good_combo_link(:,1)),...
        repmat('-',size(reportTable(1).good_combo_link,1),1),...
        num2str(reportTable(1).good_combo_link(:,2),2)]));
    xtickangle(35);
    xlabel('Channel #');
    ylabel('#Scans');
    legend(['HQ above ',num2str(master_threshold*100),'% of time']);

else
    report_mat=zeros(1,size(reportTable(1).MeasList,1)/2);
    master_threshold = reportTable(1).thresholds.quality;
    
    report_mat = reportTable.good_combo_link(:,3);
    
    %report_mat = report_mat./max(report_mat);
    f=figure('Name','Channel-level Report','NumberTitle','off');
    bar(report_mat);
    title(reportTable.scanInfo);
    xticks(1:size(reportTable(1).MeasList,1)/2);
    xticklabels(string([num2str(reportTable(1).good_combo_link(:,1)),...
        repmat('-',size(reportTable(1).good_combo_link,1),1),...
        num2str(reportTable(1).good_combo_link(:,2),2)]));
    xtickangle(35);
    xlabel('Channel #');
    ylabel('%Quality');
    %legend(['HQ above ',num2str(master_threshold*100),'% of time']);
end

end
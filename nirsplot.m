function report_table = nirsplot(raw,fcut,window,overlap,lambda_mask)

% raw: raw data in Homer format. Ex: raw = load('nirx_sample.nirs','-mat') 
% fcut: 1x2 array [fmin fmax] representing the bandpass of the cardiac pulsation (default [0.5 2.5])
% window: length in seconds of the window to partition the signal with (defaut: 5)
% overlap: fraction overlap (0..0.99) between adjacent windows (default: 0, no overlap)
% lambda_mask: binary array mapping the selected two wavelenth to correlate
% (default: [1 1 ...], the first two encountered, no matter how many there are)

close all

plot_quality_plot = 1;
plot_bad = 0;
save_report_table = 0;

if nargin < 4
    overlap = 0;
end

if nargin < 3
    window = 5;
end

if nargin < 2
    fcut = [0.5 2.5];
end

% if license('test','Distrib_Computing_Toolbox')
%     parpool
% end


%% Set the bandpass filter parameters
fs = 1/mean(diff(raw.t));
fcut_min = fcut(1);
fcut_max = fcut(2);
if fcut_max >= (fs)/2
    fcut_max = (fs)/2 - eps;
    warning(['The highpass cutoff has been reduced from ' num2str(fcut(2)) ' Hz to ' num2str(fcut_max) ' Hz to satisfy the Nyquist sampling criterion']);
end
[B1,A1]=butter(1,[fcut_min*(2/fs) fcut_max*(2/fs)]);

report_table = table();

for i = 1:length(raw)   % For each scan/file
    lambdas = unique(raw(i).SD.MeasList(:,4));
    if nargin < 5
        lambda_mask = zeros(length(lambdas),1);
        lambda_mask(1)=1;
        lambda_mask(2)=1;
    end
    nirs_data = zeros(length(lambdas),size(raw(i).d,1),size(raw(i).d,2)/length(lambdas));
    cardiac_data = zeros(length(lambdas),size(raw(i).d,1),size(raw(i).d,2)/length(lambdas)); % Lambdas x time x channels
    for j = 1:length(lambdas)
        % Filter everything but the cardiac component
        idx = find(raw(i).SD.MeasList(:,4) == lambdas(j));
        nirs_data(j,:,:) = raw(i).d(:,idx);
        filtered_nirs_data=filtfilt(B1,A1,squeeze(nirs_data(j,:,:)));
        cardiac_data(j,:,:)=filtered_nirs_data./repmat(std(filtered_nirs_data,0,1),size(filtered_nirs_data,1),1); % Normalized heartbeat
    end
    overlap_samples = floor(window*fs*overlap);
    window_samples = floor(window*fs);
    n_windows = floor((size(cardiac_data,2)-overlap_samples)/(window_samples-overlap_samples));
    cardiac_data = cardiac_data(find(lambda_mask),:,:);
    sci_array = zeros(size(cardiac_data,3),n_windows);    % Number of optode is from the user's layout, not the machine
    power_array = zeros(size(cardiac_data,3),n_windows);
    fpower_array = zeros(size(cardiac_data,3),n_windows);
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
            similarity = xcorr(squeeze(cardiac_window(1,:,k)),squeeze(cardiac_window(2,:,k)),'unbiased');  %cross-correlate the two wavelength signals - both should have cardiac pulsations
            similarity = length(squeeze(cardiac_window(1,:,k)))*similarity./sqrt(sum(abs(squeeze(cardiac_window(1,:,k))).^2)*sum(abs(squeeze(cardiac_window(2,:,k))).^2));  % this makes the SCI=1 at lag zero when x1=x2 AND makes the power estimate independent of signal length, amplitude and Fs
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
        fpower_array(:,j) = fpower_array_channels;
    end
    
    %% Summary analysis
    mean_sci_link  = mean(sci_array,2);
    std_sci_link  = std(sci_array,0,2);
    good_sci_link = sum(sci_array>0.8,2)/size(sci_array,2);
    mean_sci_window  = mean(sci_array,1);
    std_sci_window  = std(sci_array,0,1);
    good_sci_window = sum(sci_array>0.8,1)/size(sci_array,1);
    
    mean_power_link  = mean(power_array,2);
    std_power_link  = std(power_array,0,2);
    good_power_link = sum(power_array>0.1,2)/size(power_array,2);
    mean_power_window  = mean(power_array,1);
    std_power_window  = std(power_array,0,1);
    good_power_window = sum(power_array>0.1,1)/size(power_array,1);
    
    combo_array = (sci_array >= 0.8) & (power_array >= 0.10);
    mean_combo_link  = mean(combo_array,2);
    std_combo_link  = std(combo_array,0,2);
    good_combo_link = sum(combo_array>0.1,2)/size(combo_array,2);
    mean_combo_window  = mean(combo_array,1);
    std_combo_window  = std(combo_array,0,1);
    good_combo_window = sum(combo_array>0.1,1)/size(combo_array,1);
    
    %% Plot summary
    if plot_quality_plot
        %figure('Name',['File ' num2str(i) ': ' raw(i).description(idx(end)+1:end)],'NumberTitle','off')
        figure
        subplot(4,1,1)
        imagesc(sci_array,[0 1])
%         axis equal
%         axis tight
        h = colorbar;
        h.Ticks=[0 1];
        colormap(gca,gray)
        h.FontSize = 9;
        set(get(h,'label'),'string','SCI');
        ylabel('Channel #')
        set(gca,'XTick',[])
        
        subplot(4,1,2)
        %imagesc(power_array,[0 max(max(power_array))])
        imagesc(power_array,[0 0.3])
%         axis equal
%         axis tight
        h = colorbar;
        %tick_interval = max(max(power_array))/4; 
        %h.Ticks = [0:round(tick_interval,2):round(max(max(power_array)),2)];
        h.Ticks = [0:0.1:0.3];
        colormap(gca,gray)
        h.FontSize = 9;
        set(get(h,'label'),'string','Peak Power');
        ylabel('Channel #')
        set(gca,'XTick',[])
       
        subplot(4,1,3)
        imagesc(combo_array,[0 1])
%         axis equal
%         axis tight
        h = colorbar;
        h.Ticks=[0 1];
        colormap(gca,[0 0 0;1 1 1])
        h.FontSize = 9;
        set(get(h,'label'),'string','Overall Quality');
        ylabel('Channel #')
        set(gca,'XTick',[])
%         xlabel('Window #')
        
        sbplt = subplot(4,1,4);
        tail = length(raw.t) - n_windows*window_samples;
        for ch = 1:size(raw.d,2)
           p(ch) = plot(raw.t(1:end-tail),raw.d(1:end-tail,ch));
           hold on
        end
        hold off
        for cond = 1:size(raw.s,2)
            idx = find(raw.s(:,cond)==1);
            idx(idx>n_windows*window_samples)=[];
            for j=1:size(idx)
                p(size(raw.d,2)+cond) = line([raw.t(idx(j)) raw.t(idx(j))],get(gca,'YLim'));
            end
            hold on
        end
        axis tight
        h = colorbar;
        cbar_pos = get(h,'Position');
        set(h,'Visible','off');
        for cond = 1:size(raw.s,2)
            legend(p(size(raw.d,2)+cond),['cond' num2str(cond)])
        end
        leg_pos = get(legend(sbplt),'Position');
        set(legend(sbplt),'Position',[cbar_pos(1) cbar_pos(2)+cbar_pos(4)-leg_pos(4) leg_pos(3) leg_pos(4)]);
        
        xlabel('Time (s)')
        ylabel('Raw Data')
        
        drawnow
    end
    
    %% Detect artifacts and bad links
    bad_links = find(good_combo_link<0.85);
    bad_windows = find(good_combo_window<0.85);
    
    new_table_row = table(i, {bad_links'}, {bad_windows});
    report_table = [report_table; new_table_row];
        
    if plot_bad
        n_col = ceil(length(bad_links)/5);
        if ~isempty(bad_links)
           figure('Name',['File ' num2str(i) ' - Bad Links'],'NumberTitle','off') 
           for j = 1:length(bad_links)  
               subplot(5,n_col,j)
               hold on
               plot(cardiac_data(1,:,bad_links(j)),'b')
               plot(cardiac_data(2,:,bad_links(j)),'r')
               axis tight
               ylim([-4 4])
           end
        end

        n_col = ceil(length(bad_windows)/5);
        if ~isempty(bad_windows)
           figure('Name',['File ' num2str(i) ' - Bad Windows'],'NumberTitle','off') 
           for j = 1:length(bad_windows)
               subplot(5,n_col,j)
               idx = find(~combo_array(:,bad_windows(1)));
               interval = (bad_windows(j)-1)*window_samples-(bad_windows(j)-1)*(overlap_samples)+1 : bad_windows(j)*window_samples-(bad_windows(j)-1)*(overlap_samples);
               hold on
               plot(cardiac_data(1,interval,idx(1)),'b')
               plot(cardiac_data(2,interval,idx(1)),'r')
               axis tight
               ylim([-4 4])
           end
        end
        close all
    end

end
report_table.Properties.VariableNames = {'file_idx','Bad_Links','Bad_Windows'};

%% Save on Excel file

for i=1:size(report_table,1)
    a = report_table.Bad_Links{i};
    b = report_table.Bad_Windows{i};
    a1 = num2str(a);
    b1 = num2str(b);
    report_table.Bad_Links{i} = a1;
    report_table.Bad_Windows{i} = b1;
end

if save_report_table
    writetable(report_table,'Quality_Report.xls');
end


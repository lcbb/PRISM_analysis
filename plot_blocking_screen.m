function plot_blocking_screen(imgPath,plateInfo)
% plot bar plots for blocking screen results
%%%%%%%%%%%%%%%%%%%
global fig_num
% imgPath = 'E:\data\Neuron\cortical\14days\BlockingScreen_061015';
[uniqueConds, uniq_ind] = unique({plateInfo.comboLabelNoRep});
plateFieldNames = {'PermCond','PermConc','PermTime','BlockCond','BlockConc',...
    'BlockTime','BlueChan','GreenChan','RedChan','FarRedChan','Replicate',...
    'RowLabels','ColLabels'};
x_label = [{plateInfo(uniq_ind).PermConc}', {plateInfo(uniq_ind).PermCond}', {plateInfo(uniq_ind).PermTime}', {plateInfo(uniq_ind).BlockConc}', {plateInfo(uniq_ind).BlockCond}'] ;
block_ind = cellfun(@(x)(isempty(x)), strfind(x_label(:,2), 'Saponin'))& cellfun(@(x)(isempty(x)), strfind(x_label(:,2), 'Digitonin')) ;

none_ind = cellfun(@(x)(~isempty(x)), strfind(x_label(:,4), 'none')) ;
x_label(block_ind, 1:4) = cell(size(x_label(block_ind, 1:4))) ;
x_label(none_ind, 4:5) = cell(size(x_label(none_ind, 4:5))) ;
empty_ind = cellfun(@(x)(isempty(x)), x_label) ;
control_ind = all(empty_ind,2) ;
x_label(control_ind, 2) = {'Control'} ;
fig_path = fullfile(imgPath,'ParsedImgStacks'); 
% fig_num  = 1 ;
%%
for condNum = 1:numel(uniqueConds)
    if ~block_ind(condNum)
        switch x_label{condNum,1}

            case '1'
                x_label{condNum,1} = 'High ' ;
            case '2'
                x_label{condNum,1} = 'Low ' ;

            otherwise
        end
    end
    x_label_all{condNum} = [x_label{condNum,:}] ;
end
%%
bar_ind = 0 ;
%%
for condNum = 1:numel(uniqueConds)    
    file_name = fullfile(imgPath,'ParsedImgStacks',[uniqueConds{condNum},'.mat']) ;
    if ~isempty(strfind(uniqueConds{condNum},'Bassoon'))||~isempty(strfind(uniqueConds{condNum},'Homer')) %if synapsin 1

        disp(['process ' file_name '...'])
        load(file_name);

    %         mkdir(plateInfo.fig_path) ;

    %         ab_params = cell2mat(ab_params_cell) ;
    %         Iab_nuc = cat(1,ab_params.InucMeanNormal) ;
        Iab_nuc = cat(1,ab_params.InucMean) ;
        Iab_syn = cat(1, ab_params.ISynMean) ;
        Imarker_syn = cat(1, marker_params.ISynMean) ;

    %         ab_params_cell = struct2cell(ab_params) ;
    %         Iab_nuc = cell2mat(ab_params_cell(3,:)) ;

    
    
        if ~isempty(strfind(uniqueConds{condNum},'Homer'))
            ind = 1 ;
            plot_ind(condNum) = 1 ;
            bar_ind = bar_ind +1 ;
        else
            ind = 2 ;
            plot_ind(condNum) = 0 ;
        end
        Iab_nuc_syn_rt = Iab_nuc./Iab_syn ;
        
        Iab_nuc_syn_rt_mean(bar_ind,ind) = mean(Iab_nuc_syn_rt) ;
        Iab_nuc_syn_rt_std(bar_ind,ind) = std(Iab_nuc_syn_rt) ;
        Iab_nuc_syn_rt_ste(bar_ind,ind) = std(Iab_nuc_syn_rt)/sqrt(numel(Iab_nuc_syn_rt)) ;
        [~, ~, Ci,~] = normfit(Iab_nuc_syn_rt,0.05) ; 
        Iab_nuc_syn_rt_lci(bar_ind,ind) = Iab_nuc_syn_rt_mean(bar_ind,ind)- Ci(1) ;
        Iab_nuc_syn_rt_hci(bar_ind,ind) = Ci(2)-Iab_nuc_syn_rt_mean(bar_ind,ind) ;  
        
        
        Iab_nuc_mean(bar_ind,ind) = mean(Iab_nuc) ;
        Iab_nuc_std(bar_ind,ind) = std(Iab_nuc) ;
        Iab_nuc_ste(bar_ind,ind) = std(Iab_nuc)/sqrt(numel(Iab_nuc)) ;

        Iab_syn_mean(bar_ind,ind) = mean(Iab_syn) ;
        Iab_syn_std(bar_ind,ind) = std(Iab_syn) ;
        Iab_syn_ste(bar_ind,ind) = std(Iab_syn)/sqrt(numel(Iab_syn)) ;

        Imarker_syn_mean(bar_ind,ind) = mean(Imarker_syn) ;
        Imarker_syn_std(bar_ind,ind) = std(Imarker_syn) ;
        Imarker_syn_ste(bar_ind,ind) = std(Imarker_syn)/sqrt(numel(Imarker_syn)) ;
%             bar_ind = bar_ind +1 ;
    else
        plot_ind(condNum) = 0 ;
    end
end
Iab_nuc_syn_rt_ci = cat(3, Iab_nuc_syn_rt_lci, Iab_nuc_syn_rt_hci) ; 

% [hpt p_val] = ttest2(Dis1,Dis2) ;



%%
figure(1)  
% barmap=parula(5);
[hb he] =  barwitherr_v2(Iab_nuc_ste, Iab_nuc_mean, 1,'LineWidth', 2);    % Plot with errorbars

% barmap=[0.5 0.5  0.5]; %[0.7 0.7 0.7] is grey, [ 0.05 .45 0.1] is green
% colormap(barmap);

set(gca,'XTickLabel',x_label_all(logical(plot_ind)))
% set(gca, 'ylim', [0 1])
  legend('Homer-DNA','bassoon-DNA')

ylabel('Nuclear intensity')
  format_fig2(1)  
pos = get(gca,'Position');
set(gca,'Position',[pos(1), .3, pos(3) .65])
set(gca,'LineWidth', 2)
set(gca,'TickLength'  , [.02 .02])

Xt = [1:size(Iab_nuc_mean,1)];
% 
% If you want to set x-axis limit, uncomment the following two lines of 
% code and remove the third
%Xl = [1 6]; 
%set(gca,'XTick',Xt,'XLim',Xl);
set(gca,'XTick',Xt);
% ensure that each string is of the same length, using leading spaces

ax = axis; % Current axis limits
axis(axis); % Set the axis limit modes (e.g. XLimMode) to manual
Yl = ax(3:4); % Y-axis limits

% Remove the default labels
set(gca,'XTickLabel','')

% Place the text labels
t = text(Xt,Yl(1)*ones(1,length(Xt)),x_label_all(logical(plot_ind)));
set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
'Rotation',45, 'Fontsize', 14);
%%
figure(2)  
%   bar(mp_count);    % Plot with errorbars
[hb he] =  barwitherr_v2(Imarker_syn_ste, Imarker_syn_mean, 1,'LineWidth', 2);    % Plot with errorbars

% barmap=[0.5 0.5  0.5]; %[0.7 0.7 0.7] is grey, [ 0.05 .45 0.1] is green
% colormap(barmap);

set(gca,'XTickLabel',x_label_all(logical(plot_ind)))
% set(gca, 'ylim', [0 1])
  legend('Synapsin-1(SC)','Synapsin-1(MP)')

ylabel('Synaptic intensity')
  format_fig2(1)  
pos = get(gca,'Position');
set(gca,'Position',[pos(1), .3, pos(3) .65])
set(gca,'LineWidth', 2)
set(gca,'TickLength'  , [.02 .02])

Xt = [1:size(Iab_nuc_mean,1)];
% 
% If you want to set x-axis limit, uncomment the following two lines of 
% code and remove the third
%Xl = [1 6]; 
%set(gca,'XTick',Xt,'XLim',Xl);
set(gca,'XTick',Xt);
% ensure that each string is of the same length, using leading spaces

ax = axis; % Current axis limits
axis(axis); % Set the axis limit modes (e.g. XLimMode) to manual
Yl = ax(3:4); % Y-axis limits

% Remove the default labels
set(gca,'XTickLabel','')

% Place the text labels
t = text(Xt,Yl(1)*ones(1,length(Xt)),x_label_all(logical(plot_ind)));
set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
'Rotation',45, 'Fontsize', 14);
%%
figure(3)  
%   bar(mp_count);    % Plot with errorbars
[hb he] =  barwitherr_v2(Iab_syn_ste, Iab_syn_mean, 1,'LineWidth', 2);    % Plot with errorbars

% barmap=[0.5 0.5  0.5]; %[0.7 0.7 0.7] is grey, [ 0.05 .45 0.1] is green
% colormap(barmap);

set(gca,'XTickLabel',x_label_all(logical(plot_ind)))
% set(gca, 'ylim', [0 1])
  legend('Homer-DNA','bassoon-DNA')

ylabel('Synaptic intensity')
  format_fig2(1)  
pos = get(gca,'Position');
set(gca,'Position',[pos(1), .3, pos(3) .65])
set(gca,'LineWidth', 2)
set(gca,'TickLength'  , [.02 .02])

Xt = [1:size(Iab_nuc_mean,1)];
% 
% If you want to set x-axis limit, uncomment the following two lines of 
% code and remove the third
%Xl = [1 6]; 
%set(gca,'XTick',Xt,'XLim',Xl);
set(gca,'XTick',Xt);
% ensure that each string is of the same length, using leading spaces

ax = axis; % Current axis limits
axis(axis); % Set the axis limit modes (e.g. XLimMode) to manual
Yl = ax(3:4); % Y-axis limits

% Remove the default labels
set(gca,'XTickLabel','')

% Place the text labels
t = text(Xt,Yl(1)*ones(1,length(Xt)),x_label_all(logical(plot_ind)));
set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
'Rotation',45, 'Fontsize', 14);
%%
figure(4)  
%   bar(mp_count);    % Plot with errorbars
[hb he] =  barwitherr_v2(Iab_nuc_syn_rt_ci, Iab_nuc_syn_rt_mean, 1,'LineWidth', 2);    % Plot with errorbars

% barmap=[0.5 0.5  0.5]; %[0.7 0.7 0.7] is grey, [ 0.05 .45 0.1] is green
% colormap(barmap);

set(gca,'XTickLabel',x_label_all(logical(plot_ind)))
% set(gca, 'ylim', [0 1])
  legend('Homer-DNA','bassoon-DNA')

ylabel('Nuclear/Synaptic Intensity')
  format_fig2(1)  
pos = get(gca,'Position');
set(gca,'Position',[pos(1), .3, pos(3) .65])
set(gca,'LineWidth', 2)
set(gca,'TickLength'  , [.02 .02])

Xt = [1:size(Iab_nuc_mean,1)];
% 
% If you want to set x-axis limit, uncomment the following two lines of 
% code and remove the third
%Xl = [1 6]; 
%set(gca,'XTick',Xt,'XLim',Xl);
set(gca,'XTick',Xt);
% ensure that each string is of the same length, using leading spaces

ax = axis; % Current axis limits
axis(axis); % Set the axis limit modes (e.g. XLimMode) to manual
Yl = ax(3:4); % Y-axis limits

% Remove the default labels
set(gca,'XTickLabel','')

% Place the text labels
t = text(Xt,Yl(1)*ones(1,length(Xt)),x_label_all(logical(plot_ind)));
set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
'Rotation',45, 'Fontsize', 14);
  %%
%   print(gcf ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ;
print(1 ,'-dpdf','-r300', fullfile(fig_path,[num2str(fig_num, '%03d'),'.pdf'])) ; 
print(1 ,'-dpng','-r300', fullfile(fig_path,[num2str(fig_num, '%03d'),'.png'])) ; fig_num = fig_num +1 ;
print(2 ,'-dpdf','-r300', fullfile(fig_path,[num2str(fig_num, '%03d'),'.pdf'])) ; 
print(2 ,'-dpng','-r300', fullfile(fig_path,[num2str(fig_num, '%03d'),'.png'])) ; fig_num = fig_num +1 ;
print(3 ,'-dpdf','-r300', fullfile(fig_path,[num2str(fig_num, '%03d'),'.pdf'])) ; 
print(3 ,'-dpng','-r300', fullfile(fig_path,[num2str(fig_num, '%03d'),'.png'])) ; fig_num = fig_num +1 ;
print(4 ,'-dpdf','-r300', fullfile(fig_path,[num2str(fig_num, '%03d'),'.pdf'])) ; 
print(4 ,'-dpng','-r300', fullfile(fig_path,[num2str(fig_num, '%03d'),'.png'])) ; fig_num = fig_num +1 ;


%         save(file_name, 'nuc_params_cell', 'ab_params_cell', '-append');

end
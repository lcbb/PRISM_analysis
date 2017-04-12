function cluster_synapse(dir_cell)
% clustering analysis of all synpses from all replicates
global fig_num fig_path work_path;
     IabInMarkerAll = [] ;
     IabWithMarkerAll = [] ;
     AabWithMarkerAll = [] ;
    for k = 1:numel(dir_cell)
        imgCondiFolders = readDirSubfolders(dir_cell{k},'all');
        for j = 1:numel(imgCondiFolders)
            load(fullfile(dir_cell{k}, imgCondiFolders{j}, 'mask.mat'))            
            IabInMarker = cell2mat(IabInMarkerCell') ;
            IabWithMarker = cell2mat(IabWithMarkerCell') ;
            AabWithMarker = cell2mat(AabWithMarkerCell') ;
            IabInMarkerAll = [IabInMarkerAll; IabInMarker] ;
            IabWithMarkerAll = [IabWithMarkerAll; IabWithMarker] ;
            AabWithMarkerAll = [AabWithMarkerAll; AabWithMarker] ;            
        end
    end
    %%
    ch_name = imgRoundNamesNoWO ;
    nonsyn_ind_bi = cellfun(@(x) ~isempty(regexpi(x,'tuj-1|MAP2', 'match')), imgRoundNamesNoWO,'uniformoutput',0) ; 
    nonsyn_ind_bi = cat(1,nonsyn_ind_bi{:}) ;
    marker2nd_ind_bi = cellfun(@(x) ~isempty(regexpi(x,'synapsin-1\(2\)', 'match')), imgRoundNamesNoWO,'uniformoutput',0) ; % exclude the second synapsin-1 
    marker2nd_ind_bi = cat(1,marker2nd_ind_bi{:}) ;
%     DataMat = [IabWithMarkerAll'; AabWithMarkerAll']  ;
    IabWithMarkerAll(:,nonsyn_ind_bi) = IabInMarkerAll(:,nonsyn_ind_bi) ;  % replace intensities of cytoskeletal targets that fail to be segmented properly with intensities measured within synapsin-1 mask  
    DataMat = [IabWithMarkerAll(:,~marker2nd_ind_bi)'; AabWithMarkerAll(:,~nonsyn_ind_bi & ~marker2nd_ind_bi)']  ;
    % DataMat = [IabWithMarker(:,3:end)'; AabWithMarker(:,3:end)']  ;
    % DataMat = zscore(DataMat,0,2) ;
    % DataMat = DataMat./repmat(mean(DataMat,2), [1 size(DataMat,2)]) ;
    rowlabI = cellfun(@(c2)(['I(',c2, ')']), ch_name(~marker2nd_ind_bi),'UniformOutput', false);
    rowlabA = cellfun(@(c2)(['A(',c2, ')']), ch_name(~nonsyn_ind_bi & ~marker2nd_ind_bi),'UniformOutput', false);
    rowlab = [rowlabI; rowlabA] ;
    %% Plot distributions of features 
    figure(103)
    distributionPlot(IabWithMarkerAll,'xyOri','flipped','histOri','right','showMM',5,'color', [0.7 0.7 0.7],'histOpt',2)
    set(gca,'FontSize', 10)
    set(gca,'LineWidth', 2)
    set(gca,'TickLength'  , [.02 .02])
    set(gca,'yticklabel'  ,ch_name)
    xl = xlabel('Intensity (AU)') ; 
    % % ylabel('log(Count)')
    ylabel('Count') 
    format_fig2(5)
    hold off
    print(103 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; 
    print(103 ,'-dpdf','-r300', [fig_path,num2str(fig_num, '%03d'),'.pdf']) ; fig_num = fig_num +1 ; save([work_path, 'startup.mat'], 'fig_num');
        %%
    figure(104)
    distributionPlot(AabWithMarkerAll,'xyOri','flipped','histOri','right','showMM',5,'color', [0.7 0.7 0.7],'histOpt',2)
    set(gca,'FontSize', 10)
    set(gca,'LineWidth', 2)
    set(gca,'TickLength'  , [.02 .02])
    set(gca,'yticklabel'  ,ch_name)
    xl = xlabel('Intensity (AU)') ; 
    % % ylabel('log(Count)')
    ylabel('Count') 
    format_fig2(5)
    hold off
    print(104 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; 
    print(104 ,'-dpdf','-r300', [fig_path,num2str(fig_num, '%03d'),'.pdf']) ; fig_num = fig_num +1 ; save([work_path, 'startup.mat'], 'fig_num');
%%
    hier_clustering(DataMat, rowlab)
end
function hier_clustering(DataMat, rowlab)
    %% clustering analysis
    global fig_num fig_path work_path;
    for j = 1:size(DataMat,1)
        DataMatS = DataMat(j,:) ;
        DataMatS_nonzero =  DataMatS(DataMatS>0) ;
        out_ind = left_outlier_1d(DataMatS_nonzero, -2.5) ;
        DataMatS_nonzero(out_ind) = 0 ;
        disp(sum(out_ind))
        DataMatS(DataMatS>0) = DataMatS_nonzero ;
%         DataMatS =  log(DataMatS+1) ; %log transformation
        DataMatS = DataMatS./std(DataMatS(DataMatS>0)) ; % exclude 0 values from std estimation
%         DataMatS = DataMatS/std(DataMatS) ;
        DataMatS(DataMatS>0) = DataMatS(DataMatS>0) - min(DataMatS(DataMatS>0)) ;
        DataMat(j,:) =  DataMatS ;
%         DataMat(j,:) =  log(DataMatS+1) ; %log transformation
    end    
%%
    figure(105)
    distributionPlot(DataMat','xyOri','flipped','histOri','right','showMM',5,'color', [0.7 0.7 0.7],'histOpt',2)
    set(gca,'FontSize', 10)
    set(gca,'LineWidth', 2)
    set(gca,'TickLength'  , [.02 .02])
    set(gca,'yticklabel'  ,rowlab)
    xl = xlabel('Intensity (AU)') ; 
    % % ylabel('log(Count)')
    ylabel('Count') 
    format_fig2(5)
    hold off
    print(105 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; 
    print(105 ,'-dpdf','-r300', [fig_path,num2str(fig_num, '%03d'),'.pdf']) ; fig_num = fig_num +1 ; save([work_path, 'startup.mat'], 'fig_num');
    %% find optimal number of clusters using silhouette distance
    clusfunc = @(X,K)(cluster(linkage(X,'ward'),'maxclust',K));
    % E = evalclusters(DataMat','linkage','silhouette','klist',[1:size(DataMat,2)]) ;
    % E = evalclusters(DataMat,clusfunc,'silhouette','klist',[1:size(DataMat,1)]) ;
    EvRow = evalclusters(DataMat,clusfunc,'silhouette','klist',[1:20]) ;
    EvCln = evalclusters(DataMat',clusfunc,'silhouette','klist',[1:10]) ;
    %%%
    figure(60)
    hER = plot(EvRow) ;
    set(hER, 'Marker', 'none', 'MarkerSize',10, 'LineWidth',2,'Color', lines(1))
    set(gca,'LineWidth',2, 'ylim', [0 1]) 
    format_fig2(2)

    figure(61)
    hEC = plot(EvCln) ;
    set(hEC, 'Marker', 'none', 'MarkerSize',10, 'LineWidth',2,'Color', lines(1))
    set(gca,'LineWidth',2, 'ylim', [0 1]) 
    format_fig2(2)
    %%%
    print(60 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; 
    print(60 ,'-dpdf','-r300', [fig_path,num2str(fig_num, '%03d'),'.pdf']) ; fig_num = fig_num +1 ; save([work_path, 'startup.mat'], 'fig_num');
    print(61 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; 
    print(61 ,'-dpdf','-r300', [fig_path,num2str(fig_num, '%03d'),'.pdf']) ; fig_num = fig_num +1 ; save([work_path, 'startup.mat'], 'fig_num');
    %%
    treeRow = linkage(DataMat,'ward') ;
    treeCln = linkage(DataMat','ward') ;
    DenThr(1) = treeRow(end-EvRow.OptimalK+2,3)-eps;
    DenThr(2) = treeCln(end-EvCln.OptimalK+2,3)-eps;
    % DenThr(2) = treeCln(end-2+2,3)-eps;
%     DenThr(2) = 45;
    %%
    cgo = clustergram(DataMat,'Linkage','ward') ;
    %%%
    cmap = french(512, 3) ;
    cmap = cmap(end:-1:1,:) ;
    cmapr = cmap(257:end,:) ;
    set(cgo,'RowLabels',rowlab, 'Colormap', cmapr, 'Dendrogram', DenThr,'Symmetric', 'false','DisplayRange', 6)
    % set(cgo,'RowLabels',rowlab, 'Colormap', cmap, 'DisplayRange', 3)
    % figure(60)
    h_cgo = plot(cgo) ;
    %%
    print(gcf ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; 
    print(gcf ,'-dpdf','-r300', [fig_path,num2str(fig_num, '%03d'),'.pdf']) ; fig_num = fig_num +1 ; save([work_path, 'startup.mat'], 'fig_num');
    % print(gcf ,'-depsc','-r300', [fig_path,num2str(fig_num, '%03d'),'.eps']) ;
    % addYLabel(cgo , rowlab, 'FontSize', 8)

end

function out_ind = left_outlier_1d(x, thresh)
% find the outliers that are on the left side of the input 1D distribution only
x_med = median(x) ;
x_dev = x-x_med ;
x_abs_dev = abs(x_dev) ;
x_abs_dev_med = median(x_abs_dev) ;
z_score_modi = 0.6745 * x_dev/x_abs_dev_med ;
out_ind = z_score_modi < thresh ;
end
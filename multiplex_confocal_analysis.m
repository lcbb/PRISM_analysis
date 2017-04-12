function multiplex_confocal_analysis(dir_cell)
global fig_num fig_path work_path;
RepNameCell = {'Rep1', 'Rep2', 'Rep3'} ;
pccCellName = {'pccIabInMarkerMean', 'pccIabWithMarkerMean', 'pccAabWithMarkerMean'} ;
net_thresh = [0.35:0.01:0.4] ; % threshold for showing connections in the network representation
for j = 1:numel(dir_cell)    
    pccIabInMarker = [] ;
    pccIabWithMarker = [] ; 
    pccAabWithMarker = [] ;
    IabInMarkerMean = [] ;            
    IabInMarker = [] ;
    IabWithMarkerMean = [] ;            
    IabWithMarker = [] ;
    AabWithMarkerMean = [];            
    AabWithMarker = [] ;
    Nsynapse = [] ;
    
    imgCondiFolders = readDirSubfolders(dir_cell{j},'all');
    
    for k = 1:numel(imgCondiFolders)    
        load(fullfile(dir_cell{j}, imgCondiFolders{k}, 'mask.mat'))
        ch_name = imgRoundNamesNoWO ;
        pccIabInMarker = cat(3, pccIabInMarker, pccIabInMarkerAll) ;
        pccIabWithMarker = cat(3, pccIabWithMarker, pccIabWithMarkerAll) ; 
        pccAabWithMarker = cat(3, pccAabWithMarker, pccAabWithMarkerAll) ; 
        
        IabInMarkerMeanCell = cellfun(@(x)mean(x,1), IabInMarkerCell, 'uniformoutput',0) ;
        IabInMarkerMean = cat(1,IabInMarkerMean, cat(1,IabInMarkerMeanCell{:})) ;            
        IabInMarker = cat(1,IabInMarker, cat(1,IabInMarkerCell{:})) ;
        NsynapseCell = cellfun(@(x)size(x,1), IabInMarkerCell, 'uniformoutput',0) ;
        Nsynapse = cat(1,Nsynapse, cat(1, NsynapseCell{:})) ;

        IabWithMarkerMeanCell = cellfun(@(x)nonzeromean(x,1), IabWithMarkerCell, 'uniformoutput',0) ;
%         IabWithMarkerMeanCell = cellfun(@(x)mean(x,1), IabWithMarkerCell, 'uniformoutput',0) ;        
        IabWithMarkerMean = cat(1,IabWithMarkerMean, cat(1,IabWithMarkerMeanCell{:})) ;            
        IabWithMarker = cat(1,IabWithMarker, cat(1,IabWithMarkerCell{:})) ;

        AabWithMarkerMeanCell = cellfun(@(x)nonzeromean(x,1), AabWithMarkerCell, 'uniformoutput',0) ;
        AabWithMarkerMean = cat(1,AabWithMarkerMean, cat(1,AabWithMarkerMeanCell{:})) ;            
        AabWithMarker = cat(1,AabWithMarker, cat(1,AabWithMarkerCell{:})) ;
        
    end
    [pccIabInMarkerMean, pccIabInMarkerStd, pccIabInMarkerSte,  pccIabInMarkerCi] = pcc_stats(pccIabInMarker) ;
    [pccIabWithMarkerMean, pccIabWithMarkerStd, pccIabWithMarkerSte,  pccIabWithMarkerCi] = pcc_stats(pccIabWithMarker) ;
    [pccAabWithMarkerMean, pccAabWithMarkerStd, pccAabWithMarkerSte,  pccAabWithMarkerCi] = pcc_stats(pccAabWithMarker) ;
    pccCell(j,:) = {pccIabInMarkerMean, pccIabWithMarkerMean, pccAabWithMarkerMean} ;
    pccCiCell(j,:) = {pccIabInMarkerCi, pccIabWithMarkerCi, pccAabWithMarkerCi} ;
    IabInMarkerMeanAll{j} = IabInMarkerMean ;
    IabInMarkerAll{j} = IabInMarker ;
    IabWithMarkerMeanAll{j} = IabWithMarkerMean ;
    IabWithMarkerAll{j} = IabWithMarker ;
    AabWithMarkerMeanAll{j} = AabWithMarkerMean ;
    AabWithMarkerAll{j} = AabWithMarker ;
    
    %%
    tlt = ['Correlation of synaptic intensity in synapsin (' RepNameCell{j} ')'];
    cmap = french(512, 2) ;
    cmap = cmap(end:-1:1,:) ;
    triangleHeatMap(pccIabInMarkerMean, pccIabInMarkerCi, tlt,ch_name,ch_name, cmap, 0)
    % triangleHeatMap(pccIabInMarker, pccIabInMarkerCi, tlt,ch_name,ch_name, cmap, 0)
    % end
    tlt = ['Correlation of synaptic intensity with synapsin (' RepNameCell{j} ')'];
    triangleHeatMap(pccIabWithMarkerMean, pccIabWithMarkerCi, tlt,ch_name,ch_name, cmap, 0)
    tlt = ['Correlation of synapse size (' RepNameCell{j} ')'];
    triangleHeatMap(pccAabWithMarkerMean, pccAabWithMarkerCi, tlt,ch_name,ch_name, cmap, 0)    

    marker2nd_ind_bi = cellfun(@(x) ~isempty(regexpi(x,'\(2\)', 'match')), imgRoundNamesNoWO,'uniformoutput',0) ; 
    marker2nd_ind_bi = cat(1,marker2nd_ind_bi{:}) ;
%     marker2nd_ind = find(marker2nd_ind_bi, 1, 'first') ; % pick the first marker channel if there are more than 1
    
    cd(params.outputImgsPath)
    for m = 1: numel(pccCellName)               
        generateGraphMLInput.adjMtx = pccCell{j,m};
        generateGraphMLInput.adjMtx(marker2nd_ind_bi,:) =  [] ; % delete the duplicated synapsin channel
        generateGraphMLInput.adjMtx(:,marker2nd_ind_bi) =  [] ;
        generateGraphMLInput.adjMtx = generateGraphMLInput.adjMtx - net_thresh(1) ;
        generateGraphMLInput.adjMtx(generateGraphMLInput.adjMtx<0)=0 ; % scale for better visualization         
        generateGraphMLInput.nodeLabels = ch_name(~marker2nd_ind_bi);
        generateGraphMLInput.filename = [params.Condi, '_' , pccCellName{m},'_thresh=0o', num2str(net_thresh(1))] ; 
%         adjMtx2GraphMLFixedNodes(generateGraphMLInput)
        adjMtx2GraphML(generateGraphMLInput)
    end        
end
disp([num2str(mean(Nsynapse)), '+-' ,num2str(std(Nsynapse)), ' synapsese per field of view'])

%% network of mean correlation matrices
for j = 1:numel(net_thresh)
    for m = 1: numel(pccCellName)      
        pccCellMean = mean(cat(3,pccCell{:,m}),3) ;
        generateGraphMLInput.adjMtx = pccCellMean;
        generateGraphMLInput.adjMtx(marker2nd_ind_bi,:) =  [] ; % delete the duplicated synapsin channel
        generateGraphMLInput.adjMtx(:,marker2nd_ind_bi) =  [] ;
        generateGraphMLInput.adjMtx = generateGraphMLInput.adjMtx - net_thresh(j) ;
        generateGraphMLInput.adjMtx(generateGraphMLInput.adjMtx<0)=0 ; % scale for better visualization         
        generateGraphMLInput.nodeLabels = ch_name(~marker2nd_ind_bi);
        generateGraphMLInput.filename = ['Mean_' , pccCellName{m}, '_thresh=0o', num2str(net_thresh(j)) ] ; 
%         adjMtx2GraphMLFixedNodes(generateGraphMLInput)
        adjMtx2GraphML(generateGraphMLInput)
    end
end


%%
GroupDistributionPlot(IabInMarkerAll, ch_name, 'synaptic intensity in synapsin (A.U.)', RepNameCell); close(2)
GroupDistributionPlot(IabWithMarkerAll, ch_name, 'synaptic intensity with synapsin (A.U.)', RepNameCell); close(2)
GroupDistributionPlot(AabWithMarkerAll, ch_name, 'synaptic size (pixels)', RepNameCell); close(2)
%%
GroupDistributionPlot(IabInMarkerMeanAll, ch_name, 'synaptic intensity in synapsin (A.U.)', RepNameCell); close(2)
GroupDistributionPlot(IabWithMarkerMeanAll, ch_name, 'synaptic intensity with synapsin (A.U.)', RepNameCell); close(2)
GroupDistributionPlot(AabWithMarkerMeanAll, ch_name, 'synaptic size (pixels)', RepNameCell); 


%% scattered plots for intensities 
IabInMarkerAll = cat(1,IabInMarkerCell{:}) ;
IabLabelAll = ones(size(IabInMarkerAll,1),1) ;

    figure(50)
    [h,ax,bigax] = gplotmatrix(IabInMarkerAll,[],IabLabelAll,...
        [],[],[1 1],'off','hist',ch_name, ch_name);
        set(ax, 'TickLabelInterpreter','none')
        set(ax, 'FontSize',6)
    print(50 ,'-dpng','-r600', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ; save([work_path, 'startup.mat'], 'fig_num');
    close(50)
%%
cmap = french(512, 1) ;
cmap = cmap(end:-1:1,:) ;
pccCellDiff =cellfun(@(a,b)(b-a), pccCell(1,:), pccCell(2,:), 'uniformoutput',0) ;
tlt = ['Difference in correlation of synaptic intensity(',RepNameCell{2}, '-' RepNameCell{1}, ')' ]  ;
triangleHeatMap(pccCellDiff{1}, [], tlt,ch_name,ch_name, cmap, 0)
%%
pccCellDiff =cellfun(@(a,b)(b-a), pccCell(1,:), pccCell(3,:), 'uniformoutput',0) ;
tlt = ['Difference in correlation of synaptic intensity(',RepNameCell{3}, '-' RepNameCell{1}, ')' ]  ;
triangleHeatMap(pccCellDiff{1}, [], tlt,ch_name,ch_name, cmap, 0)

end
function [pccIabInMarkerMean, pccIabInMarkerStd, pccIabInMarkerSte,  pccIabInMarkerCi] = pcc_stats(pccIabInMarkerAll)
    pccIabInMarkerMean = mean(pccIabInMarkerAll,3) ;
    pccIabInMarkerStd = std(pccIabInMarkerAll,0,3) ;
    pccIabInMarkerSte = pccIabInMarkerStd/sqrt(size(pccIabInMarkerAll,3)) ;% Standard Error          
    ts = tinv([0.025  0.975],size(pccIabInMarkerAll,3)-1);      % T-Score
    pccIabInMarkerCi = ts(2)*pccIabInMarkerSte;                      % Confidence Intervals
end

function  [DStats, markerPairFrc] = synapse_mask_analysis_confocal(imgsProjAligned, markerMaskCell, abMaskCell, params, plot_ind, ch_name, int_opt) 
% int_opt: 1 or 0, 1 to re-intergrate the synaptic intensity within the circles centered at synapsin-1 punta

global fig_num fig_path work_path;
close all
maskfile = fullfile(params.outputImgsPath,params.Condi, 'mask.mat') ;
load(maskfile);
ch_name = {'Synapsin-1', 'ARPC2','PSD95','NMDAR', 'Homer-1bc','Tuj-1','Actin','Bassoon','Shank3','Synapsin-1(2)','Cortactin'} ; %rename the channel name
marker_ind = params.marker_ind ;
ab_ind_all = params.ab_ind_all ;
uniqueFields = [1:numel(imgsProjAligned)] ;
pccIabInMarkerAll = [] ;
pccIabWithMarkerAll = [] ;
pccAabWithMarkerAll = [] ;
% for f=2:2
for f = 1:numel(uniqueFields)
    disp(['process field ', num2str(f)])
        sizeI = size(imgsProjAligned{f}{1,1});
        Idapi = imgsProjAligned{f}{2,1};
        IMAP2 = imgsProjAligned{f}{1,1};
        ItargetStack = cat(3,imgsProjAligned{f}{3,plot_ind});
    %     ItargetStack = cat(3,ItargetStack,imgsProjAligned{f}{[2,1],1});
        if params.tophat 
            for i = 1:size(ItargetStack,3)
                ItargetStack(:,:,i) = imtophat(ItargetStack(:,:,i),strel('disk',params.disksize));
            end
        end        
ItargetStack = ItargetStack(:,:,plot_ind) ;    
% ItargetStack_rgb = imstack2RGB(ItargetStack , [0.5 0.995]) ;
Imarker = ItargetStack(:,:,marker_ind) ;
pixsize = 0.18733333 ; % check this later

syn_img_size_um = [3 3] ; % size of the the cropped image size in unit of micron, note the circle for synaptic intensity integration is half of this size 
syn_img_size = round(syn_img_size_um/pixsize) ;
delta_r = [0:syn_img_size(1)-1]*pixsize ;
rccf_base_ref = 2; %unit of um
tlt_name = maskfile(1:end-4) ;
os.plot = 0 ;
   
%%
markerMask = markerMaskCell{f} ;
abMaskAll =  abMaskCell{f} ;
params.minObjSize = 12 ;
markerMask = bwareaopen(markerMask, params.minObjSize); % remove small objects 
markerMask = imclearborder(markerMask) ;
% markerStats = regionprops(markerMask,'Area', 'Centroid');
markerStats = regionprops(markerMask,Imarker,'Area', 'WeightedCentroid','MeanIntensity');
% markerCen = cat(1, markerStats.Centroid);
markerCen = cat(1, markerStats.WeightedCentroid);
markerArea = cat(1, markerStats.Area);
nan_ind = any(isnan(markerCen),2) ;
markerCen(nan_ind, :) = [] ;
markerStats(nan_ind) = [] ;
markerArea(nan_ind) = [] ; 
markerAbmask  = zeros([size(markerMask), 3], 'uint8') ;% show overlay of marker and ab masks after distance filtering 
markerAbmask_raw  = zeros([size(markerMask), 3], 'uint8') ;% % show overlay of marker and ab masks 
markerAbmask(:,:,1) = uint8(markerMask*255) ;
markerAbmask_raw(:,:,1) = uint8(markerMask*255) ;
abStatsAll = cell(1, numel(ab_ind_all)) ;

%%
n_obj = size(markerCen,1) ;

% figure(41)
% imshow(markerMask,[], 'InitialMagnification', 'fit');
% tlt_name = ch_name{marker_ind} ;  
% ht = title(tlt_name,'FontSize',20) ;
% hold on
% 
% 
% for j =1:n_obj
%     plot(markerCen(j,1), markerCen(j,2), 'r+')
%     text(markerCen(j,1)+6,markerCen(j,2)+6,num2str(j),'FontSize',10,'color','r')   
% end
% hold off
% print(41 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ; save([work_path, 'startup.mat'], 'fig_num');

%% integrate intensity within the circle centered at the synapsin puncta for all channels
if int_opt ==1
    % ind = 0 ;
    IabInMarker = zeros(n_obj, numel(plot_ind)) ;
    IabInMarkerTtl = IabInMarker ;
    cir_size = 0.5*syn_img_size ;
%     parpool(20)
    for k = 1:numel(plot_ind)
        disp(k)
    %     ind = ind+1 ;

    %     abMask = logical(abMaskAll(:,:,k)) ;
        Iab = ItargetStack(:,:,k);
        AbInMarkerStats = regionprops(markerMask,Iab, 'MeanIntensity');
        IabInMarker(:,k) =  cat(1, AbInMarkerStats.MeanIntensity) ;        
        IabInMarkerTtl(:,k) = IabInMarker(:,k).*markerArea ;
%         parfor j = 1:n_obj
%             disp(j)
%             cen = markerCen(j,2:-1:1) ;%be careful about the reverse oder here!!         
%             cir_pos = markerCen(j,:)-round(0.5*cir_size) ;
%             cir_pos = [cir_pos, round(cir_size)] ;        
%             figure(410+j)
%             imshow(markerMask,[], 'InitialMagnification', 'fit');
%             hEllipse = imellipse(gca,cir_pos); % Second argument defines ellipse shape and position.
%             cir_mask = hEllipse.createMask(); % Create a binary image ("mask") from the ROI object.
%             close(410+j)
%     %         imshow(cir_mask,[], 'InitialMagnification', 'fit');
%     %         hold on
%     %         plot(markerCen(j,1), markerCen(j,2), 'r+')
%     %         text(markerCen(j,1)+6,markerCen(j,2)+6,num2str(j),'FontSize',10,'color','r') 
%     %         hold off
%             IabInMarker(j,k) = sum(Iab(cir_mask)) ;
%         end                              
    end
%     delete(gcp)
    IabInMarkerCell{f}= IabInMarker ;     
else
end
%% puncta distantace to marker, colocalization fraction, puncta intensity and area
cut_off_v = [0:50:1000] ;
markerPairFrcV = zeros(1,numel(cut_off_v))  ;
abPairFrcV = zeros(1,numel(cut_off_v))  ;
markerPairIndAll = [] ;
IabWithMarker = zeros(size(markerStats,1), numel(plot_ind)) ;
AabWithMarker = zeros(size(markerStats,1), numel(plot_ind)) ;
ImarkerSyn = cat(1,markerStats.MeanIntensity) ;
AmarkerSyn = cat(1,markerStats.Area) ;
% for j = 1:numel(cut_off_v)
clear DStats
for k = 1:numel(ab_ind_all)
% for k = 1:1
    %%%
    ab_ind = ab_ind_all(k) ;   
    abMask = logical(abMaskAll(:,:,k)) ;
    Iab = ItargetStack(:,:,ab_ind); 

%%%   
    abStats = regionprops(abMask, Iab, 'Area', 'WeightedCentroid', 'MeanIntensity');      
    abCen = cat(1, abStats.WeightedCentroid);   
    
%     [IDX,D]=knnsearch(markerCen,abCen,'k',5);    
    [IDX,D]=knnsearch(abCen,markerCen,'k',5);    
    % bin_l = [0:10:max(trj_roi_cell{j,k}(:,1))] ;
    %
    n_bin = 40 ;
    D = D*pixsize*1000; %convert distance to nm
    %%%
    figure(50) % show histograms of distances to neareast puncta
    h1 = histogram(D(:,1),n_bin) ;
%     xlabel('D(nm)')
%     ylabel('count')
%     format_fig2(2)
   hold all
    
%     figure(51) % show histograms of distances to the second neareast puncta
    h2 = histogram(D(:,2),n_bin) ;
    xlabel('D(nm)')
    ylabel('count')
    legend('D1', 'D2')
    format_fig2(2)
    hold off
%     print(50 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ; save([work_path, 'startup.mat'], 'fig_num');
%%%       
%     pairIndBi = D(:,1)<cut_off_v(j) ; % set cut-off distance 400 nm to define colocolization
    markerPairIndBi = D(:,1)< 1000; % set cut-off distance 400 nm to define colocolization
    markerUnpairIndBi = ~markerPairIndBi ;
%     markerPairInd = find(markerPairIndBi) ;
    markerPairInd = IDX(:,1) ; % For each marker punctum the index of paired ab punctum. Same size as markerCenUnpaired. Unpaired maker puncta are assgined "0"  
     markerPairInd(markerUnpairIndBi) = 0  ;
%     markerUnpairInd = find(markerUnpairIndBi) ;    
    markerCenPair = markerCen(markerPairIndBi,:) ;
    
%     ImarkerSyn = ImarkerSyn(markerPairIndBi) ;
    
    
%%%
    abPairInd = unique(IDX(markerPairIndBi,1)) ;
%     abStats.pairInd = abPairInd ;
    abPairIndNonUni = IDX(markerPairIndBi,1) ; % non-unique ab paired indices
    abCenPair = abCen(abPairIndNonUni,:) ; 
    abPairFrc(k) = numel(abPairInd)/size(abCen,1) ;
    markerPairFrc(k) = sum(markerPairIndBi)/size(markerCen,1) ;
    %%%
    % intensity correlation
%     IabSyn = cat(1,abStats.MeanIntensity) ; 
%     IabSyn = IabSyn(abPairIndNonUni) ;
    IabSyn = markerPairInd ;
    IabSynAll = cat(1,abStats.MeanIntensity) ;
%     IabSynAll = IabSynAll.* cat(1,abStats.Area) ;
    IabSyn(IabSyn>0) = IabSynAll(abPairIndNonUni) ; % assign zero intensity to unpaired puncta
    IabWithMarker(:,ab_ind) = IabSyn ;
    
    AabSyn = markerPairInd ;
    AabSynAll = cat(1,abStats.Area) ;    
    AabSyn(AabSyn>0) = AabSynAll(abPairIndNonUni) ; % assign zero intensity to unpaired puncta
    AabWithMarker(:,ab_ind) = AabSyn ;
 
%%%    
%     [DStats.mean(k) DStats.median(k) DStats.std(k) DStats.lq(k) DStats.hq(k)] = map_stac(D(pairInd)) ;
    [DStats.mean(k), DStats.std(k), Ci,~] = normfit(D(markerPairIndBi,1),0.05) ; %
    DStats.meanCi(k) = (Ci(2) -Ci(1))/2 ;
%     abInd = pairInd ;
    markerPairIndAll = [markerPairIndAll, markerPairInd] ;
   
    
    abLbl = bwlabel(abMask) ;    
    abPairMask = ismember(abLbl, abPairInd) ; % remove puncta that are not paired with marker
    
    abStatsAll{k} = abStats;
%%%
    markerAbmask(:,:,2) = uint8(abPairMask*255) ;
    markerAbmask(:,:,3) = 0 ;
    markerAbmask_raw(:,:,2) = uint8(abMask*255) ;  
    markerAbmask_raw(:,:,3) = 0 ;  
        
    figure(52)
    imshow(markerAbmask,[], 'InitialMagnification', 'fit');
    tlt_name = ['{\color{red}' ch_name{marker_ind}  '}, ' '{\color{green}' ch_name{ab_ind} '}'] ;  
    ht = title(tlt_name,'FontSize',20) ;
%     print(52 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ; save([work_path, 'startup.mat'], 'fig_num');

    figure(53)
    imshow(markerAbmask_raw,[], 'InitialMagnification', 'fit');    
    hold on
    
%     for j =1:size(abCenPair,1)
% %         plot(abCenPair(j,1), abCenPair(j,2), 'g+')
%         line([abCenPair(j,1) markerCenPair(j,1)],[abCenPair(j,2) markerCenPair(j,2)],...
%             'Marker','+', 'LineStyle','-', 'color','w')
%     end
        
    tlt_name = ['{\color{red}' ch_name{marker_ind}  '}, ' '{\color{green}' ch_name{ab_ind} '}'] ;  
    ht = title(tlt_name,'FontSize',20) ;
%     print(53 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ; save([work_path, 'startup.mat'], 'fig_num');
    hold off
%%%  
%     markerPairFrcV(j) = markerPairFrc(k)  ;
%     abPairFrcV(j) = abPairFrc(k)  ;
%     a = mean(IabSyn)
%     b = mean(Iab(:))
end
%%%

markerPairIndAllCell = mat2cell(markerPairIndAll, ones(1,size(markerPairIndAll,1)), size(markerPairIndAll,2)) ;
 [markerStats.pairInd] = markerPairIndAllCell{:} ;
%% 
IabWithMarker(:,marker_ind) = ImarkerSyn ;
AabWithMarker(:,marker_ind) = AmarkerSyn ;
IabWithMarkerCell{f}= IabWithMarker ;  
AabWithMarkerCell{f}= AabWithMarker ;  
[pccIabWithMarker, pccIabWithMarkerPval, pccIabWithMarkerLc, pccIabWithMarkerUc]  = corrcoef(IabWithMarker) ;
[pccAabWithMarker, pccAabWithMarkerPval, pccAabWithMarkerLc, pccAabWithMarkerUc]  = corrcoef(AabWithMarker) ;
pccIabWithMarkerAll = cat(3, pccIabWithMarkerAll, pccIabWithMarker) ;
pccAabWithMarkerAll = cat(3, pccAabWithMarkerAll, pccAabWithMarker) ;
%% Heat Map of PCC of integrated intensity + 95% confidence interval 
[pccIabInMarker, pccIabInMarkerPval, pccIabInMarkerLc, pccIabInMarkerUc]  = corrcoef(IabInMarker) ;
pccIabInMarkerCi = (pccIabInMarkerUc - pccIabInMarkerLc)/2 ;
pccIabInMarkerAll = cat(3, pccIabInMarkerAll, pccIabInMarker) ;
end
%%
save(maskfile, 'IabInMarkerCell','IabWithMarkerCell','AabWithMarkerCell', '-append');
%%
pccIabInMarkerMean = mean(pccIabInMarkerAll,3) ;
pccIabInMarkerStd = std(pccIabInMarkerAll,0,3) ;
pccIabInMarkerSte = pccIabInMarkerStd/sqrt(size(pccIabInMarkerAll,3)) ;% Standard Error          
ts = tinv([0.025  0.975],size(pccIabInMarkerAll,3)-1);      % T-Score
pccIabInMarkerCi = ts(2)*pccIabInMarkerSte;                      % Confidence Intervals

pccIabWithMarkerMean = mean(pccIabWithMarkerAll,3) ;
pccIabWithMarkerStd = std(pccIabWithMarkerAll,0,3) ;
pccIabWithMarkerSte = pccIabWithMarkerStd/sqrt(size(pccIabWithMarkerAll,3)) ;% Standard Error          
ts = tinv([0.025  0.975],size(pccIabWithMarkerAll,3)-1);      % T-Score
pccIabWithMarkerCi = ts(2)*pccIabWithMarkerSte;                      % Confidence Intervals

pccAabWithMarkerMean = mean(pccAabWithMarkerAll,3) ;
pccAabWithMarkerStd = std(pccAabWithMarkerAll,0,3) ;
pccAabWithMarkerSte = pccAabWithMarkerStd/sqrt(size(pccAabWithMarkerAll,3)) ;% Standard Error          
ts = tinv([0.025  0.975],size(pccAabWithMarkerAll,3)-1);      % T-Score
pccAabWithMarkerCi = ts(2)*pccAabWithMarkerSte;                      % Confidence Intervals
save(maskfile, 'pccIabInMarkerMean','pccIabWithMarkerMean','pccAabWithMarkerMean',...
    'pccIabInMarkerCi','pccIabWithMarkerCi','pccAabWithMarkerCi', '-append');

%%
% figure(60)
% plot(cut_off_v, markerPairFrcV, 'linewidth', 2)
% hold all
% plot(cut_off_v, abPairFrcV, 'linewidth', 2)
% xlabel('Cut-off distance (nm)')
% ylabel('Colocalization fraction')
% legend('C_{A}', 'C_{B}')
% tlt_name = ['A:', ch_name{marker_ind}, ', B:', ch_name{ab_ind} ] ;  
% ht = title(tlt_name,'FontSize',20) ;
% % set(gca, 'xscale', 'log')
% format_fig2(2)
% print(60 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ; save([work_path, 'startup.mat'], 'fig_num');

%%
DmeanSymCell = mat2cell(DStats.mean, ones(1,size(DStats.mean,1)),  ones(1,size(DStats.mean,2))) ;
DmeanCiSymCell = mat2cell(DStats.meanCi, ones(1,size(DStats.meanCi,1)),  ones(1,size(DStats.meanCi,2))) ;

DmeanSymCell = cellfun(@(x)num2str(x,'%0.0f'),DmeanSymCell,'UniformOutput', false);
DmeanCiSymCell = cellfun(@(x)num2str(x,'%0.0f'),DmeanCiSymCell,'UniformOutput', false);
% DstdSymCell = cellfun(@(x)num2str(x,'%0.2f'),DstdSymCell,'UniformOutput', false);
% DmeanSymCell = cellfun(@(c1,c2)sprintf([c1, '\n¡À', c2]), DmeanSymCell, DstdSymCell,'UniformOutput', false);
DmeanSymCell = cellfun(@(c1,c2)sprintf([c1, '\n¡À', c2]), DmeanSymCell, DmeanCiSymCell,'UniformOutput', false);
%%%
cmin = min(DStats.mean(:)) ;
cmax = max(DStats.mean(:)) ;

%%
figure(54)
heatmap(DStats.mean, ch_name(ab_ind_all), ch_name(marker_ind), DmeanSymCell,...
    'TickAngle', 45,'TickFontsize',15,'Fontsize', 15,'Colormap', cmapr,...
    'Textcolor', [0 0 0],'MinColorValue', cmin, 'MaxColorValue',cmax,'TickTexInterpreter',1) ;% 'RowLabelsRotate', 45, 'ColumnLabelsRotate', 45) ;
% colormap summer
axis image
colorbar
title('Distance to synapsin-1 (nm)', 'Fontsize', 20)
print(54 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ; save([work_path, 'startup.mat'], 'fig_num');
%%%
figure(55)
heatmap(markerPairFrc, ch_name(ab_ind_all), ch_name(marker_ind),'%0.2f',...
     'TickTexInterpreter',1, 'Colormap', cmapr,'MinColorValue', 0, 'MaxColorValue',1,...
     'TickAngle', 45,'TickFontsize',15,'Fontsize', 15) ;% 'RowLabelsRotate', 45, 'ColumnLabelsRotate', 45) ;
% colormap summer
axis image
colorbar
title('Fraction of colocalized synapsin-1 puncta', 'Fontsize', 20)
print(55 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ;
print(55 ,'-depsc','-r300', [fig_path,num2str(fig_num, '%03d'),'.eps']) ; fig_num = fig_num +1 ; save([work_path, 'startup.mat'], 'fig_num');
%%%
figure(56)
heatmap(abPairFrc, ch_name(ab_ind_all), ch_name(marker_ind),'%0.2f',...
     'TickTexInterpreter',1, 'Colormap', cmapr,'MinColorValue', 0, 'MaxColorValue',1,...
     'TickAngle', 45,'TickFontsize',15,'Fontsize', 15) ;% 'RowLabelsRotate', 45, 'ColumnLabelsRotate', 45) ;
% colormap summer
axis image
colorbar
title('Fraction of puncta colocalized with synapsin-1', 'Fontsize', 20)
print(56 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ;
print(56 ,'-depsc','-r300', [fig_path,num2str(fig_num, '%03d'),'.eps']) ; fig_num = fig_num +1 ; save([work_path, 'startup.mat'], 'fig_num');


%%
abAbPairFrc = zeros(numel(ab_ind_all)) ;
DStats.mean = zeros(numel(ab_ind_all)) ;
DStats.meanCi = zeros(numel(ab_ind_all)) ;
DStats.median = zeros(numel(ab_ind_all)) ;
DStats.std = zeros(numel(ab_ind_all)) ;
DStats.lq = zeros(numel(ab_ind_all)) ;
DStats.hq = zeros(numel(ab_ind_all)) ;
for k = 1:numel(ab_ind_all)
% for k = 4:4
    for m = 1:numel(ab_ind_all)
       
%         ab_ind1 = k ;
%         ab_ind2 = m ;
        ab_ind1 = ab_ind_all(k) ;
        ab_ind2 = ab_ind_all(m) ;   
        abMask1 = logical(abMaskAll(:,:,k)) ;
        abMask2 = logical(abMaskAll(:,:,m)) ;
        abLbl1 = bwlabel(abMask1) ;
        abLbl2 = bwlabel(abMask2) ; 
        abStats1 = abStatsAll{k} ;
        abStats2 = abStatsAll{m} ;
%         pairIndBi1 = cat(1,abStats1.pairIndBi) ;
%         pairIndBi2 = cat(1,abStats2.pairIndBi) ;
        markerPairInd = cat(1, markerStats.pairInd) ;
        markerPairInd1 = markerPairInd(:,k) ; % ab indices paired with each marker
        markerPairInd2 = markerPairInd(:,m) ;
%%%        
%         pairInd1 = abStats1.pairInd ;
%         pairInd2 = abStats2.pairInd ;
        pairInd1 = markerPairInd1(markerPairInd1>0) ;
        pairInd2 =  markerPairInd2(markerPairInd2>0) ; 
%         pairMask1 = ismember(abLbl1, pairInd1) ; % remove puncta that are not paired with marker
%         pairMask2 = ismember(abLbl2, pairInd2) ; % remove puncta that are not paired with marker
        abCen1 = cat(1, abStats1.WeightedCentroid) ;
%         abCen1 = abCen1(pairInd1,:) ;
        abCen2 = cat(1, abStats2.WeightedCentroid) ;
%         abCen2 = abCen2(pairInd2,:) ;
        abPairInd1 = [] ;
        abPairInd2 = [] ;
        abPairIndBi = (markerPairInd1>0)&(markerPairInd2>0) ; 
        abPairInd1 = markerPairInd1(abPairIndBi) ;    
        abPairInd2 = markerPairInd2(abPairIndBi) ;   
        abCenPair1 = abCen1(abPairInd1,:) ;
        abCenPair2 = abCen2(abPairInd2,:) ;
        
%         if isempty(abCen1)||isempty(abCen2)
        if sum(abPairIndBi) == 0
            abAbPairFrc(k,m) = 0 ;
            DStats.mean(k,m) = 0 ;
            DStats.meanCi(k,m) = 0 ;
        else
 
%             [IDX,D]=knnsearch(abCen1,abCen2,'k',5);
%             D = D*os.um_per_px*1000 ;
%             abPairIndBi2 = D(:,1)<400 ;
%             abPairIndNonUni1 = IDX(abPairIndBi2,1) ;       
%             abPairInd1 = unique(IDX(abPairIndBi2,1)) ;        
%     %         abPairInd1 = find(abPairIndBi1) ;
            
            

        %     abAbPairFrc(k,m) = numel(markerPairInd)/size(markerCen,1) ;     
            abAbPairFrc(k,m) = sum(abPairIndBi)/sum(markerPairInd2>0) ;
            D = norm(abCenPair1-abCenPair2) ;
    %         [DStats.mean(k,m), DStats.median(k,m), DStats.std(k,m), DStats.lq(k,m), DStats.hq(k,m)] = map_stac(D(abPairInd2)) ;
            [DStats.mean(k,m), DStats.std(k,m), Ci,~] = normfit(D,0.05) ;
            DStats.meanCi(k,m) = (Ci(2) -Ci(1))/2 ;
        end
  %%%
%         markerAbmask_raw(:,:,2) = uint8(abMask*255) ;
%         abPairLbl1 = bwlabel(pairMask1) ;
%         abPairLbl2 = bwlabel(pairMask2) ; 
        IabSyn1 = cat(1,abStats1.MeanIntensity) ;
        IabSyn2 = cat(1,abStats2.MeanIntensity) ;
        IabSyn1 = IabSyn1(abPairInd1) ;
        IabSyn2 = IabSyn2(abPairInd2) ;
        [R,P,RL,RU] = corrcoef(IabSyn1, IabSyn2) ;
        pccIabAbSyn(k,m) = R(1,2) ;
        pccIabAbSynPval(k,m) = P(1,2) ;
        pccIabAbSynLc(k,m) = RL(1,2) ;
        pccIabAbSynUc(k,m) = RU(1,2) ;            

        abPairMask1 = ismember(abLbl1, abPairInd1) ; % remove puncta that are not paired with marker
        abPairMask2 = ismember(abLbl2, abPairInd2) ; % remove puncta that are not paired with marker
        markerAbmask(:,:,2) = uint8(abPairMask1*255) ;
        markerAbmask(:,:,3) = uint8(abPairMask2*255) ;        

%         figure(530)
%         imshow(markerAbmask,[], 'InitialMagnification', 'fit');    
%         hold on        
%         tlt_name = ['{\color{red}', ch_name{marker_ind}  '}, ', '{\color{green}', ch_name{ab_ind1}, '},',...
%              '{\color{blue}', ch_name{ab_ind2}, '}'] ;  
%         ht = title(tlt_name,'FontSize',30) ;   
%         print(530 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ; save([work_path, 'startup.mat'], 'fig_num');
%         if isempty(abCen1)||isempty(abCen2)
%         else
%             for j =1:size(abCenPair1,1)
%             %         plot(abCenPair(j,1), abCenPair(j,2), 'g+')
%                 line([abCenPair1(j,1) abCenPair2(j,1)],[abCenPair1(j,2) abCenPair2(j,2)],'Marker','+','LineStyle','-', 'color','w')
%             end 
%         end
%          
%         hold off                     
        
    end
end
%%
pccIabAbSynCi = (pccIabAbSynUc - pccIabAbSynLc)/2 ;
pccIabAbSynCell = mat2cell(pccIabAbSyn, ones(1,size(pccIabAbSyn,1)),  ones(1,size(pccIabAbSyn,2))) ;
pccIabAbSynCiCell = mat2cell(pccIabAbSynCi, ones(1,size(pccIabAbSynCi,1)),  ones(1,size(pccIabAbSynCi,2))) ;
pccIabAbSynCell = cellfun(@(x)num2str(x,'%0.2f'),pccIabAbSynCell,'UniformOutput', false);
pccIabAbSynCiCell = cellfun(@(x)num2str(x,'%0.2f'),pccIabAbSynCiCell,'UniformOutput', false);
pccIabAbSynCiCell = cellfun(@(c1,c2)sprintf([c1, '\n¡À', c2]), pccIabAbSynCell, pccIabAbSynCiCell,'UniformOutput', false);
% Deleted the upper triangle of the matrix
pccIabAbSynTri = pccIabAbSyn ;
nanL = tril(NaN(numel(ab_ind), numel(ab_ind)), -1)' ;
pccIabAbSynTri = pccIabAbSynTri + nanL ;
for j = 1:numel(ab_ind)
    for k = j+1:numel(ab_ind)
        pccIabAbSynCiCell{j,k} = ' ' ;
    end
end
%
pccIabAbSynCiTri = pccIabAbSynCi ;
% nanL = tril(NaN(numel(ab_ind_all), numel(ab_ind_all)), -1)' ;
pccIabAbSynCiTri = pccIabAbSynCiTri + nanL ;

%%
cmap = french(512, 3) ;
cmap = cmap(end:-1:1,:) ;
figure(49)
% heatmap(pccIabAbSyn, ch_name(ab_ind_all), ch_name(ab_ind_all),pccIabAbSynCiCell,'TickAngle', 45,'TickFontsize',15,'Fontsize', 15) ;% 'RowLabelsRotate', 45, 'ColumnLabelsRotate', 45) ;
% heatmap(pccIabAbSynTri, ch_name(ab_ind_all), ch_name(ab_ind_all),pccIabAbSynCiCell,...
heatmap(pccIabAbSynTri, ch_name, ch_name,pccIabAbSynCiCell,...
'TickAngle', 45,'TickFontsize',12,'Fontsize', 15,'Colormap', cmap,...
'NaNColor', [1 1 1], 'Textcolor', [0 0 0],'MinColorValue', -1, 'MaxColorValue',1, 'TickTexInterpreter',1) ;% 'RowLabelsRotate', 45, 'ColumnLabelsRotate', 45) ;
caxis([-1 1]) ;
colormap(cmap)
axis image
set(gca, 'box', 'off')
colorbar
title('Correlation of synaptic intensity', 'Fontsize', 20)
print(49 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ; save([work_path, 'startup.mat'], 'fig_num'); 
%% average asymetric distances
DmeanSym = (DStats.mean' + DStats.mean)/2 ;
DstdSym = (DStats.std' + DStats.std)/2 ;
DmeanCiSym = (DStats.meanCi' + DStats.meanCi)/2 ;
%%%
DmeanSymCell = mat2cell(DmeanSym, ones(1,size(DmeanSym,1)),  ones(1,size(DmeanSym,2))) ;
DmeanCiSymCell = mat2cell(DmeanCiSym, ones(1,size(DmeanCiSym,1)),  ones(1,size(DmeanCiSym,2))) ;
DstdSymCell = mat2cell(DstdSym, ones(1,size(DstdSym,1)),  ones(1,size(DstdSym,2))) ;
DmeanSymCell = cellfun(@(x)num2str(x,'%0.0f'),DmeanSymCell,'UniformOutput', false);
DmeanCiSymCell = cellfun(@(x)num2str(x,'%0.0f'),DmeanCiSymCell,'UniformOutput', false);
DstdSymCell = cellfun(@(x)num2str(x,'%0.0f'),DstdSymCell,'UniformOutput', false);
% DmeanSymCell = cellfun(@(c1,c2)sprintf([c1, '\n¡À', c2]), DmeanSymCell, DstdSymCell,'UniformOutput', false);
DmeanSymCell = cellfun(@(c1,c2)sprintf([c1, '\n¡À', c2]), DmeanSymCell, DmeanCiSymCell,'UniformOutput', false);
%% raw asymetric distances
DmeanCell = mat2cell(DStats.mean, ones(1,size(DStats.mean,1)),  ones(1,size(DStats.mean,2))) ;
DmeanCiCell = mat2cell(DStats.meanCi, ones(1,size(DStats.meanCi,1)),  ones(1,size(DStats.meanCi,2))) ;

DmeanCell = cellfun(@(x)num2str(x,'%0.0f'),DmeanCell,'UniformOutput', false);
DmeanCiCell = cellfun(@(x)num2str(x,'%0.0f'),DmeanCiCell,'UniformOutput', false);
DmeanCell = cellfun(@(c1,c2)sprintf([c1, '\n¡À', c2]), DmeanCell, DmeanCiCell,'UniformOutput', false);

%% Deleted the upper triangle of the matrix

nanL = tril(NaN(numel(ab_ind_all), numel(ab_ind_all)), -1)' ;
DmeanSymTri = DmeanSym + nanL ;
for j = 1:numel(ab_ind_all)
    for k = j+1:numel(ab_ind_all)
        DmeanSymCell{j,k} = ' ' ;
    end
end
%%
figure(57)
heatmap(abAbPairFrc, ch_name(ab_ind_all), ch_name(ab_ind_all),'%0.2f',...
     'TickTexInterpreter',1, 'Colormap', cmapr,'MinColorValue', 0, 'MaxColorValue',1,...
     'TickAngle', 45,'TickFontsize',15,'Fontsize', 15) ;% 'RowLabelsRotate', 45, 'ColumnLabelsRotate', 45) ;
% colormap summer
axis image
colorbar
title('Fraction of colocalized puncta', 'Fontsize', 20)
print(57 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ; 
print(57 ,'-depsc','-r300', [fig_path,num2str(fig_num, '%03d'),'.eps']) ; fig_num = fig_num +1 ; save([work_path, 'startup.mat'], 'fig_num'); 
%%
figure(58)
heatmap(DmeanSymTri , ch_name(ab_ind_all), ch_name(ab_ind_all), DmeanSymCell,...
'TickAngle', 45,'TickFontsize',15,'Fontsize', 15, 'TickTexInterpreter',1, 'Colormap', cmapr,...
'NaNColor', [1 1 1], 'Textcolor', [0 0 0]) ;
axis image
colorbar
title('Distance between colocalized puncta (\mum)', 'Fontsize', 20)
print(58 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ; save([work_path, 'startup.mat'], 'fig_num');

figure(59)
heatmap(DStats.mean, ch_name(ab_ind_all), ch_name(ab_ind_all), DmeanCell,...
'TickAngle', 45,'TickFontsize',15,'Fontsize', 15, 'TickTexInterpreter',1, 'Colormap', cmapr,...
'NaNColor', [1 1 1], 'Textcolor', [0 0 0]) ;
axis image
colorbar
title('Distance between colocalized puncta (\mum)', 'Fontsize', 20)
print(59 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ; save([work_path, 'startup.mat'], 'fig_num');
%% clustering analysis
% DataMat = [IabInMarker(:,3:end)'; AabWithMarker(:,3:end)']  ;
DataMat = [IabWithMarker(:,3:end)'; AabWithMarker(:,3:end)']  ;
% DataMat = zscore(DataMat,0,2) ;
% DataMat = DataMat./repmat(mean(DataMat,2), [1 size(DataMat,2)]) ;
for j = 1:size(DataMat,1)
    DataMatS = DataMat(j,:) ;
    DataMatS = DataMatS./std(DataMatS(DataMatS>0)) ;
    DataMatS(DataMatS>0) = DataMatS(DataMatS>0) - min(DataMatS(DataMatS>0)) ;
    DataMat(j,:) =  DataMatS ;
end
%%
rowlabI = cellfun(@(c2)(['I(',c2, ')']), ch_name(3:end),'UniformOutput', false);
rowlabA = cellfun(@(c2)(['A(',c2, ')']), ch_name(3:end),'UniformOutput', false);
rowlab = [rowlabI, rowlabA] ;

%% find optimal number of clusters using silhouette distance
clusfunc = @(X,K)(cluster(linkage(X,'ward'),'maxclust',K));
% E = evalclusters(DataMat','linkage','silhouette','klist',[1:size(DataMat,2)]) ;
% E = evalclusters(DataMat,clusfunc,'silhouette','klist',[1:size(DataMat,1)]) ;
EvRow = evalclusters(DataMat,clusfunc,'silhouette','klist',[1:12]) ;
EvCln = evalclusters(DataMat',clusfunc,'silhouette','klist',[1:30]) ;
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
print(60 ,'-depsc','-r300', [fig_path,num2str(fig_num, '%03d'),'.eps']) ; fig_num = fig_num +1 ; save([work_path, 'startup.mat'], 'fig_num');
print(61 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; 
print(61 ,'-depsc','-r300', [fig_path,num2str(fig_num, '%03d'),'.eps']) ; fig_num = fig_num +1 ; save([work_path, 'startup.mat'], 'fig_num');
%%
treeRow = linkage(DataMat,'ward') ;
treeCln = linkage(DataMat','ward') ;
DenThr(1) = treeRow(end-EvRow.OptimalK+2,3)-eps;
% DenThr(2) = treeCln(end-EvCln.OptimalK+2,3)-eps;
% DenThr(2) = treeCln(end-2+2,3)-eps;
DenThr(2) = 45;
%%
cgo = clustergram(DataMat,'Linkage','ward') ;
%%%
cmap = french(512, 3) ;
cmap = cmap(end:-1:1,:) ;
set(cgo,'RowLabels',rowlab, 'Colormap', cmapr, 'Dendrogram', DenThr,'Symmetric', 'false','DisplayRange', 6)
% set(cgo,'RowLabels',rowlab, 'Colormap', cmap, 'DisplayRange', 3)
% figure(60)
h_cgo = plot(cgo) ;
%%
print(gcf ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; 
print(gcf ,'-dpdf','-r300', [fig_path,num2str(fig_num, '%03d'),'.pdf']) ; fig_num = fig_num +1 ; save([work_path, 'startup.mat'], 'fig_num');
print(gcf ,'-depsc','-r300', [fig_path,num2str(fig_num, '%03d'),'.eps']) ;
% addYLabel(cgo , rowlab, 'FontSize', 8)


end


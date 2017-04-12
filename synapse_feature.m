function  synapse_feature(params)
% Each synaptic marker punctum defines one synapse 
% extract synapse features using the segmented masks, at the moment 3
% features are extracted for each channel:
% (1) mean intenisties within each segmented marker punctum for each channel (IabInMarker) 
% (2) mean intenisties within each segmented target punctum that are colocalized with a marker punctum for each channel (IabwithMarker)
% (3) area of each segmented target punctum that are colocalized with a marker punctum for each channel (AabwithMarker)
    global fig_num fig_path work_path;
    imgCondiFolders = readDirSubfolders(params.outputImgsPath,'all');
    for j = 1:numel(imgCondiFolders)        
        if params.tophat == 1            
            load(fullfile(params.outputImgsPath, imgCondiFolders{j}, 'MultiplexImageAlignedTopHat.mat'),...
                'imgsProjTopHat','imgRoundNames','params') ;   
            imgsProjAligned = imgsProjTopHat ;
            clear imgsProjTopHat          
        else
            load(fullfile(params.outputImgsPath,  imgCondiFolders{j},'MultiplexImageDataAligned.mat'),...
                'imgsProjAligned','imgRoundNames','params');                      
        end                            

        close all
        maskfile = fullfile(params.outputImgsPath,imgCondiFolders{j}, 'mask.mat') ;
        load(maskfile);
        ch_name = imgRoundNamesNoWO ; %the channel name
        marker_ind_bi = cellfun(@(x) ~isempty(regexpi(x,'synapsin', 'match')), imgRoundNamesNoWO,'uniformoutput',0) ; % find marker channel (synapsin-1)
        marker_ind_bi = cat(1,marker_ind_bi{:}) ;
        marker_ind = find(marker_ind_bi, 1, 'first') ; % pick the first marker channel if there are more than 1
        ab_ind_all = [1:numel(imgRoundNamesNoWO)] ; % indices of other antibodies
        ab_ind_all(marker_ind) = [] ;
        uniqueFields = [1:numel(imgsProjAligned)] ;
        pccIabInMarkerAll = [] ; % Pearson correlation matrix of features 
        pccIabWithMarkerAll = [] ; 
        pccAabWithMarkerAll = [] ;
    % for f=2:2
        for f = 1:numel(uniqueFields)
            disp(['process ', params.Condi,  ' field ', num2str(f), '...'])          
            ItargetStack = cat(3,imgsProjAligned{f}{3,targets});
            ItargetStack = cat(3,ItargetStack,imgsProjAligned{f}{4,1}); % Add VGlutT channel to target stack   
            ItargetStack = cat(3,ItargetStack,imgsProjAligned{f}{1,1}); % Add MAP2 channel to target stack            
            Imarker = ItargetStack(:,:,marker_ind) ;
            pixsize = 0.18733333 ;   % pixel size in um                           
            %%
            markerMask = synMask{f}(:,:,marker_ind) ;
            abMaskAll =  synMask{f}(:,:,ab_ind_all) ;
            params.minObjSize = 12 ;
            markerMask = bwareaopen(markerMask, params.minObjSize); % remove small objects that are mostly likely noise
            markerMask = imclearborder(markerMask) ;            
            markerStats = regionprops(markerMask,Imarker,'Area', 'WeightedCentroid','MeanIntensity');            
            markerCen = cat(1, markerStats.WeightedCentroid); % get marker centroids
            markerArea = cat(1, markerStats.Area);% get marker centroids
            nan_ind = any(isnan(markerCen),2) ;
            markerCen(nan_ind, :) = [] ;
            markerStats(nan_ind) = [] ;
            markerArea(nan_ind) = [] ; 
            markerAbmask  = zeros([size(markerMask), 3], 'uint8') ;% preallocate overlaid marker and ab masks after distance filtering 
            markerAbmask_raw  = zeros([size(markerMask), 3], 'uint8') ;
            markerAbmask(:,:,1) = uint8(markerMask*255) ;
            markerAbmask_raw(:,:,1) = uint8(markerMask*255) ;
            abStatsAll = cell(1, numel(ab_ind_all)) ;

            %% show marker mask and centers 
            n_obj = size(markerCen,1) ;
            if params.plot_opt == 1
                figure(41)
                imshow(markerMask,[], 'InitialMagnification', 'fit');
                tlt_name = ch_name{marker_ind} ;  
                ht = title(tlt_name,'FontSize',20) ;
                hold on

                for j =1:n_obj
                    plot(markerCen(j,1), markerCen(j,2), 'r+')
                    text(markerCen(j,1)+6,markerCen(j,2)+6,num2str(j),'FontSize',10,'color','r')   
                end
                hold off
                print(41 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ; save([work_path, 'startup.mat'], 'fig_num');
            end
            %% measure intensity within synapsin mask for all channels
            IabInMarker = zeros(n_obj, numel(imgRoundNamesNoWO)) ;
            IabInMarkerTtl = IabInMarker ;                 
            for k = 1:numel(imgRoundNamesNoWO)                            
                Iab = ItargetStack(:,:,k);
                AbInMarkerStats = regionprops(markerMask,Iab, 'MeanIntensity');
                IabInMarker(:,k) =  cat(1, AbInMarkerStats.MeanIntensity) ;        
                IabInMarkerTtl(:,k) = IabInMarker(:,k).*markerArea ;                         
            end
            IabInMarkerCell{f}= IabInMarker ;     

            %% calculate punctae distantace to marker, keep  punctae that are colocalized with the markers only
            cut_off_v = [0:50:1000] ;
            markerPairFrcV = zeros(1,numel(cut_off_v))  ;
            abPairFrcV = zeros(1,numel(cut_off_v))  ;
            markerPairIndAll = [] ;
            IabWithMarker = zeros(size(markerStats,1), numel(imgRoundNamesNoWO)) ;
            AabWithMarker = zeros(size(markerStats,1), numel(imgRoundNamesNoWO)) ;
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
                abStats = regionprops(abMask, Iab, 'Area', 'WeightedCentroid', 'MeanIntensity');      
                abCen = cat(1, abStats.WeightedCentroid);               
                [IDX,D]=knnsearch(abCen,markerCen,'k',5); % find the first 5 closest target punctae to each marker punctum                    
                n_bin = 40 ;
                D = D*pixsize*1000; %convert distance to nm
                %%%
                figure(50) % show histograms of distances to neareast puncta
                h1 = histogram(D(:,1),n_bin) ;            
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
            
                markerPairIndBi = D(:,1)< 1000; % set cut-off distance 1 um to define colocolization, based on the above distance histograms
                markerUnpairIndBi = ~markerPairIndBi ;
            %     markerPairInd = find(markerPairIndBi) ;
                markerPairInd = IDX(:,1) ; % For each marker punctum find the index of paired ab punctum. Same size as markerCenUnpaired. Unpaired maker puncta are assgined "0"  
                markerPairInd(markerUnpairIndBi) = 0  ;             
                markerCenPair = markerCen(markerPairIndBi,:) ;            
            %%%
                abPairInd = unique(IDX(markerPairIndBi,1)) ;%unique ab paired indices
            %     abStats.pairInd = abPairInd ;
                abPairIndNonUni = IDX(markerPairIndBi,1) ; % non-unique ab paired indices
                abCenPair = abCen(abPairIndNonUni,:) ; 
                abPairFrc(k) = numel(abPairInd)/size(abCen,1) ;
                markerPairFrc(k) = sum(markerPairIndBi)/size(markerCen,1) ;
                %%%
           
                IabSyn = markerPairInd ; % preallocate intensity of the colocalized target punctae
                IabSynAll = cat(1,abStats.MeanIntensity) ;
            %     IabSynAll = IabSynAll.* cat(1,abStats.Area) ;
                IabSyn(IabSyn>0) = IabSynAll(abPairIndNonUni) ; % assign zero intensity to unpaired puncta
                IabWithMarker(:,ab_ind) = IabSyn ;

                AabSyn = markerPairInd ;
                AabSynAll = cat(1,abStats.Area) ;    
                AabSyn(AabSyn>0) = AabSynAll(abPairIndNonUni) ; % assign zero area to unpaired puncta
                AabWithMarker(:,ab_ind) = AabSyn ;            %%%    
            
                [DStats.mean(k), DStats.std(k), Ci,~] = normfit(D(markerPairIndBi,1),0.05) ; %
                DStats.meanCi(k) = (Ci(2) -Ci(1))/2 ;            
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
    %             print(52 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ; save([work_path, 'startup.mat'], 'fig_num');

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
            % replace intensities of cytoskeletal targets that cannot be segmented properly with intensity measured within marker masks 
            nonsyn_ind_bi = cellfun(@(x) ~isempty(regexpi(x,'tuj-1|MAP2', 'match')), imgRoundNamesNoWO,'uniformoutput',0) ; 
            nonsyn_ind_bi = cat(1,nonsyn_ind_bi{:}) ;
            IabWithMarker(:,nonsyn_ind_bi) = IabInMarker(:,nonsyn_ind_bi) ;
            [pccIabWithMarker, pccIabWithMarkerPval, pccIabWithMarkerLc, pccIabWithMarkerUc]  = corrcoef(IabWithMarker) ; % calculate correlation matrices
            [pccAabWithMarker, pccAabWithMarkerPval, pccAabWithMarkerLc, pccAabWithMarkerUc]  = corrcoef(AabWithMarker) ;
            pccIabWithMarkerAll = cat(3, pccIabWithMarkerAll, pccIabWithMarker) ;
            pccAabWithMarkerAll = cat(3, pccAabWithMarkerAll, pccAabWithMarker) ;
            %% Heat Map of PCC of integrated intensity + 95% confidence interval 
            [pccIabInMarker, pccIabInMarkerPval, pccIabInMarkerLc, pccIabInMarkerUc]  = corrcoef(IabInMarker) ;            
            pccIabInMarkerAll = cat(3, pccIabInMarkerAll, pccIabInMarker) ;
        end
    %%
        save(maskfile, 'IabInMarkerCell','IabWithMarkerCell','AabWithMarkerCell',...
            'pccIabInMarkerAll', 'pccIabWithMarkerAll', 'pccAabWithMarkerAll', '-append');              
        %%
%         figure(103)
%         distributionPlot(IabInMarkerAll,'xyOri','flipped','histOri','right','showMM',6,'color', [0.7 0.7 0.7],'histOpt',2)
%         set(gca,'FontSize', 20)
%         set(gca,'LineWidth', 2)
%         set(gca,'TickLength'  , [.02 .02])
%         set(gca,'yticklabel'  ,ch_name)
%         xl = xlabel('Intensity (AU)') ; 
%         % % ylabel('log(Count)')
%         ylabel('Count') 
%         format_fig2(5)
%         hold off
%         %%
%         figure(104)
%         distributionPlot(AabWithMarkerAll,'xyOri','flipped','histOri','right','showMM',6,'color', [0.7 0.7 0.7],'histOpt',2)
%         set(gca,'FontSize', 20)
%         set(gca,'LineWidth', 2)
%         set(gca,'TickLength'  , [.02 .02])
%         set(gca,'yticklabel'  ,ch_name)
%         xl = xlabel('Intensity (AU)') ; 
%         % % ylabel('log(Count)')
%         ylabel('Count') 
%         format_fig2(5)
%         hold off
        
        %%
%             print(103 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; 
%             print(103 ,'-dpdf','-r300', [fig_path,num2str(fig_num, '%03d'),'.pdf']) ; fig_num = fig_num +1 ; save('.\startup.mat', 'fig_num'); 
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


    end
end




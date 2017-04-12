function multi_condition_confocal_plot(params, CondiInd, PlotRoundInd, ChanInd, sat_frc)
global fig_num fig_path;
%-----------------------------------------------------------------------------------------------------------------------
% opt_corp: export the full field of view plus manually cropped regions (1) or the full feild of views only (0)

% Path of the parent folder for current experiment and analysis
% params.parentFolderForAnalysis = 'E:\Dropbox\Dropbox (MIT)\Neurons\Nucleic Acid-based Multiplexing Analysis\051316' ;
% 
% 
% %-----------------------------------------------------------------------------------------------------------------------
% %Set the name of the folder that contains the folders of images of all fields per target. 
% params.inputImgsPath = 'C:\Users\Simon\Desktop\LNA Phenix Images 051316' ;
% 
% %-----------------------------------------------------------------------------------------------------------------------
% %Set the type of projection to do on the z-stack for each target and field of view, as 'max' or 'average'
% params.projType = 'average' ;
% 
% %-----------------------------------------------------------------------------------------------------------------------
% %Set the name of the output folder that will contains the pre-processed images. Within this folder will be folders per target, and within per-target folders,
% %folders for different fields of view
% params.outputImgsPath = 'E:\Dropbox\Dropbox (MIT)\Neurons\Nucleic Acid-based Multiplexing Analysis\051316\post_processing' ;
% 
% %---------------------------------------------------------------------------------------------------
% %Set the channel names of the nucleus, MAP2, and LNA target, in order of channel number.
% params.channels.cellimgs = {'MAP2','DAPI','LNA'} ;
% addpath(genpath(fullfile(params.parentFolderForAnalysis,'Scripts'))) ;
    
    load(fullfile(params.outputImgsPath, params.Condi,'MultiplexImageData.mat'),'imgsProjAligned',...
        'imgRoundNames');                
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %PLOTTING

    close all;
    params.pixsize = 0.187 ; % pixel size in um
   [~,plate.Condi] = xlsread(fullfile(params.parentFolderForAnalysis,params.plateMapFile),'Condi');                        
    uniqueCondibi = cellfun(@(x) strcmp(x,'N/A'),plate.Condi,'uniformoutput',0) ;
    uniqueCondibi = ~cell2mat(uniqueCondibi) ;
%             uniqueCondiInd = find(uniqueCondibi) ;
    uniqueCondi = plate.Condi(uniqueCondibi) ;
    uniqueFields = [1:numel(imgsProjAligned)] ;         
    uniqueCondi = uniqueCondi(CondiInd) ;    
    [imgTitles, ~] = ParseImgName(imgRoundNames) ;
    
    sizeI = size(imgsProjAligned{1}{1,1}); % get the size of the first image
    % automatically determine the contrast
    [imgs_all, cmin_condi, cmax_condi] = contrast_est(params, uniqueCondi, PlotRoundInd, ChanInd, sat_frc,sizeI) ;
    [uniqueCondi, ~] = ParseAbName(uniqueCondi) ;
%     for f = 1:numel(uniqueFields)
    for f = 3:3
%         tempImgName = fullfile(params.outputImgsPath,['field_',num2str(uniqueFields(f)),'_Processed.pdf']);
        for r = 1:numel(PlotRoundInd)
            tempImgName = fullfile(params.outputImgsPath,['field_',num2str(uniqueFields(f)),'_',...
                imgRoundNames{r},'_Processed.png']);
            figure(100*f+r)
            pos = get(0,'screensize') ;
            pos(4) = pos(4)*0.8 ;
            set(gcf,'position',pos);
            set(gcf,'color','k')        
            for c = 1:numel(uniqueCondi)
                imgs = imgs_all{c} ;                
                IStack = imgs{f}(1:3,PlotRoundInd(r));     
                IStack = permute(IStack, [2 3 1]) ;
                IStack = cell2mat(IStack) ;
                 if isempty(IStack)
                    IStack = zeros(sizeI(1), sizeI(2), 3) ;                
                 end
                chan_invis = setdiff([1:3],ChanInd) ; % make unselected channel invisible               
                IStack(:,:,chan_invis) = 0 ; % hide MAP2                
                IStack = IStack(:,:,[3,1,2]) ;                   
                %%
                sat_frc = [cmin_condi(r,:)', cmax_condi(r,:)'] ;
%                 sat_frc = [1.0e-03 *[0.7172;    0.3510;    0.3662],[0.0595;    0.1221;    0.1217]];
                cmap = hsv(size(IStack, 3)) ;
                imr2 = imstack2RGB(IStack,sat_frc, {'limit', 'limit', 'limit'}, cmap);               
                figure(100*f+r)
                subplot_tight(2,ceil(numel(uniqueCondi)/2),c,[0.05 0.01])
                imshow(imr2)
                hsb = scalebar('scalelength_um', 20, 'pixelsize', params.pixsize, 'location', 'southeast') ;   
                title(uniqueCondi(c),'color','w', 'fontsize', 20)
%                 format_fig2(9)
                       
            end
            ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],...
            'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
            text(0.5, 1,['\bf', imgTitles{r}],'HorizontalAlignment'... 
            ,'center','VerticalAlignment', 'top', 'color','w', 'fontsize', 22)                

%             suptitle(imgTitles{r})
            %Export
            figure(100*f+r)
            export_fig(tempImgName,'-painters','-append') 
            close gcf;
        end
    end        
end

function [imgs_all, cmin_condi, cmax_condi] = contrast_est(params, uniqueCondi, PlotRoundInd, ChanInd, sat_frc, sizeI)
% estimate contrast range for displaying images
    cmin_arr = zeros(numel(uniqueCondi), numel(PlotRoundInd),numel(ChanInd)) ;
    cmax_arr = zeros(numel(uniqueCondi), numel(PlotRoundInd),numel(ChanInd)) ;
    imgs_all = cell(numel(uniqueCondi),1) ;
    f = 1; % use the first field to determine the contrast
    for c = 1:numel(uniqueCondi)
        if params.tophat == 1
            load(fullfile(params.outputImgsPath, uniqueCondi{c},'MultiplexImageTopHat.mat'),'imgsProjTopHat');
            imgs = imgsProjTopHat ;
            clear imgsProjTopHat            
        else
            load(fullfile(params.outputImgsPath, uniqueCondi{c},'MultiplexImageData.mat'),'imgsProjAligned');
            imgs = imgsProjAligned ;
            clear imgsProjAligned         
        end
        imgs_all{c} = imgs ;
        for r = 1:numel(PlotRoundInd)                                                 
            IStack = imgs{f}(1:3,PlotRoundInd(r));       
            IStack = permute(IStack, [2 3 1]) ;
            IStack = cell2mat(IStack) ;
            IStack = IStack(:,:,[3,1,2]) ;
            if isempty(IStack)
                IStack = zeros(sizeI(1), sizeI(2), 3) ;
                cmin_arr(c,r,:) = NaN;
                cmax_arr(c,r,:) = NaN;
            else
                for z = 1:size(IStack,3)                    
                    Low_High = stretchlim(IStack(:,:,z), sat_frc(z,:)) ;    
                    cmin_arr(c,r,z) = Low_High(1);
                    cmax_arr(c,r,z) = Low_High(2);
                end
            end
        end
    end
    cmin_condi = squeeze(min(cmin_arr, [], 1)) ;
    cmax_condi = squeeze(max(cmax_arr, [], 1)) ;    
end

function [AbNames, targets]= ParseAbName(AbNames)
% (1)skip washout rounds, basline, duplicates
% (2)replace seqence names with protein names
    seq_target_cell = {'actin-p2','p2';'Tuj-1-p3','p3';'cortactin-p4','p4';...
        'SHANK3-p6','p6';'ARPC2-p7','p7'; 'bassoon-p8','p8'; 'synapsin-1-p9(2)','p9-t2'; 'synapsin-1-p9','p9';...
        'Homer-1b/c-p10','p10';'NR2B-p12','p12';'PSD-95-p1','p1-'} ;
    targets = [] ;
    
    for j = 1:numel(AbNames)
        if isempty(regexpi(AbNames{j}, '(wo|baseline|x)', 'match'))
            targets = [targets, j] ;                    
        end                       
        for k =  1:size(seq_target_cell,1)            
           if ~isempty(regexpi(AbNames{j},seq_target_cell{k,2}, 'match'))
                AbNames{j} = seq_target_cell{k,1} ;
            end 
        end
    end
end

function [imgRoundNames, targets]= ParseImgName(imgRoundNames)
% (1)skip washout rounds, basline, duplicates
% (2)replace seqence names with protein names
    seq_target_cell = {'actin','p2';'Tuj-1','p3';'cortactin','p4';...
        'SHANK3','p6';'ARPC2','p7'; 'bassoon','p8'; 'synapsin-1(2)','p9-t2'; 'synapsin-1','p9';...
        'Homer-1b/c','p10';'NR2B','p12';'PSD-95','p1-'} ;
    targets = [] ;
    
    for j = 1:numel(imgRoundNames)
        if isempty(regexpi(imgRoundNames{j}, '(wo|baseline|x)', 'match'))
            targets = [targets, j] ;
            for k =  1:size(seq_target_cell,1)            
               if ~isempty(regexpi(imgRoundNames{j},seq_target_cell{k,2}, 'match'))
                    imgRoundNames{j} = [seq_target_cell{k,2}, ' imaging probe'] ;
               end 
            end                    
        end                               
    end
end     


            
            
            
            
 

        
        
        
        
        
        
        
        
        
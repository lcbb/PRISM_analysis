function multiplex_confocal_plot(params,sat_frc, opt_crop)
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
imgCondiFolders = readDirSubfolders(params.outputImgsPath,'all');    
    for j = 1:numel(imgCondiFolders)
        params.tophat = 1  ; % redefine tophat option as it gets overwritten in each round
        if params.tophat == 1            
            load(fullfile(params.outputImgsPath, imgCondiFolders{j}, 'MultiplexImageAlignedTopHat.mat'),...
                'imgsProjTopHat','imgRoundNames','params') ;   
            imgsProjAligned = imgsProjTopHat ;
            clear imgsProjTopHat          
        else
            load(fullfile(params.outputImgsPath,  imgCondiFolders{j},'MultiplexImageDataAligned.mat'),...
                'imgsProjAligned','imgRoundNames','params');                      
        end                            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %PLOTTING
        params.pixsize = 0.187 ; % pixel size in um
        uniqueFields = [1:numel(imgsProjAligned)] ;
        uniqueChannels = size(imgsProjAligned{1},1) ;       
        close all;
        % parse image round names
        [imgRoundNames, targets]= ParseImgName(imgRoundNames) ;
        if uniqueChannels>3
            imgRoundNames = [imgRoundNames(targets(1:end-1)); {'VGluT'; 'MAP2';'DAPI'}; imgRoundNames(targets(end))] ;
        else
            imgRoundNames = [imgRoundNames(targets(1:end-1)); { 'MAP2';'DAPI'};imgRoundNames(targets(end))] ;
        end
%         targets  = [1:numel(imgRoundNames)] ;        
        
%         for f = 1:numel(uniqueFields)
        for f = 3:3            
            sizeI = size(imgsProjAligned{f}{1,1});            
            ItargetStack = cat(3,imgsProjAligned{f}{3,targets});
            if uniqueChannels>3
                ItargetStack = cat(3,ItargetStack(:,:,1:end-1),imgsProjAligned{f}{[4,1,2],1}, ItargetStack(:,:,end));            
            else
                ItargetStack = cat(3,ItargetStack(:,:,1:end-1),imgsProjAligned{f}{[1,2],1}, ItargetStack(:,:,end));       
            end
            cmap = hsv(size(ItargetStack, 3)-1) ;
            cmap = cmap([1:9,13, 11:12, 10, 5],:) ; % change the color order
            imr2 = imstack2RGB(mat2gray(ItargetStack),sat_frc,{'fraction'}, cmap);            
            imr2 = mat2gray(imr2) ;            
            
            figure(1);
%             set(gcf,'position',get(0,'screensize'));
            
            subplot_tight(3,5,1,[0.05 0.01])
            imshow(imr2)
            hsb = scalebar('scalelength_um', 20, 'pixelsize', params.pixsize, 'location', 'southeast') ;   
            title('Composite','color','w', 'fontsize', 20)
            
             for i = 1:size(ItargetStack,3)                
                I_temp = mat2gray(ItargetStack(:,:,i));
                Iplot = zeros(sizeI(1),sizeI(2),numel(targets));
                Iplot(:,:,i) = I_temp;
                imr = imstack2RGB(Iplot,sat_frc,{'fraction'}, cmap);                
                figure(1);
                subplot_tight(3,5,i+1,[0.05 0.01])
                imshow(imr)
                hsb = scalebar('scalelength_um', 20, 'pixelsize', params.pixsize, 'location', 'southeast') ;   
                title(imgRoundNames(i),'color','w', 'fontsize', 20)                
             end
            
            figure(1);            
            set(gcf,'position',[ 100,100,1500, 1000]); 
            set(gcf,'color','k')
            tightfig();
            
            if opt_crop == 1
                figure(400)    
                imshow(imr2, 'InitialMagnification', 'fit')           
                caxis auto                
                format_fig2(1)            
                set(gcf,'units','normalized','outerposition',[0 0 1 1])
                title('Please selet ROI to crop...','BackgroundColor','r','FontSize',10)%% crop image
                waitfor(gcf, 'selectiontype', 'alt')
                rect_sup = getrect;
    %                     rect_sup(3:4) = max(rect_sup(3:4))*[1,1] ; %make ROI square 
                rectangle('Position',rect_sup,'LineStyle','--','LineWidth',3, 'edgeColor','w');
                imr2_cropped = imcrop(imr2, rect_sup);                
    %                     ItargetStack_cropped = imsequence_crop(ItargetStack, rect_sup);

                figure(2);
                subplot_tight(3,5,1, [0.05 0.01])
                imshow(imr2_cropped)
                title('Composite','color','w', 'fontsize', 20)
                hsb = scalebar('scalelength_um', 5, 'pixelsize', params.pixsize, 'location', 'southeast') ;   
                for i = 1:size(ItargetStack,3)

                    I_temp = mat2gray(ItargetStack(:,:,i));
                    Iplot = zeros(sizeI(1),sizeI(2),numel(targets));
                    Iplot(:,:,i) = I_temp;
                    imr = imstack2RGB(Iplot,sat_frc,{'fraction'}, cmap);
                    imr_cropped = imcrop(imr, rect_sup);  
%                     imr_cropped = imsequence_crop(Iplot, rect_sup);
%                     imr_cropped = imstack2RGB(imr_cropped,sat_frc,{'fraction'}, cmap);
                    
                    figure(2);
                    subplot_tight(3,5,i+1, [0.05 0.01])
                    imshow(imr_cropped)
                    hsb = scalebar('scalelength_um', 5, 'pixelsize', params.pixsize, 'location', 'southeast') ;   
                    title(imgRoundNames(i),'color','w', 'fontsize', 20)
                end
                figure(2);
                tightfig();
                set(gcf,'position',get(0,'screensize'));
                set(gcf,'color','k')
            else
            end                                               
             %%
            %Export
            tempImgName = fullfile(params.outputImgsPath,params.Condi,['field_',num2str(uniqueFields(f)),'_Processed.pdf']);
            tempImgName2 = fullfile(params.outputImgsPath,params.Condi,['field_',num2str(uniqueFields(f)),'_Processed.png']);
            figHandles = findall(0,'Type','figure'); % get all figure handles
            for numFigs = 1:numel(figHandles)
                export_fig(tempImgName,'-painters','-append') 
                export_fig(tempImgName2,'-painters','-append') 
                close gcf;
            end
        end
    end
end

function [imgRoundNames, targets]= ParseImgName(imgRoundNames)
% (1)skip washout rounds, basline, duplicates
% (2)replace seqence names with protein names
%     seq_target_cell = {'actin','p2';'Tuj-1','p3';'cortactin','p4';...
%         'Shank3','p6';'ARPC2','p7'; 'bassoon','p8'; 'synapsin-1','p9';...
%         'Homer-1b/c','p10';'NR2B','p12';'PSD-95','p1-'} ;
    seq_target_cell = {'actin','p2';'Tuj-1','p3';'cortactin','p4';...
        'SHANK3','p6';'ARPC2','p7'; 'bassoon','p8'; 'synapsin-1(2)','p9-t2'; 'synapsin-1','p9';...
        'Homer-1b/c','p10';'NR2B','p12';'PSD-95','p1-'} ;
    targets = [] ;
    
    for j = 1:numel(imgRoundNames)
        if isempty(regexpi(imgRoundNames{j}, '(wo|baseline|all|x)', 'match'))
            targets = [targets, j] ;                    
        end                       
        for k =  1:size(seq_target_cell,1)            
           if ~isempty(regexpi(imgRoundNames{j},seq_target_cell{k,2}, 'match'))
                imgRoundNames{j} = seq_target_cell{k,1} ;
            end 
        end
    end
end     

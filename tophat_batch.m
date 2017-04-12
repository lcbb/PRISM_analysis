function tophat_batch(params)
% iterate through each folder and apply top-hat filter    
    imgCondiFolders = readDirSubfolders(params.outputImgsPath,'all');    
    for j = 1:numel(imgCondiFolders)
        load(fullfile(params.outputImgsPath, imgCondiFolders{j}, 'MultiplexImageData.mat'));
        uniqueFields = [1:numel(imgsProj)] ;
        uniqueChannels = [1:size(imgsProj{1},1)] ;
        params.Condi = imgCondiFolders{j} ;    
        params.disksize = 100;                   
        imgsProjTopHat =  cell(numel(uniqueChannels),numel(imgRoundNames)) ; % generate 3 masks for each round % 
        imgsProjTopHat = repmat({imgsProjTopHat},numel(uniqueFields),1);                               

        for f = 1:numel(uniqueFields)
            disp(['process ', params.Condi,  ' field ', num2str(f), '...'])
            for roundInd = 1:numel(imgRoundNames)
                sizeI = size(imgsProj{f}{1,roundInd});
                Idapi = imgsProj{f}{2,roundInd};
                IMAP2 = imgsProj{f}{1,roundInd};            
                ILNA = imgsProj{f}{3,roundInd};
                if numel(uniqueChannels)>3
                    IVGluT = imgsProj{f}{4,roundInd};
                end
                if isempty(Idapi)||isempty(IMAP2)...
                        ||~any(Idapi(:))||~any(IMAP2(:)) % if any channel is empty or all zero, leave the mask empty
                else
%                     figure(1)
%                     space_h = 0.01 ;   space_v = 0.01 ;
%                     set(gcf,'position',[100,100,1500, 1000]); 
                    for c = 1:numel(uniqueChannels)
                        I = imgsProj{f}{c,roundInd};
                        I_bg = illum_prof_c{c} ;                                                                                
                        I_norm = single(I)./single(I_bg)*100 ;                                            
                        I_tophat = imtophat(I,strel('disk',params.disksize));
                        I_tophat_norm = single(I_tophat)./single(I_bg)*100 ;
                        if any(any(I_norm>65535))||any(any(I_tophat_norm>65535))
                            warning('Image intenisty values exceed uint16 limit. Reduce the scaling factor')
                        end
                        I_norm = uint16(I_norm) ;
                        I_tophat_norm = uint16(I_tophat_norm) ;
                        imgsProj{f}{c,roundInd} = I_norm ;
                        imgsProjTopHat{f}{c,roundInd} = I_tophat_norm;                        
%                         subplot_tight(numel(uniqueChannels),5,5*(c-1)+1,[space_v space_h])
%                         imshow(imadjust(I),[], 'InitialMagnification', 'fit');
%                         subplot_tight(numel(uniqueChannels),5,5*(c-1)+2,[space_v space_h])
%                         imshow(I_bg,[], 'InitialMagnification', 'fit');
%                         subplot_tight(numel(uniqueChannels),5,5*(c-1)+3,[space_v space_h])
%                         imshow(imadjust(mat2gray(I_norm)),[], 'InitialMagnification', 'fit');
%                         subplot_tight(numel(uniqueChannels),5,5*(c-1)+4,[space_v space_h])
%                         imshow(imadjust(I_tophat),[], 'InitialMagnification', 'fit');
%                         subplot_tight(numel(uniqueChannels),5,5*(c-1)+5,[space_v space_h])
%                         imshow(imadjust(mat2gray(I_tophat_norm)),[], 'InitialMagnification', 'fit');
%                         drawnow
                    end
                end       
            end
        end
        params.tophat = 1 ;    
        outputfile = fullfile(params.outputImgsPath,imgCondiFolders{j}, 'MultiplexImageTopHat.mat') ;
        save(outputfile, 'imgsProjTopHat', 'params', 'imgRoundNames','illum_prof_c')
        params.tophat = 0 ; 
        outputfile = fullfile(params.outputImgsPath,imgCondiFolders{j}, 'MultiplexImageData.mat') ;
        save(outputfile, 'imgsProj', 'params', 'imgRoundNames','illum_prof_c', '-append')
    end
    %%
    %  figure(5)
    % for  roundInd =2:4
    % subplot_tight(2,3,roundInd-1,[0.05 0.01])
    % imshow(imgsProj{f}{3,roundInd},[], 'InitialMagnification', 'fit');
    % end
    %   
end

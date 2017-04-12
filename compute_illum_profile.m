function compute_illum_profile(params)
    imgCondiFolders = readDirSubfolders(params.outputImgsPath,'all');
    load(fullfile(params.outputImgsPath, imgCondiFolders{1}, 'MultiplexImageData.mat'));
    uniqueFields = [1:numel(imgsProj)] ;
    uniqueChannels = [1:size(imgsProj{1},1)] ;    
    illum_prof_confc = cell(numel(imgCondiFolders),numel(uniqueFields),numel(uniqueChannels)) ; % illumination profile for each channel
    illum_prof_conc = cell(numel(imgCondiFolders),numel(uniqueChannels)) ; % illumination profile for each field and channel
    illum_prof_c = cell(numel(uniqueChannels),1) ; % illumination profile for each channel
    for j = 1:numel(imgCondiFolders) ;              
        load(fullfile(params.outputImgsPath, imgCondiFolders{j}, 'MultiplexImageData.mat'));                             
        
%         figure(4)
%         space_h = 0.01 ;
%         space_v = 0.01 ;
%         ind = 1 ;
        for f = 1:numel(uniqueFields)
            for ch = 1:numel(uniqueChannels)
                 for roundFolder = 1:numel(imgRoundNames)
                    imgsProj{f}{ch,roundFolder} = imopen(imgsProj{f}{ch,roundFolder},strel('disk',params.disksize)) ;
                 end
                illum_prof_confc{j,f,ch} = uint16(mean(cat(3,imgsProj{f}{ch,:}),3)) ;
%                 subplot_tight(6,4,ind,[space_v space_h])
%                 imshow(imadjust(illum_prof_confc{j,f,ch}),[], 'InitialMagnification', 'fit');
%                 ind = ind+1 ;
            end
        end
    end
    figure(5)
    space_h = 0.01 ;
    space_v = 0.01 ;
    for ch = 1:numel(uniqueChannels)
        for j = 1:numel(imgCondiFolders) ; 
            illum_prof_conc{j,ch} = uint16(mean(cat(3,illum_prof_confc{j,:,ch}),3)) ;            
        end
        illum_prof_c{ch} = uint16(mean(cat(3,illum_prof_conc{:,ch}),3)) ;
%         illum_prof_c{ch} = imopen(illum_prof_c{ch},strel('disk',params.disksize)) ;
        subplot_tight(2,2,ch,[space_v space_h])
        imshow(illum_prof_c{ch},[], 'InitialMagnification', 'fit');
    end                
    
    for j = 1:numel(imgCondiFolders) ; 
        save(fullfile(params.outputImgsPath, imgCondiFolders{j},'MultiplexImageData.mat'),'illum_prof_c','-append');
    end
end
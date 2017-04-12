function align_images(params)
    tophat_opt = params.tophat ;
    imgCondiFolders = readDirSubfolders(params.outputImgsPath,'all');    
    for j = 1:numel(imgCondiFolders) ;
        if tophat_opt == 1
            load(fullfile(params.outputImgsPath, imgCondiFolders{j}, 'MultiplexImageTopHat.mat'));
            uniqueFields = [1:numel(imgsProjTopHat)] ;            
            imgsProjAligned = corr_align(imgsProjTopHat,imgRoundNames) ;
            save(fullfile(params.outputImgsPath, params.Condi,'MultiplexImageTopHatAligned.mat'),'imgsProjAligned',...
                'imgRoundNames','params');
        else
            load(fullfile(params.outputImgsPath, imgCondiFolders{j}, 'MultiplexImageData.mat'));
            uniqueFields = [1:numel(imgsProj)] ;            
            imgsProjAligned = corr_align(imgsProj,imgRoundNames) ;
            save(fullfile(params.outputImgsPath, params.Condi,'MultiplexImageDataAligned.mat'),'imgsProjAligned',...
                'imgRoundNames','params');
        end
    end 
end
function imgsProjAligned = corr_align(imgsProj,imgRoundFolders)    
%% Perform alignment of projected images based on cross-correlation of MAP2 
    uniqueFields = [1:numel(imgsProj)] ;  
    uniqueChannels = [1:size(imgsProj{1},1)] ; 
    corrMtxZpeakShift = repmat({repmat({zeros(2,1)},numel(imgRoundFolders),1)},numel(uniqueFields),1);
    for f = 1:numel(uniqueFields)
        round_start = 1 ;
        imA = imgsProj{f}{1,round_start};
        while ~any(imA(:)) % find the first non-zeros round as the alignment reference
            corrMtxZpeakShift{f}{round_start,1}(1) = 0;
            corrMtxZpeakShift{f}{round_start,1}(2) = 0;
            round_start = round_start+1 ;
            imA = imgsProj{f}{1,round_start};
        end
        for roundFolder = round_start+1:numel(imgRoundFolders)
            imA = imgsProj{f}{1,1};
            imB = imgsProj{f}{1,roundFolder};
            % No alignment if any of the images is all zeros or empty 
            if ~any(imA(:))||~any(imB(:))||isempty(imA)||isempty(imB)
                corrMtxZpeakShift{f}{roundFolder,1}(1) = 0;
                corrMtxZpeakShift{f}{roundFolder,1}(2) = 0;
                if ~any(imB(:))
                    imgsProj{f}{1,roundFolder} = imB+1 ; % add dummy value to all-zero images
                end
            else
                C = normxcorr2(imA,imB); % cross-correlation function
                Cmax = max(C(:)) ;
                if Cmax <0.1
                    warning(['correlation between field ', num2str(f), ' ',...
                        imgRoundNames{1}, ' and ', imgRoundNames{roundFolder},...
                        ' < 0.1, no registration is performed'])
                else                    
                    [ypeak,xpeak] = find(C==Cmax);
                    % save the shift of CCF peak
                    corrMtxZpeakShift{f}{roundFolder,1}(1) = xpeak - size(imA,2);
                    corrMtxZpeakShift{f}{roundFolder,1}(2) = ypeak - size(imA,1);
                end
            end
        end
        disp(['Field ',num2str(f),' of ',num2str(numel(uniqueFields)),' complete.']);
    end
    %%
    %Align the images based on the shift of CCF peak from the 1st round
    imgsProjAligned = imgsProj;
    for f = 1:numel(uniqueFields)
        %%           
        for roundFolder = 1:numel(imgRoundFolders)
            imA = imgsProj{f}{1,1};
            imB = imgsProj{f}{1,roundFolder};
            [nRows,nCols] = size(imA); 
            r = corrMtxZpeakShift{f}{roundFolder,1}(2);
            c = corrMtxZpeakShift{f}{roundFolder,1}(1);
            if r > 0 %if padded at the bottom
                typePadRows = 'post';
                idxRemRows = 1:r; %if padding at the bottom, remove same number of image's pixels from top
            elseif r < 0 %if padded at the top
                typePadRows = 'pre';
                idxRemRows = (nRows-abs(r)+1):nRows; %if padding at the top, remove same number of image's pixels from bottom
            end
            if c > 0 %if padded at the right
                typePadCols = 'post';
                idxRemCols = 1:c; %if padding at the right, remove same number of image's pixels from left
            elseif c < 0 %if padded at the left
                typePadCols = 'pre';
                idxRemCols = (nCols-abs(c)+1):nCols; %if padding at the left, remove same number of image's pixels from right
            end

            for chan = 1:numel(uniqueChannels)

                Itemp = imgsProj{f}{chan,roundFolder};

                if r~=0
                    Itemp = padarray(Itemp,[abs(r) 0],typePadRows);
                end
                if c~=0
                    Itemp = padarray(Itemp,[0 abs(c)],typePadCols);
                end

                [nRows,nCols] = size(Itemp); 
                if r > 0 %if padded at the bottom
                    idxRemRows = 1:r; %if padding at the bottom, remove same number of image's pixels from top
                elseif r < 0 %if padded at the top
                    idxRemRows = (nRows-abs(r)+1):nRows; %if padding at the top, remove same number of image's pixels from bottom
                end
                if c > 0 %if padded at the right
                    idxRemCols = 1:c; %if padding at the right, remove same number of image's pixels from left
                elseif c < 0 %if padded at the left
                    idxRemCols = (nCols-abs(c)+1):nCols; %if padding at the left, remove same number of image's pixels from right
                end

                if r~=0
                    Itemp(idxRemRows,:) = [];
                end
                if c~=0
                    Itemp(:,idxRemCols) = [];
                end

%                     figure;
%                     subplot(1,2,1)
%                     imshowpair(imB,imA)
%                     subplot(1,2,2)
%                     imshowpair(Itemp,imA)
%                    
                %Save aligned image
                imgsProjAligned{f}{chan,roundFolder} = Itemp;
            end
        end

        %Crop the aligned images so they are all aligned
        Icombo = cat(3,imgsProjAligned{f}{1,:});
        IcomboMin = min(Icombo(:,:,2:end),[],3);
        r0 = find(sum(IcomboMin,2)==0);
        c0 = find(sum(IcomboMin,1)==0)';

        for roundFolder = 1:numel(imgRoundFolders)
            for chan = 1:numel(uniqueChannels)
                if ~isempty(imgsProjAligned{f}{chan,roundFolder})
                    imgsProjAligned{f}{chan,roundFolder}(r0,:) = [];
                    imgsProjAligned{f}{chan,roundFolder}(:,c0) = [];
                end
            end
        end
    end
end
function Parse_align_images()

        
 %-----------------------------------------------------------------------------------------------------------------------
% opt_corp: export the full field of view plus manually cropped regions (1) or the full feild of views only (0)

% Path of the parent folder for current experiment and analysis
params.parentFolderForAnalysis = 'E:\data\Neuron\cortical\Broad_HCS\14days\20160802-CrossTalk-neurons' ;
cd(params.parentFolderForAnalysis) ;

%-----------------------------------------------------------------------------------------------------------------------
%Set the name of the folder that contains the folders of images of all fields per target. 
params.inputImgsPath = 'E:\data\Neuron\cortical\Broad_HCS\14days\20160802-CrossTalk-neurons' ;

%-----------------------------------------------------------------------------------------------------------------------
%Set the type of projection to do on the z-stack for each target and field of view, as 'max' or 'average'
params.projType = 'average' ;

%-----------------------------------------------------------------------------------------------------------------------
%Set the name of the output folder that will contains the pre-processed images. Within this folder will be folders per target, and within per-target folders,
%folders for different fields of view
params.outputImgsPath = 'E:\data\Neuron\cortical\Broad_HCS\14days\20160802-CrossTalk-neurons\post_processing' ;

%---------------------------------------------------------------------------------------------------
%List the numbers of the channels that are imaged. For instance, if 4
%channel are images set to '0123'. If 3 channels are imaged set to
%'012', etc. Check the raw image file names to see the channel label
%numbers output by Cellomics. These numbers are based on Cellomics
%designation (i.e. channels start at zero).

channelLabels = '1234';
numPlanes = 2 ;

%Set the channel names of the nucleus, MAP2, and LNA target, in order of channel number.
params.channels.cellimgs = {'MAP2','DAPI','LNA'} ;
    %Add path of the folder that contains all relevant scripts
        addpath(genpath('E:\Google Drive\my code\Phenix_analysis'));

    %Read in the parameter file 
%         params = readParamsFromFile(paramsFilePath);
    
    %Create output folder path if it doesn't exist, that is called "Raw
    %Images" where corresponding images in each field of view will be
    %placed in their own folders
        if ~exist(params.outputImgsPath,'dir')
            mkdir(params.outputImgsPath);
        end
        
    %Read in all the folders of images
        imgRoundFolders = readDirSubfolders(params.inputImgsPath,'all');
        
    %Extract the imaging round names from the folders
        toksRoundNames = cellfun(@(x) regexp(x,'LL20160427-21DIV-DNA-PAINT-\d+-(.*)__2016','tokens'),...
            imgRoundFolders,'uniformoutput',0);
        nonempty_ind = cellfun(@(x) ~isempty(x), toksRoundNames,'uniformoutput',0) ;
        nonempty_ind = cat(1,nonempty_ind{:}) ;
        toksRoundNames = toksRoundNames(nonempty_ind) ;
        imgRoundNames = cellfun(@(x) cat(1,x{1}{1}),toksRoundNames,'uniformoutput',0);
        nonwash_ind = cellfun(@(x) isempty(strfind(x,'wash')), imgRoundNames,'uniformoutput',0) ;        
        nonwash_ind = cat(1,nonwash_ind{:}) ;
        imgRoundFolders =imgRoundFolders(nonempty_ind) ;
        imgRoundFolders =imgRoundFolders(nonwash_ind) ;
        imgRoundNames =imgRoundNames(nonwash_ind) ;
        imgRoundFolders(1) = [] ; % delete the "no-probe" round
        imgRoundNames(1) = [] ;
        
        
    %Create field folders and write images after projection
        roundFolder = 1 ;
        roundFolderImgs = readDirImages(fullfile(params.inputImgsPath,imgRoundFolders{roundFolder},'Images'),'TIFF',0);
            %Extract: field, z-stack number, and channel number for each
            %image in current imaging round folder
%             toksCellImgs = cellfun(@(x) regexp(x,'f(.*)p(.*)-ch(.*)sk','tokens'),...
%                 roundFolderImgs.filenames,'uniformoutput',0);
            toksCellImgs = cellfun(@(x) regexp(x,['r(\d+)c(\d+)f(\d+)p(\d+)-ch([',channelLabels,'])\w*.tiff'],'tokens'),....
                roundFolderImgs.filenames, 'uniformoutput',0);  %% filename parser for Phenix
            imgInfo.rowLabel = cellfun(@(x) cat(1,x{1}{1}),toksCellImgs,'uniformoutput',0);
            imgInfo.colLabel = cellfun(@(x) cat(1,x{1}{2}),toksCellImgs,'uniformoutput',0) ;
            uniqueRowLabels = unique(imgInfo.rowLabel);
            uniqueColLabels = unique(imgInfo.colLabel);
            
            %Read in the Excel plate layout and name designations, and assign to each
            %image the corresponding labels.
            plateFieldNames = {'Condi','BlueChan','GreenChan','RedChan','FarRedChan',...
                'RowLabels','ColLabels'};
            for i = 1:numel(plateFieldNames)
                eval(['[~,plate.',plateFieldNames{i},'] = xlsread(fullfile(params.parentFolderForAnalysis,''plateLayout_MATLAB.xlsx''),''',plateFieldNames{i},''');']);
            end
for r = 1:numel(uniqueRowLabels)
    for c = 2:numel(uniqueColLabels)
        Condi = plate.Condi{str2num(uniqueRowLabels{r}), str2num(uniqueColLabels{c})};
        params.Condi = Condi ;
        parpool('local',10)
        for roundFolder = 1:numel(imgRoundFolders)
            %% Read in all image names in current imaging round folder
            roundFolderImgs = readDirImages(fullfile(params.inputImgsPath,imgRoundFolders{roundFolder},'Images'),'TIFF',0);
            %Extract: field, z-stack number, and channel number for each
            %image in current imaging round folder
%             toksCellImgs = cellfun(@(x) regexp(x,'f(.*)p(.*)-ch(.*)sk','tokens'),...
%                 roundFolderImgs.filenames,'uniformoutput',0);
            toksCellImgs = cellfun(@(x) regexp(x,['r(\d+)c(\d+)f(\d+)p(\d+)-ch([',channelLabels,'])\w*.tiff'],'tokens'),....
                roundFolderImgs.filenames, 'uniformoutput',0);  %% filename parser for Phenix
            imgInfo.rowLabel = cellfun(@(x) cat(1,x{1}{1}),toksCellImgs,'uniformoutput',0);
            imgInfo.colLabel = cellfun(@(x) cat(1,x{1}{2}),toksCellImgs,'uniformoutput',0) ;
            imgInfo.fieldNum = cellfun(@(x) cat(1,x{1}{3}),toksCellImgs,'uniformoutput',0);
            imgInfo.planeNum = cellfun(@(x) cat(1,x{1}{4}),toksCellImgs,'uniformoutput',0);
            imgInfo.channelNum = cellfun(@(x) cat(1,x{1}{5}),toksCellImgs,'uniformoutput',0);
            
            for i = 1:numel(roundFolderImgs.filenames)
                imgInfo.comboLabel{i,1} = [imgInfo.rowLabel{i,1},'_',...
                     imgInfo.colLabel{i,1},'_',imgInfo.fieldNum{i,1}];
            end                       
            
            %% Create a folder for each field, if it doesn't exist
            uniqueFields = unique(imgInfo.fieldNum) ;
            uniqueChannels = unique(imgInfo.channelNum);
            
            if roundFolder == 1
                imgsProj = cell(numel(params.channels.cellimgs),numel(imgRoundFolders)) ;
                imgsProj = repmat({imgsProj},numel(uniqueFields),1); %# fields, rows = channels, cols = acquisition
%                 imgsProj = repmat({imgsProj},numel(uniqueRowLabels),numel(uniqueColLabels));
                imgsCombo = imgsProj ;
            else
            end
            
            
                    Condi = plate.Condi{str2num(uniqueRowLabels{r}), str2num(uniqueColLabels{c})};
                        if ~exist(fullfile(params.outputImgsPath, Condi),'dir')
                            mkdir(fullfile(params.outputImgsPath, Condi));
                        end
   
                        for f = 1:numel(uniqueFields)
                            field = ['field_',uniqueFields{f}];
                            if ~exist(fullfile(params.outputImgsPath,Condi, field),'dir')
                                mkdir(fullfile(params.outputImgsPath,Condi, field));
                            end
                        end
                      
            %% Load in all images in current imaging round folder
            loc =  all([ismember(imgInfo.rowLabel,uniqueRowLabels{r}),...
                                ismember(imgInfo.colLabel,uniqueColLabels{c}),...                                
                                ],2);
            roundFolderImgsPathCondi = roundFolderImgs.filenamesWithPath(loc);
            fieldNumCondi = imgInfo.fieldNum(loc) ;
            channelNumCondi = imgInfo.channelNum(loc) ;
            imgs = cell(size(roundFolderImgsPathCondi)) ;
            
            parfor pathInd = 1:numel(roundFolderImgsPathCondi)
                imgs{pathInd} = imread(roundFolderImgsPathCondi{pathInd}) ;
            end

%             imgs = cellfun(@(x) loadtiff(x),roundFolderImgsPathCondi,'uniformoutput',0);
            %Make image projections (cell array of arrays, each for a given
            %channel)
            
            
                    for f = 1:numel(uniqueFields)
                        for ch = 1:numel(uniqueChannels)
                            loc = all([ismember(fieldNumCondi,uniqueFields{f}),...
                                ismember(channelNumCondi,uniqueChannels{ch}),...
                                ],2);
                            imgsCombo{f}{ch,roundFolder} = imgs(loc);
                            if strcmp(params.projType,'max')
                                imgsProj{f}{ch,roundFolder} = max(cat(3,imgs{loc}),[],3);    
                            elseif strcmp(params.projType,'average')
                                imgsProj{f}{ch,roundFolder} = uint16(mean(cat(3,imgs{loc}),3)); 
                            end
                        end
                    end                    
        end
        poolobj = gcp('nocreate');
        delete(poolobj);
%         save(fullfile(params.outputImgsPath, Condi,'MultiplexImageData.mat'),...
%             'imgsCombo','imgsProj','imgRoundNames','params');
        %Save the projected image to the correct field folder
        %and channel/target name as multiplage tiff. Each tiff is all
        %rounds for given channel and field of view
         
   %%             
                for f = 1:numel(uniqueFields)
                    for ch = 1:numel(uniqueChannels)
                        tempImgName = fullfile(params.outputImgsPath,Condi ,['field_',uniqueFields{f}],['I_',params.channels.cellimgs{ch},'_proj2.tif']);
                        options.big = true; % Use BigTIFF format
                        options.message = true;
                        saveastiff(cat(3,imgsProj{f}{ch,:}), tempImgName, options); % "saveastiff" is faster than "imwrite"
                    end
                end
          
        
        
        %% Perform alignment of projected images based on cross-correlation of MAP2 
        corrMtxZpeakShift = repmat({repmat({zeros(2,1)},numel(imgRoundFolders),1)},numel(uniqueFields),1);
        for f = 1:numel(uniqueFields)
            for roundFolder = 2:numel(imgRoundFolders)
                imA = imgsProj{f}{1,1};
                imB = imgsProj{f}{1,roundFolder};
                % add 1 to the image if it's all zeros 
                if max(imA(:)) == 0 ;
                    imgsProj{f}{1,1}= imA +1 ;
                end
                if max(imB(:)) == 0 ;
                    imgsProj{f}{1,roundFolder} = imB + 1 ;
                end
                C = normxcorr2(imA,imB);
                Cmax = max(C(:)) ;
                if Cmax <0.1
                    warning(['correlation between field ', num2str(f), ' ',...
                        imgRoundNames{1}, ' and ', imgRoundNames{roundFolder},...
                        ' < 0.1, no registration is performed'])
                else                    
                    [ypeak,xpeak] = find(C==Cmax);
                    %Shift of peak
                    corrMtxZpeakShift{f}{roundFolder,1}(1) = xpeak - size(imA,2);
                    corrMtxZpeakShift{f}{roundFolder,1}(2) = ypeak - size(imA,1);
                end
            end
            disp(['Field ',num2str(f),' of ',num2str(numel(uniqueFields)),' complete.']);
        end
        %%
        %Align the images based on 1st round
        imgsProjAligned = imgsProj;
        for f = 1:numel(uniqueFields)
            %%
            if f == 24
                pause(1);
            end
            for roundFolder = 2:numel(imgRoundFolders)
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
                    imgsProjAligned{f}{chan,roundFolder}(r0,:) = [];
                    imgsProjAligned{f}{chan,roundFolder}(:,c0) = [];
                end
            end

        end
        
        save(fullfile(params.outputImgsPath, params.Condi,'MultiplexImageData.mat'),'imgsProjAligned',...
            'imgRoundNames','corrMtxZpeakShift','params');
        multiplex_confocal_plot(params,0)
        end
    end
end
        
       
        

            
            
            
            
 

        
        
        
        
        
        
        
        
        
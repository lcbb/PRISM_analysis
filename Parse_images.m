function Parse_images(params)
% parse images in Phenix format accroding to the plate map specified ('plateLayout_MATLAB.xlsx')
% images from the same plate well are stacked together
    channelLabels = params.channelLabels ;
    numPlanes = params.numPlanes ;

    %Add path of the folder that contains all relevant scripts
    addpath(genpath('E:\Google Drive\my code\Phenix_analysis'));

    %Read in the parameter file 
    %         params = readParamsFromFile(paramsFilePath);

    %Create output folder path if it doesn't exist
    if ~exist(params.outputImgsPath,'dir')
    mkdir(params.outputImgsPath);
    end

    %Read in all the folders of images
    imgRoundFolders = readDirSubfolders(params.inputImgsPath,'all');

    %Extract the imaging round names from the folders    
    toksRoundNames = cellfun(@(x) regexp(x,params.imgRoundFoldersRegExp,'tokens'),...
    imgRoundFolders,'uniformoutput',0);
    nonempty_ind = cellfun(@(x) ~isempty(x), toksRoundNames,'uniformoutput',0) ;
    nonempty_ind = cat(1,nonempty_ind{:}) ; % remove the folders without matched names
    toksRoundNames = toksRoundNames(nonempty_ind) ;
    imgRoundFolders = imgRoundFolders(nonempty_ind) ;
    imgRoundNames = cellfun(@(x) cat(1,x{1}{2}),toksRoundNames,'uniformoutput',0);
    imgRoundOrder = cellfun(@(x) cat(1,x{1}{1}),toksRoundNames,'uniformoutput',0);
    imgRoundOrder = cellfun(@(x) str2double(x),imgRoundOrder,'uniformoutput',0);
    imgRoundOrder = cell2mat(imgRoundOrder) ;
    [imgRoundOrderSorted, ind] = sort(imgRoundOrder) ;   
    imgRoundNames = imgRoundNames(ind) ; % sort the image round name according to the image order
    imgRoundFolders = imgRoundFolders(ind) ;    
    sub_ind = [1:numel(imgRoundNames)] ;
    imgRoundNames = imgRoundNames(sub_ind) ; % sort the image round name according to the image order
    imgRoundFolders = imgRoundFolders(sub_ind) ;      

    %Read in the Excel plate layout and name designations, and assign to each
    %image the corresponding labels.
    plateFieldNames = {'Condi','BlueChan','GreenChan','RedChan','FarRedChan',...
        'RowLabels','ColLabels'};
    for i = 1:numel(plateFieldNames)
        eval(['[~,plate.',plateFieldNames{i},'] = xlsread(fullfile(params.parentFolderForAnalysis,params.plateMapFile),','plateFieldNames{i}',');']);
    end

    uniqueCondibi = cellfun(@(x) strcmp(x,'N/A'),plate.Condi,'uniformoutput',0) ;
    uniqueCondibi = ~cell2mat(uniqueCondibi) ;
    [uniqueCondiIndR, uniqueCondiIndC] = find(uniqueCondibi) ;

    poolobj = gcp('nocreate') ;
    if isempty(poolobj) % start parpool if it doen't exist 
        parpool('local',params.n_jobs)
    else
    end
    
    for j = 1:numel(uniqueCondiIndR)
%     for c = 1:numel(uniqueColLabels)
        r = uniqueCondiIndR(j) ;
        c = uniqueCondiIndC(j) ;
        Condi = plate.Condi{r, c};
        params.Condi = Condi ;
        if ~exist(fullfile(params.outputImgsPath, Condi),'dir')
            mkdir(fullfile(params.outputImgsPath, Condi));
        end        
         for roundFolder = 1:numel(imgRoundFolders)
            %% Read in all image names in current imaging round folder
            roundFolderImgs = readDirImages(fullfile(params.inputImgsPath,imgRoundFolders{roundFolder},'Images'),'TIFF',0);
            %Extract: field, z-stack number, and channel number for each
            %image in current imaging round folder
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
            
            for f = 1:numel(uniqueFields)
                field = ['field_',uniqueFields{f}];
                if ~exist(fullfile(params.outputImgsPath,Condi, field),'dir')
                    mkdir(fullfile(params.outputImgsPath,Condi, field));
                end
            end                                                                 
            %% Load in all images in current imaging round folder
            loc =  all([ismember(imgInfo.rowLabel,plate.RowLabels(r)),...
                                ismember(imgInfo.colLabel, plate.ColLabels(c)),...                                
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
        %Save the projected image to the correct field folder
        %and channel/target name as multiplage tiff. Each tiff is all
        %rounds for given channel and field of view                        
%         for f = 1:numel(uniqueFields)
%             for ch = 1:numel(uniqueChannels)
%                 tempImgName = fullfile(params.outputImgsPath,Condi ,['field_',uniqueFields{f}],['I_',params.channels.cellimgs{ch},'_proj2.tif']);
%                 options.big = true; % Use BigTIFF format
%                 options.message = true;
%                 saveastiff(cat(3,imgsProj{f}{ch,:}), tempImgName, options); % "saveastiff" is faster than "imwrite"
%             end
%         end
        save(fullfile(params.outputImgsPath, params.Condi,'MultiplexImageData.mat'),'imgsProj',...
                'imgRoundNames','params');
    end            
    poolobj = gcp('nocreate');
    delete(poolobj); % shut down parpool 
end




        
       
        

            
            
            
            
 

        
        
        
        
        
        
        
        
        
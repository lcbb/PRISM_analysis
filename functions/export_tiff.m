cd(params.outputImgsPath) ;
pccCellName = {'pccIabInMarkerMean', 'pccIabWithMarkerMean', 'pccAabWithMarkerMean'} ;
dirList = dir(params.outputImgsPath) ;
isub = [dirList(:).isdir]; % returns logical vector
dirList = {dirList(isub).name}';
dirList(ismember(dirList,{'.','..'})) = [];
normBool = 1;
thsize = 100;
for k = 1:numel(dirList)    
    load(fullfile(params.outputImgsPath, dirList{k}, 'MultiplexImageData.mat'))           
        uniqueFields = [1:numel(imgsProjAligned)] ;
        imgRoundNames = [imgRoundNames; {'DAPI'; 'MAP2'}] ;
        targets  = [1:numel(imgRoundNames)] ;
        
        if ~exist(fullfile(params.outputImgsPath,params.Condi, 'TIFF'),'dir')
            mkdir(fullfile(params.outputImgsPath,params.Condi, 'TIFF'));
        end
        for f = 1:numel(uniqueFields)            
            sizeI = size(imgsProjAligned{f}{1,1});
            Idapi = cat(3,zeros(sizeI),zeros(sizeI),imadjust(imgsProjAligned{f}{2,1}));
            IMAP2 = cat(3,zeros(sizeI),imadjust(imgsProjAligned{f}{1,1}),zeros(sizeI));
            ItargetStack = cat(3,imgsProjAligned{f}{3,targets(1:end-2)});
            ItargetStack = cat(3,ItargetStack,imgsProjAligned{f}{[2,1],1});
            if normBool 
                for i = 1:size(ItargetStack,3)
                    ItargetStack(:,:,i) = imtophat(ItargetStack(:,:,i),strel('disk',thsize));
                end
            end                                                             
            tempImgName = fullfile(params.outputImgsPath,params.Condi ,'TIFF',['field_',num2str(uniqueFields(f)),'_aligned.tif']);
            options.big = true; % Use BigTIFF format
            options.message = true;
            saveastiff(ItargetStack, tempImgName, options); % "saveastiff" is faster than "imwrite"
        end
end
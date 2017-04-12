function cell_mask_batch(params)
    imgRoundFolders = readDirSubfolders(params.outputImgsPath,'all');
    plot_opt = params.plot_opt ;
    for j = 1:numel(imgRoundFolders)
        if params.tophat == 1
            load(fullfile(params.outputImgsPath,imgRoundFolders{j}, 'MultiplexImageTopHatAligned.mat')) ;           
            params.Condi = imgRoundFolders{j} ;
            params.plot_opt =  plot_opt;
            cell_mask(imgsProjAligned, params, imgRoundNames) 
        else
            load(fullfile(params.outputImgsPath,imgRoundFolders{j}, 'MultiplexImageDataAligned.mat')) ;             
            params.Condi = imgRoundFolders{j} ;
            params.plot_opt = plot_opt ;
            cell_mask(imgsProjAligned, params, imgRoundNames) 
        end            
    end
end
%-----------------------------------------------------------------------------------------------------------------------
% Path of the parent folder for current experiment and analysis
% params.parentFolderForAnalysis = 'E:\data\Neuron\cortical\Broad_HCS\14days\20161027-Rep3 and Rep4\Rep4' ;
params.parentFolderForAnalysis = 'E:\data\Neuron\cortical\Broad_HCS\14days\20161021-crotalk in cells-t9-Rep2' ;
cd(params.parentFolderForAnalysis) ;

%-----------------------------------------------------------------------------------------------------------------------
%Set the name of the folder that contains the folders of images of all fields per target. 
params.inputImgsPath = params.parentFolderForAnalysis ;

%-----------------------------------------------------------------------------------------------------------------------
%Set the type of projection to do on the z-stack for each target and field of view, as 'max' or 'average'
params.projType = 'average' ;

%-----------------------------------------------------------------------------------------------------------------------
%Set the name of the output folder that will contains the pre-processed images. Within this folder will be folders per target, and within per-target folders,
%folders for different fields of view
params.outputImgsPath = fullfile(params.inputImgsPath, 'post_processing_multicolor') ;
params.plateMapFile = 'plateLayout_multicolor.xlsx' ;
% params.outputImgsPath = fullfile(params.inputImgsPath, 'post_processing_crosstalk') ;
% params.plateMapFile = 'plateLayout_crosstalk.xlsx' ;
%---------------------------------------------------------------------------------------------------
%List the numbers of the channels that are imaged. For instance, if 4
%channel are images set to '0123'. If 3 channels are imaged set to
%'012', etc. Check the raw image file names to see the channel label
%numbers output by Cellomics. These numbers are based on Cellomics
%designation (i.e. channels start at zero).

params.channelLabels = '1234';

%Set the channel names of the nucleus, MAP2, and LNA target, in order of channel number.
params.channels.cellimgs = {'MAP2','DAPI','LNA','vGlut'} ;
% params.channels.cellimgs = {'DAPI','LNA'} ;
%---------------------------------------------------------------------------------------------------
% set the regular expression for parsing the image folder names. Tokens are 
% LNA sequences "p#" 
% params.imgRoundFoldersRegExp = '20160802-384-DNA-PAINT-(.*)-(\d+)__2016'
% params.imgRoundFoldersRegExp = 'LL20160902-DNA-PAINT-(\d+)-(.*)__2016'    
% params.imgRoundFoldersRegExp = 'LL20160909-p1-(\d+)-(.*)__2016'    
% params.imgRoundFoldersRegExp = 'LL20161027-DNA-PAINT-N2-(\d+)-(.*)__2016' ;
params.imgRoundFoldersRegExp = 'LL20161021-DNA-PAINT-(\d+)-(.*)__2016' ;
params.n_jobs = 20 ; % # of paralle jobs
params.pixsize = 0.187 ; % pixel size in um
params.disksize = 100;  % disk size for illumination fuction estimation and tophat filtering   
%%---------------------------------------------------------------------------------------------------
% parse the images from the same well into individual mat file 
Parse_images(params)
% compute the illuminatin profile
compute_illum_profile(params)
% correct uneven illumination and appy tophat filter
tophat_batch(params)
% align images 
params.tophat = 0 ;
align_images(params)
%% display images
% params.plateMapFile = 'plateLayout_crosstalk.xlsx' ;
PlotRoundInd = [1:28] ;
% PlotRoundInd = [1:2] ;
% CondiInd = [1:3, 6:7, 10:14] ;
% sat_frc = [0/65535,1000/65535; 0/65535,13000/65535; 0/65535,2500/65535] ;
% CondiInd = [1:5] ;
% CondiInd = [5:8] ;
CondiInd = [1:10] ;
ChanInd = [1,2,3] ;
% sat_frc = [0.1, 0.97; 0.1, 0.95; 0.1, 0.95] ; % pixel saturation fraction for RGB channels
sat_frc = [0.1, 0.99] ; % pixel saturation fraction for RGB channels
params.tophat = 1 ; % apply tophat filter if 1
params.align = 1 ;
%
multiplex_confocal_plot(params,sat_frc, 1)
% multi_condition_confocal_plot(params, CondiInd, PlotRoundInd, ChanInd, sat_frc)
%% synapse segmentation 
params.plot_opt = 1 ;
cell_mask_batch(params)
%% feature extraction 
synapse_feature(params)

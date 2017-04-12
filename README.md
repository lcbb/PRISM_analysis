PRISM analysis - program to process and analyze [PRISM](http://www.biorxiv.org/content/early/2017/02/25/111625) images of neuronal synapses 

This program takes the PRISM images acquired using Perkin Elmer Opera Phenix High-Content Screening System as input, 
1. performs drift and flatfield corrections
2. parse images from the same well and save as individual MAT files
3. segment the images and save the masks
4. extract features of synapses 

Example image data will be available soon

# Usage
Change "params.parentFolderForAnalysis" in "main.m" to the location of the image folders, and run "main.m".
Processed images will be saved at "params.outputImgsPath"

######## Output Data structure ############
Each "Rep#" folder contains one biological replicate (neuronal culture, plate, or batch), total 3 batches. 
Each sub-folder "RepX-Y" contains one replicate (well) within each biological replicate, total 3 wells. 
	MultiplexImageData.mat: parsed, unprocessed images
	
		imgsProj: collection of aligned, top-hat filtered images from a single well, 
		has the structure imgsProjTopHat{field of view}{channel, imaging round}, 
		channel order:'MAP2','DAPI','LNA','vGlut'
		imaging round: see variable "imgRoundNames"
		
		imgRoundNames: list of names of incubated probe sequences (p#), wash-out (WO) or not, exposure time (#S). 
		
		List of correspoding protein targets: {'actin','p2';'Tuj-1','p3';'cortactin','p4';...
        'Shank3','p6';'ARPC2','p7'; 'bassoon','p8'; 'synapsin-1(2)','p9-t2'; 'synapsin-1','p9';...
        'Homer-1b/c','p10';'NR2B','p12';'PSD-95','p1-'} ;
		
	MultiplexImageDataAligned.mat: processed images after drift correction and flatfiel correction.
	
	MultiplexImageAlignedTopHat.mat: same as "MultiplexImageDataAligned.mat" but with top-hat filtering.
PRISM analysis - program to process and analyze [PRISM](http://www.biorxiv.org/content/early/2017/02/25/111625) images of neuronal synapses 

This program takes the PRISM images acquired using Perkin Elmer Opera Phenix High-Content Screening System as input, 
1. performs drift and flatfield corrections
2. parse images from the same well and save as individual MAT files
3. segment the images and save the masks
4. extract features of synapses 

Example image data will be available soon

# Usage
change "params.parentFolderForAnalysis" in "main.m" to the location of the image folders, and run "main.m"
processed images will be saved at "params.outputImgsPath"
#
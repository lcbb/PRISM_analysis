function  cell_mask(imgsProjAligned, params, imgRoundNames) 
% gernerate sub cellular region masks for PRISM data

global fig_num fig_path;
close all
uniqueFields = [1:numel(imgsProjAligned)] ;
ch_name = imgRoundNames ;        
maskfile = fullfile(params.outputImgsPath,params.Condi, 'mask.mat') ;
cellmask =  cell(5,numel(imgRoundNames)) ; % generate masks for different cell regions for each round 
cellmask = repmat({cellmask},numel(uniqueFields),1); 
synMask = cell(numel(uniqueFields),1); % synapse masks
[imgRoundNames, targets]= ParseImgName(imgRoundNames) ;
% imgRoundNames = [imgRoundNames(targets); {'VGluT';'DAPI'; 'MAP2'}] ;
imgRoundNamesNoWO = [imgRoundNames(targets); 'VGLUT1'; 'MAP2'] ;

%%
% params.maxObjSize = round(15^2*pi);
params.nucRm = 1 ; %remove the nuclei or not
params.pixsize = 0.187 ; % pixel size in um
for f = 1:numel(uniqueFields)    
    disp(['process ', params.Condi,  ' field ', num2str(f), '...'])    
    ItargetStack = cat(3,imgsProjAligned{f}{3,targets});
    ItargetStack = cat(3,ItargetStack,imgsProjAligned{f}{4,1}); % Add VGlutT channel to target stack
    ItargetStack = cat(3,ItargetStack,imgsProjAligned{f}{1,1}); % Add MAP2 channel to target stack
    sizeI = size(imgsProjAligned{f}{1,1});
    synMaskAll = [] ;
    params.plot_opt = 0 ; % turn off plotting 
    for roundInd = 1:numel(imgRoundNames)        
        Idapi = imgsProjAligned{f}{2,roundInd};
        IMAP2 = imgsProjAligned{f}{1,roundInd};
        IVGluT = imgsProjAligned{f}{4,roundInd};
        if isempty(Idapi)||isempty(IMAP2)||isempty(IVGluT)...
                ||~any(Idapi(:))||~any(IMAP2(:))||~any(IVGluT(:)) % if any channel is empty or all zero, leave the mask empty
        else                     
            params = segmentNuc(Idapi, params) ;
            params = segmentMAP2(IMAP2, params) ;                
            params = segmentDendrite(params) ;
            params.ch_name = 'VGluT' ;
            params = segmentSynapse(IVGluT, params) ;
            params = segmentCell(params) ;                

            cellmask{f}{1,roundInd} = params.NucMask ;
            cellmask{f}{2,roundInd} = params.MAP2Mask ;
            cellmask{f}{3,roundInd} = params.DendriteMask ;
            cellmask{f}{4,roundInd} = params.SynMask ;
            cellmask{f}{5,roundInd} = params.CellMask ;
        end
    end
    params.plot_opt = 1 ; % turn on plotting
    % Stack synapse maskes from non-wash-out rounds
    for roundInd = 1:numel(imgRoundNamesNoWO)                                                 
        params.ch_name = imgRoundNamesNoWO(roundInd);
        params = segmentSynapse(ItargetStack(:,:,roundInd), params) ;
        synMaskAll = cat(3, synMaskAll, params.SynMask) ;                           
    end
        synMask{f} = synMaskAll ;      
end
save(maskfile, 'cellmask', 'synMask','params', 'targets', 'imgRoundNames','imgRoundNamesNoWO')   
end


function params = segmentSynapse(Iab_eq, params)
%% Do Otsu thresholding on marker and ab.
global fig_num
global fig_path

if params.nucRm == 1 ;
    IabOtsuMaskNuc = params.NucMask ;        
    sat_lim = stretchlim(Iab_eq(~IabOtsuMaskNuc)) ;    
    Iab_ad = imadjust(Iab_eq, sat_lim);
    Iab_eq(IabOtsuMaskNuc) =  0 ;
else
    Iab_ad = imadjust(Iab_eq);
%     Iab_ad = Iab_eq ;
end
Iab_denoise = wiener2(Iab_eq,[5 5]);
Iab_tophat = imtophat(Iab_denoise,strel('disk',8));
IabMask = ad_thresh(Iab_tophat, params);
% IabOtsuMask = imopen(IabMask,strel('disk',8/params.imBin));
% IabOtsuMask = bwareaopen(IabOtsuMask, params.maxObjSize);

%%
IabMask_perim = bwperim(IabMask);
Iab_ad_perim = imoverlay(Iab_ad, IabMask_perim, [1 .3 .3]);

% watershed 
%%watershed doesn't work well for mask of small objects. Apply wartershed to grayscale images instead 
Iab_maxs = imextendedmax(Iab_tophat,  0.001);
Iab_overlay_max = imoverlay(Iab_ad,Iab_maxs, [0.3 1 .3]);

Iab_com = imcomplement(Iab_tophat);
Iab_com = imimposemin(Iab_com, Iab_maxs);

Iab_shed = watershed(Iab_com);
Iab_shed(Iab_shed==1) = 0; %make background zero.
Iab_shed(Iab_shed>1) = 1; %make foreground one.
Iab_shed = logical(Iab_shed); %convert watershed mask to logical.
Iab_shed = imreconstruct(Iab_maxs,Iab_shed);
abMask = Iab_shed & IabMask ;
params.SynMask = abMask ;

IabMask_perim = bwperim(abMask);
Iab_shed_perim = imoverlay(Iab_ad, IabMask_perim, [1 .3 .3]);
%%
    if params.plot_opt == 1
        figure(5)
        format_fig2(13)
        space_h = 0.01 ;
        space_v = 0.01 ;
        
        subplot_tight(2,2,1,[space_v space_h])
        imshow(Iab_ad,[], 'InitialMagnification', 'fit');
        hsb = scalebar('scalelength_um', 20, 'pixelsize', params.pixsize, 'location', 'southeast') ;
        h1 = title('Raw image' , 'fontsize', 8) ;  

        figure(5)
        subplot_tight(2,2,2,[space_v space_h])
        imshow(imadjust(Iab_tophat),[], 'InitialMagnification', 'fit');
        hsb = scalebar('scalelength_um', 20, 'pixelsize', params.pixsize, 'location', 'southeast') ;
        h1 = title('Filtered image' , 'fontsize', 8) ;  
        
        figure(7)
        imshow(IabMask,[], 'InitialMagnification', 'fit');
        h1 = title(params.ch_name , 'fontsize', 6) ;  set(h1,'interpreter','none')

        figure(9)
        % subplot(2,2,3)
        imshow(Iab_ad_perim,[], 'InitialMagnification', 'fit');

        figure(5)
        subplot_tight(2,2,3,[space_v space_h])
        imshow(abMask,[], 'InitialMagnification', 'fit');
        hsb = scalebar('scalelength_um', 20, 'pixelsize', params.pixsize, 'location', 'southeast') ;
        h1 = title('Segmented synapses' , 'fontsize', 8) ;  

        figure(5)
        subplot_tight(2,2,4,[space_v space_h])
        imshow(Iab_shed_perim,[], 'InitialMagnification', 'fit');
        hsb = scalebar('scalelength_um', 20, 'pixelsize', params.pixsize, 'location', 'southeast') ;
        h1 = title('Raw image+\color{red}segmentation' , 'fontsize', 8) ;
        
%         ha = axes('Position',[0 0 1 1.05],'Xlim',[0 1],'Ylim',[0 1],...
%             'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
%         text(0.5, 1,['\bf', params.ch_name],'HorizontalAlignment'... 
%         ,'center','VerticalAlignment', 'top', 'color','k', 'fontsize', 15)
        suptitle(params.ch_name, 12)                

%         print(5 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ; save('.\startup.mat', 'fig_num');
        print(5 ,'-dpdf','-r300', [fig_path,num2str(fig_num, '%03d'),'.pdf']) ; fig_num = fig_num +1 ; save('.\startup.mat', 'fig_num');

        figure(11)
        imshow(Iab_overlay_max,[], 'InitialMagnification', 'fit');

        figure(12)
        imshow(Iab_com,[], 'InitialMagnification', 'fit');                     

    end

end

function params = segmentCell(params)
global fig_num
global fig_path
MAP2Mask = imdilate(params.MAP2Mask,strel('disk',5)) ; 
NucMask = imdilate(params.NucMask,strel('disk',10)) ;
SynMask = imdilate(params.SynMask,strel('disk',5)) ;
CellMask = (MAP2Mask +NucMask+SynMask)>0;
params.CellMask = CellMask ;
if params.plot_opt == 1
    figure(104)
    subplot_tight(2,2,1,[0.05 0.01])
    imshow(MAP2Mask,[], 'InitialMagnification', 'fit');

    figure(104)
    subplot_tight(2,2,2,[0.05 0.01])
    % figure(61)
    imshow(NucMask,[], 'InitialMagnification', 'fit');

    figure(104)
    subplot_tight(2,2,3,[0.05 0.01])
    imshow(SynMask,[], 'InitialMagnification', 'fit');
    
    figure(104)
    subplot_tight(2,2,4,[0.05 0.01])
    imshow(CellMask,[], 'InitialMagnification', 'fit');
   
    tightfig();
    print(103 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ;
end

end
function params = segmentDendrite(params)
global fig_num
global fig_path
MAP2Mask = imdilate(params.MAP2Mask,strel('disk',5)) ; 
NucMask = imdilate(params.NucMask,strel('disk',10)) ;
DendriteMask = (MAP2Mask - NucMask)>0;
params.DendriteMask = DendriteMask ;
if params.plot_opt == 1
    figure(103)
    subplot_tight(1,3,1,[0.05 0.01])
    imshow(MAP2Mask,[], 'InitialMagnification', 'fit');

    figure(103)
    subplot_tight(1,3,2,[0.05 0.01])
    % figure(61)
    imshow(NucMask,[], 'InitialMagnification', 'fit');

    figure(103)
    subplot_tight(1,3,3,[0.05 0.01])
    imshow(DendriteMask,[], 'InitialMagnification', 'fit');
   
    tightfig();
    print(103 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ;
end

end

function params = segmentMAP2(Iab, params)
%%
global fig_num
global fig_path
% Iab_eq = adapthisteq(Iab);   
Iab_eq = imadjust(Iab);   


Tab = multithresh(Iab_eq,2);

% IabOtsuMaskNuc = im2bw(Iab_eq,1.4*Tab); %segment bright nuclei
% IabOtsuMaskNuc = im2bw(Iab_eq,Tab(1)); %segment bright nuclei
IabOtsuMaskNuc = imquantize(Iab_eq,Tab); %segment bright nuclei
IabOtsuMaskNuc = IabOtsuMaskNuc>2 ;
% IabOtsuMaskNuc = im2bw(Iab_eq,1*Tab); %segment bright nuclei
% IabOtsuMaskNuc = imopen(IabOtsuMaskNuc,strel('disk',3));
% Tab_open = graythresh(Iab_open);

% IabOtsuMaskNuc = bwareaopen(IabOtsuMaskNuc, params.maxObjSize);
% IabOtsuMaskNuc = imclose(IabOtsuMaskNuc,strel('disk',3));
% IabOtsuMaskNuc = imdilate(IabOtsuMaskNuc,strel('disk',8));
% IabOtsuMaskNuc = imfill(IabOtsuMaskNuc,'holes');

params.MAP2Mask = IabOtsuMaskNuc ;

if params.plot_opt == 1
    figure(101)
    subplot_tight(2,2,1,[0.05 0.01])
    imshow(Iab_eq,[], 'InitialMagnification', 'fit');

    figure(101)
    subplot_tight(2,2,2,[0.05 0.01])
    % figure(61)
    imshow(IabOtsuMaskNuc,[], 'InitialMagnification', 'fit');

    figure(101)
    subplot_tight(2,2,3,[0.05 0.01])
    imshow(IabOtsuMaskNuc,[], 'InitialMagnification', 'fit');

    figure(101)
    subplot_tight(2,2,4,[0.05 0.01])
    imshow(IabOtsuMaskNuc,[], 'InitialMagnification', 'fit');
    tightfig();
    print(101 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ;
end
end

function params = segmentNuc(Iab, params)
%%
global fig_num
global fig_path
% Iab_eq = adapthisteq(Iab);   
Iab_eq = imadjust(Iab);   


Tab = multithresh(Iab_eq,2);

% IabOtsuMaskNuc = im2bw(Iab_eq,1.4*Tab); %segment bright nuclei
% IabOtsuMaskNuc = im2bw(Iab_eq,Tab(1)); %segment bright nuclei
IabOtsuMaskNuc = imquantize(Iab_eq,Tab(1)); %segment bright nuclei
IabOtsuMaskNuc = IabOtsuMaskNuc>1 ;
% IabOtsuMaskNuc = im2bw(Iab_eq,1*Tab); %segment bright nuclei
IabOtsuMaskNuc = imopen(IabOtsuMaskNuc,strel('disk',3));
% Tab_open = graythresh(Iab_open);

% IabOtsuMaskNuc = bwareaopen(IabOtsuMaskNuc, params.maxObjSize);
% IabOtsuMaskNuc = imclose(IabOtsuMaskNuc,strel('disk',3));
IabOtsuMaskNuc = imdilate(IabOtsuMaskNuc,strel('disk',8));
IabOtsuMaskNuc = imfill(IabOtsuMaskNuc,'holes');

params.NucMask = IabOtsuMaskNuc ;

if params.plot_opt == 1
    figure(100)
    subplot_tight(2,2,1,[0.05 0.01])
    imshow(Iab_eq,[], 'InitialMagnification', 'fit');

    figure(100)
    subplot_tight(2,2,2,[0.05 0.01])
    % figure(61)
    imshow(IabOtsuMaskNuc,[], 'InitialMagnification', 'fit');

    figure(100)
    subplot_tight(2,2,3,[0.05 0.01])
    imshow(IabOtsuMaskNuc,[], 'InitialMagnification', 'fit');

    figure(100)
    subplot_tight(2,2,4,[0.05 0.01])
    imshow(IabOtsuMaskNuc,[], 'InitialMagnification', 'fit');
    tightfig();
    print(100 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ;
end
end

function IabMask = ad_thresh(Iab_eq, params)
global fig_num
global fig_path
tlt_name = params.ch_name ;
Iab_eq = mat2gray(Iab_eq) ;

%%
Iab_eq_no_zeros = Iab_eq ;
Iab_eq_no_zeros(Iab_eq_no_zeros==0)=Inf ;
step_size = min(Iab_eq_no_zeros(:)) ;
step_size = max(step_size, 1e-4) ;
Tab_v_diff = [] ;
for j=1:10
    Tab_vs = 2.^(j-1)*ones(1,10) ;
    Tab_v_diff = [Tab_v_diff, Tab_vs ] ;
end
Tab_v = cumsum(Tab_v_diff) ;
Tab_v = Tab_v * step_size ;
Tab_v(Tab_v>=1)= [] ;
Tab_v = [Tab_v, 1] ;
% step_size = max(0.002, step_size) ;
%%

% Tab_v = [0:step_size:1] ;
n_thresh = numel(Tab_v) ;
area_v = zeros(1, n_thresh) ;
n_obj_v = zeros(1, n_thresh) ;
for j = 1:n_thresh
    Tab = Tab_v(j) ;
    IabMask = im2bw(Iab_eq,Tab); %Corrected Otsu to reduce threshold
%     IabMask = im2bw(Iab_nobg_ad,Tab); %Corrected Otsu to reduce threshold
    abStats = regionprops(IabMask,'Area');
    abArea = cat(1, abStats.Area);
    area_v(j) = mean(abArea) ; 
    n_obj_v(j) = numel(abArea) ;
%     figure(8)
%     imshow(IabMask,[], 'InitialMagnification', 'fit');
%     title(num2str(Tab), 'fontsize', 10)
%     print(8 ,'-dpng','-r150', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ;
end
%%
Tab_v_diff = diff(Tab_v) ;
n_obj_v_dif = diff(n_obj_v)./Tab_v_diff;
area_v_dif = diff(area_v)./Tab_v_diff;
area_v_dif2 = diff(area_v_dif)./Tab_v_diff(1:end-1);
[n_obj_max, max_ind] = max(n_obj_v_dif) ;
[n_obj_min, min_ind] = min(n_obj_v_dif) ;
[area_max, area_max_ind] = max(area_v_dif2) ;
area_convg_ind = find(area_v < 30, 1, 'first') ;
% max_ind = find(n_obj_v_dif == n_obj_max, 1,'last') ;
% min_ind = find(n_obj_v_dif == n_obj_min, 1,'first') ;

%%
try
    %  [~, min_ind] = min(area_v_dif) ;
%  Tab_area = Tab_v(min_ind) ;    
%     Tab_area = Tab_v(area_max_ind+1) ;
%     Tab_area = Tab_v(area_convg_ind) ;
%     Tab = Tab_area  ;
    [pks,locs] = findpeaks(n_obj_v, Tab_v) ;
    locs(pks<100) = [] ; % remove small peaks
    Tab =  locs(end) ; % find the last peak    
    
catch
    Tab_max_obj = interp1(n_obj_v_dif(max_ind:min_ind) ,Tab_v(max_ind+1:min_ind+1),0,'linear');  
    Tab = Tab_max_obj ;
end
if isempty(Tab)
    pause(1)
end
%%
IabMask = im2bw(Iab_eq, Tab); 
% end 
%     figure(8)
%     imshow(IabMask,[], 'InitialMagnification', 'fit');
%     title(num2str(Tab), 'fontsize', 10)
%     print(8 ,'-dpng','-r150', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ;

%%
if params.plot_opt == 1
    figure(19)
    subplot(2,2,1)
    plot(Tab_v, area_v)
    xlabel('Threshold')
    ylabel('Area')
    set(gca, 'xscale', 'log')
    ylim = get(gca, 'ylim') ;
    format_fig2(1)
    h1 = title(tlt_name , 'fontsize', 6) ;  set(h1,'interpreter','none')
    line([Tab, Tab], ylim, 'LineStyle', '--', 'Color','k') 
    %
    % figure(20)
    subplot(2,2,2)
    plot(Tab_v, n_obj_v)
    xlabel('Threshold')
    ylabel('N(obj)')
    set(gca, 'xscale', 'log')
    ylim = get(gca, 'ylim') ;
    format_fig2(1)
    h1 = title(tlt_name , 'fontsize', 10) ;  set(h1,'interpreter','none')
    line([Tab, Tab], ylim, 'LineStyle', '--', 'Color','k') 

    %
    % figure(21)
    subplot(2,2,3)
    plot(Tab_v(2:end), area_v_dif)
    xlabel('Threshold')
    ylabel('Area')
    set(gca, 'xscale', 'log')
    ylim = get(gca, 'ylim') ;
    format_fig2(1)
    line([Tab, Tab], ylim, 'LineStyle', '--', 'Color','k') 

    % figure(22)
    subplot(2,2,4)
    plot(Tab_v(2:end), n_obj_v_dif)
    xlabel('Threshold')
    ylabel('N(obj)')
    set(gca, 'xscale', 'log')
    ylim = get(gca, 'ylim') ;
    format_fig2(1)
    line([Tab, Tab], ylim, 'LineStyle', '--', 'Color','k') 
    tightfig();
    % figure(23)
    % plot(Tab_v(2:end-1), area_v_dif2)
    % xlabel('Threshold')
    % ylabel('Area')
    % set(gca, 'xscale', 'log')
    % ylim = get(gca, 'ylim') ;
    % format_fig2(2)
    % line([Tab, Tab], ylim, 'LineStyle', '--', 'Color','k') 

%     print(19 ,'-dpng','-r150', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ;
   
end
end

function [imgRoundNames, targets]= ParseImgName(imgRoundNames)
% (1)skip washout rounds, basline, duplicates
% (2)replace seqence names with protein names
    seq_target_cell = {'actin','p2';'Tuj-1','p3';'cortactin','p4';...
        'SHANK3','p6';'ARPC2','p7'; 'bassoon','p8'; 'synapsin-1(2)','p9-t2'; 'synapsin-1','p9';...
        'Homer-1b/c','p10';'NR2B','p12';'PSD-95','p1-'} ;
    targets = [] ;
    
    for j = 1:numel(imgRoundNames)
        if isempty(regexpi(imgRoundNames{j}, '(wo|baseline|x)', 'match'))
            targets = [targets, j] ;                    
        end                       
        for k =  1:size(seq_target_cell,1)            
           if ~isempty(regexpi(imgRoundNames{j},seq_target_cell{k,2}, 'match'))
                imgRoundNames{j} = seq_target_cell{k,1} ;
            end 
        end
    end
end     

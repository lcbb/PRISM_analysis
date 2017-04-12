function  confocal_colocalization(imgsProjAligned, params, imgRoundNames) 
global fig_num fig_path;
close all

clim = [0.2, 0.97];
uniqueFields = [1:numel(imgsProjAligned)] ;
% imgRoundNames = [imgRoundNames; {'DAPI'; 'MAP2'}] ;
ch_name = imgRoundNames ;
plot_ind  = [1:numel(imgRoundNames)] ;
        
maskfile = fullfile(params.outputImgsPath,params.Condi, 'mask.mat') ;
marker_ind_bi = cellfun(@(x) ~isempty(strfind(x,'synapsin')), imgRoundNames,'uniformoutput',0) ; 
marker_ind_bi = cat(1,marker_ind_bi{:}) ;
marker_ind = find(marker_ind_bi, 1, 'first') ;
ab_ind_all = [1:numel(plot_ind)] ;
ab_ind_all(marker_ind) = [] ;

tlt_name = params.inputImgsPath ;
pixsize = 0.1 ; % check this
    
%%
params.tophat = 1;
params.disksize = 100;
% fraction parameters for raw images 

params.otsuCorrectFrac = 0;
params.nucOtsuCorrectFrac = params.otsuCorrectFrac + 1.5;
%Maximum continuous region size of synapses (to remove nuclei based on
%size):
params.maxObjSize = round(15^2*pi);
params.marker_ind = marker_ind ;
params.ab_ind_all =  ab_ind_all ;
params.nucRm = 1 ; %remove the bead or not
for f = 1:numel(uniqueFields)
    %         for f = 2:2
    %%
    disp(['process field ', num2str(f)])
        sizeI = size(imgsProjAligned{f}{1,1});
        Idapi = imgsProjAligned{f}{2,1};
        IMAP2 = imgsProjAligned{f}{1,1};
        ItargetStack = cat(3,imgsProjAligned{f}{3,plot_ind});
    %     ItargetStack = cat(3,ItargetStack,imgsProjAligned{f}{[2,1],1});
        if params.tophat 
            for i = 1:size(ItargetStack,3)
                ItargetStack(:,:,i) = imtophat(ItargetStack(:,:,i),strel('disk',params.disksize));
            end
        end
    params = segmentNuc(Idapi, params) ;    
    Imarker = mat2gray(ItargetStack(:,:,marker_ind)) ;
    clim = stretchlim(Imarker, clim) ;     
    params.ch_name = ch_name(marker_ind);
    markerMask = segment_synapse(Imarker, params) ;        

    %%
    abMaskAll = [] ;

    for k = 1:numel(ab_ind_all)
        ab_ind = ab_ind_all(k) ;    
        Iab = mat2gray(ItargetStack(:,:,ab_ind)); 
        
    %     ab_params = segmentNuc(Iab_eq, ab_params) ;
    %     nucMask = ab_params.NucMask ;
    %     ab_params.NucMask = nucMask ;
        params.ch_name = ch_name(ab_ind);
        abMask = segment_synapse(Iab, params) ;
        abMaskAll = cat(3, abMaskAll, abMask) ;  
    end
     markerMaskCell{f} = markerMask ;
     abMaskCell{f} = abMaskAll ;     
end
save(maskfile, 'markerMaskCell', 'abMaskCell', 'params', 'plot_ind', 'ch_name', '-v7.3')   
end


function abMask = segment_synapse(Iab_eq, params)
%% Do Otsu thresholding on marker and ab.
global fig_num
global fig_path

if params.nucRm == 1 ;
    IabOtsuMaskNuc = params.NucMask ;        
    sat_lim = stretchlim(Iab_eq(~IabOtsuMaskNuc)) ;    
    Iab_ad = imadjust(Iab_eq, sat_lim);
%     Iab_ad = adapthisteq(Iab_eq); 
%     Iab_ad = Iab_eq ;
    Iab_ad(IabOtsuMaskNuc) =  0 ;
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
%%watershed doesn't work well for mask of small objects as synapses. Apply wartershed to the original image instead 
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

IabMask_perim = bwperim(abMask);
Iab_shed_perim = imoverlay(Iab_ad, IabMask_perim, [1 .3 .3]);
%%
if params.plot_opt == 1
    figure(5)
    subplot_tight(2,2,1,[0.05 0.01])
    imshow(Iab_ad,[], 'InitialMagnification', 'fit');
    figure(5)
    subplot_tight(2,2,2,[0.05 0.01])
    imshow(imadjust(Iab_tophat),[], 'InitialMagnification', 'fit');
    figure(7)
    imshow(IabMask,[], 'InitialMagnification', 'fit');
    h1 = title(params.ch_name , 'fontsize', 6) ;  set(h1,'interpreter','none')
    
        figure(9)
    % subplot(2,2,3)
    imshow(Iab_ad_perim,[], 'InitialMagnification', 'fit');

    figure(5)
    subplot_tight(2,2,3,[0.05 0.01])
    imshow(abMask,[], 'InitialMagnification', 'fit'); h1 = title(params.ch_name , 'fontsize', 6) ;  set(h1,'interpreter','none')
    
    figure(5)
    subplot_tight(2,2,4,[0.05 0.01])
    imshow(Iab_shed_perim,[], 'InitialMagnification', 'fit');
    tightfig();
%     print(5 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ; save('.\startup.mat', 'fig_num');

    figure(11)
    imshow(Iab_overlay_max,[], 'InitialMagnification', 'fit');

    figure(12)
    imshow(Iab_com,[], 'InitialMagnification', 'fit');

    

    print(5 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ; save('.\startup.mat', 'fig_num');
    % print(13 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ; save('.\startup.mat', 'fig_num');
    % print(131 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ; save('.\startup.mat', 'fig_num');
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
    subplot(2,2,1)
    imshow(Iab_eq,[], 'InitialMagnification', 'fit');

    figure(100)
    subplot(2,2,2)
    % figure(61)
    imshow(IabOtsuMaskNuc,[], 'InitialMagnification', 'fit');

    figure(100)
    subplot(2,2,3)
    imshow(IabOtsuMaskNuc,[], 'InitialMagnification', 'fit');

    figure(100)
    subplot(2,2,4)
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
% Iab_ad = imadjust(Iab_eq);
% Iab_close = imclose(Iab_eq,strel('disk',5));
% Iab_close = imclose(Iab_ad,strel('disk',5));
% Iab_nobg = mat2gray(imopen(Iab_close,strel('disk',5))) ;
% Iab_close = mat2gray(Iab_close) ;
% Iab_close_ad = imadjust(Iab_close);
% Iab_nobg_ad = imadjust(Iab_nobg);

% figure(6)
% imshow(Iab_ad,[], 'InitialMagnification', 'fit');
% figure(7)
% imshow(Iab_nobg_ad,[], 'InitialMagnification', 'fit');
% figure(8)
% imshow(Iab_close_ad,[], 'InitialMagnification', 'fit');
% print(6 ,'-dpng','-r150', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ;
% print(7 ,'-dpng','-r150', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ;
% print(8 ,'-dpng','-r150', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ;


%%
% n_thresh = 40 ;
% 
% area_v = zeros(1, n_thresh) ;
% n_obj_v = zeros(1, n_thresh) ;
% % Tab_v = [1:n_thresh]./n_thresh ;
% Tab_v = 2.^[-n_thresh+1:0] ;
% for j = 1:n_thresh
%     Tab = Tab_v(j) ;
%     IabMask = im2bw(Iab_eq,Tab); %Corrected Otsu to reduce threshold
% %     IabMask = im2bw(Iab_nobg_ad,Tab); %Corrected Otsu to reduce threshold
%     abStats = regionprops(IabMask,'Area');
%     abArea = cat(1, abStats.Area);
%     area_v(j) = mean(abArea) ; 
%     n_obj_v(j) = numel(abArea) ;
% %     figure(8)
% %     imshow(IabMask,[], 'InitialMagnification', 'fit');
% %     title(num2str(Tab), 'fontsize', 10)
% %     print(8 ,'-dpng','-r150', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ;
% end
% %%
% n_obj_v_dif = diff(n_obj_v);
% area_v_dif = diff(area_v);
% area_v_dif2 = diff(area_v_dif);
% [n_obj_max, max_ind] = max(n_obj_v_dif) ;
% [n_obj_min, min_ind] = min(n_obj_v_dif) ;
% [area_max, area_max_ind] = max(area_v_dif2) ;
% % max_ind = find(n_obj_v_dif == n_obj_max, 1,'last') ;
% % min_ind = find(n_obj_v_dif == n_obj_min, 1,'first') ;
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
Tab_v(Tab_v>1)= [] ;
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
% IabMask = im2bw(Iab_eq,Tab_max_obj); %Corrected Otsu to reduce threshold
% if sum(n_obj_v_dif == n_obj_max) ==1
%     IabMask = im2bw(Iab_eq, max(Tab_area, Tab_max_obj)); %Corrected Otsu to reduce threshold
IabMask = im2bw(Iab_eq, Tab); %Corrected Otsu to reduce threshold
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


    print(19 ,'-dpng','-r150', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ;
    % print(20 ,'-dpng','-r150', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ;
    % print(21 ,'-dpng','-r150', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ;
    % print(22 ,'-dpng','-r150', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ;
    % print(23 ,'-dpng','-r150', [fig_path,num2str(fig_num, '%03d'),'.png']) ; fig_num = fig_num +1 ;
end
end
function [Io, clim] = imstack2RGB(im, sat_frc, opt, cmap)
if nargin <3
%     opt = {'fraction', 'fraction', 'fraction'};
    opt = repmat({'fraction'}, 1, size(im,3));
    cmap = 'gray' ;
elseif nargin <4
    cmap = 'gray' ;
end 

cmin = zeros(1, size(im,3)) ;
cmax = zeros(1, size(im,3)) ;
if size(sat_frc,1)== 1 ;
    sat_frc = repmat(sat_frc,[size(im,3), 1]) ;
end

for j = 1:size(im,3)      
   switch opt{1}
       case 'fraction'
        Low_High = stretchlim(im(:,:,j), sat_frc(j,:)) ;    
        cmin(j) = Low_High(1);
        cmax(j) = Low_High(2);
        
       case 'limit'
        cmin(j) = sat_frc(j,1) ;
        cmax(j) = sat_frc(j,2) ;
   end
end
clim = {cmin, cmax} ;

% I = cell(size(im,3), 1) ;
% for j = 1:length(I)   
%    I{j} = im(:,:,j) ;   
% end
% clear im
Io = STORMcell2img(im,'cmin',cmin,'cmax',cmax, 'colormap', cmap);
end
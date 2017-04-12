function triangleHeatMap(pccImarkerAbSyn, pccImarkerAbSynCi, tlt, xlabel, ylabel, cmap, tri_opt)
global fig_num fig_path work_path;
pccImarkerAbSynCell = mat2cell(pccImarkerAbSyn, ones(1,size(pccImarkerAbSyn,1)),  ones(1,size(pccImarkerAbSyn,2))) ;

pccImarkerAbSynCell = cellfun(@(x)num2str(x,'%0.2f'),pccImarkerAbSynCell,'UniformOutput', false);

% pccImarkerAbSynCiCell = cellfun(@(c2)sprintf(['¡À', c2]), pccImarkerAbSynCiCell,'UniformOutput', false);
% pccImarkerAbSynCiCell = cellfun(@(c1,c2)sprintf([c1, '\n', c2]), pccImarkerAbSynCell, pccImarkerAbSynCiCell,'UniformOutput', false);
% pccImarkerAbSynCiCell = cellfun(@(c1,c2)sprintf([c1, c2]), pccImarkerAbSynCell, pccImarkerAbSynCiCell,'UniformOutput', false);
if ~isempty(pccImarkerAbSynCi)
   pccImarkerAbSynCiCell = mat2cell(pccImarkerAbSynCi, ones(1,size(pccImarkerAbSynCi,1)),  ones(1,size(pccImarkerAbSynCi,2))) ;
   pccImarkerAbSynCiCell = cellfun(@(x)num2str(x,'%0.2f'),pccImarkerAbSynCiCell,'UniformOutput', false);
   pccImarkerAbSynCiCell = cellfun(@(c1,c2)sprintf([c1, '\n', '±', c2]), pccImarkerAbSynCell, pccImarkerAbSynCiCell,'UniformOutput', false);
else
   pccImarkerAbSynCiCell = pccImarkerAbSynCell ; 
end
   

% Deleted the upper triangle of the matrix
pccImarkerAbSynTri = pccImarkerAbSyn ;
if tri_opt == 1
    nanL = tril(NaN(size(pccImarkerAbSyn,1), size(pccImarkerAbSyn,1)), -1)' ;
    pccImarkerAbSynTri = pccImarkerAbSynTri + nanL ;
    for j = 1:size(pccImarkerAbSyn,1)
        for k = j+1:size(pccImarkerAbSyn,2)
            pccImarkerAbSynCiCell{j,k} = ' ' ;
        end
    end
    %
    pccImarkerAbSynCiTri = pccImarkerAbSynCi ;
    % nanL = tril(NaN(numel(ab_ind_all), numel(ab_ind_all)), -1)' ;
    pccImarkerAbSynCiTri = pccImarkerAbSynCiTri + nanL ;
else
end
%%
figure(1)
% heatmap(pccImarkerAbSynTri, xlabel, ylabel,pccImarkerAbSynCiCell,...
% 'TickAngle', 45,'TickFontsize',10,'Fontsize', 12,'Colormap', cmap,...
% 'NaNColor', [1 1 1], 'Textcolor', [0 0 0],'MinColorValue', -1, 'MaxColorValue',1, 'TickTexInterpreter',1) ;% 'RowLabelsRotate', 45, 'ColumnLabelsRotate', 45) ;
heatmap(pccImarkerAbSynTri, xlabel, ylabel,pccImarkerAbSynCiCell,...
'TickAngle', 45,'TickFontsize',12,'Fontsize', 10,'Colormap', cmap,...
'NaNColor', [1 1 1], 'Textcolor', [0 0 0],'MinColorValue', -max(pccImarkerAbSynTri(:)),...
'MaxColorValue',max(pccImarkerAbSynTri(:)), 'TickTexInterpreter',1,'ShowAllTicks', 1) ;% 'RowLabelsRotate', 45, 'ColumnLabelsRotate', 45) ;
% caxis([-1 1]) ;
% caxis([0 max(pccImarkerAbSynTri(:))]) ;
% caxis auto
colormap(cmap)
axis image
set(gca, 'box', 'off')
colorbar
title(tlt,'Fontsize', 15)
% tightfig();

print(1 ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png']) ; 
print(1 ,'-dpdf','-r300', [fig_path,num2str(fig_num, '%03d'),'.pdf']) ; fig_num = fig_num +1 ; save([work_path, 'startup.mat'], 'fig_num');     
end
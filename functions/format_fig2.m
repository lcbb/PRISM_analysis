function format_fig2(varargin)
% figure export settings (default is opt=2)
% specify the font size (), and the figure size (cm)
figsizeX = 20 ;
figsizeY = 15 ;
fig_color = [1 1 1] ;

if isempty(varargin) 
    opt = 1;
else
    opt=varargin{1};
end
% if ~any(opt==[1 2 3 4 5 6 7 8])
%     opt = 2;
% end  

if opt == 1;
        fontsize = 16;
        fontsize_leg = 12;
    elseif opt == 2
        fontsize = 24;
        fontsize_leg = 22;
    elseif opt == 3
        fontsize = 30;
        fontsize_leg = 26;
    elseif opt == 4
        fontsize = 24;
        fontsize_leg = 22;
        figsizeX = 20 ;
        figsizeY = 20 ;
    elseif opt == 5
        fontsize = 24;
        fontsize_leg = 22;
        figsizeX = 20 ;
        figsizeY = 15 ;  
    elseif opt == 6
        fontsize = 24;
        fontsize_leg = 22;
        figsizeX = 25 ;
        figsizeY = 15 ;
    elseif opt == 7
        fontsize = 24;
        fontsize_leg = 22;
        figsizeX = 25 ;
        figsizeY = 20 ;
    elseif opt == 8
        fontsize = 24;
        fontsize_leg = 22;
        figsizeX = 20 ;
        figsizeY = 25 ;
    elseif opt == 9
        fontsize = 20;
        fontsize_leg = 18;
        figsizeX = 60 ;
        figsizeY = 30 ;
    elseif opt == 10
        fontsize = 16;
        fontsize_leg = 14;
        figsizeX = 50 ;
        figsizeY = 15 ;
    elseif opt == 11
        fontsize = 24;
        fontsize_leg = 22;
        figsizeX = 9 ;
        figsizeY = 20 ;
    elseif opt == 12
        fontsize = 24;
        fontsize_leg = 22;
        figsizeX = 20 ;
        figsizeY = 20 ;
        fig_color = [0 0 0] ;
    elseif opt == 13
        fontsize = 24;
        fontsize_leg = 22;
        figsizeX = 12 ;
        figsizeY = 16 ;
        
end


pos = [5 10 figsizeX figsizeY];% figure position and size
% export paper properties
set(gcf,'PaperType','<custom>','unit','centimeters')
set(gcf,'PaperUnits','centimeters')
set(gcf,'Position',pos)
set(gcf,'PaperSize',[figsizeX figsizeY],...
    'PaperPosition', [0 0 figsizeX figsizeY], 'InvertHardcopy','off',...
    'Color',fig_color);


% legend font size
legend_handle =legend();
set(legend_handle,'FontSize', fontsize_leg)% legend font size

% axis label size
set(get(gca,'XLabel'),'FontSize',fontsize)
set(get(gca,'YLabel'),'FontSize',fontsize)

% axis font size
set(gca,'FontSize', fontsize)
set( gca                       , ...
    'FontName'   , 'Arial' );





% xSize = 8; ySize = 12;
% xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
% set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
% set(gcf,'Position',pos)
% set(gcf,'unit','centimeters','position',pos) % figure size
% set(gcf,'PaperPositionMode','auto')
% 
% % strech axis to fill window
% TI = get(gca,'TightInset');  OP = get(gca,'OuterPosition');  IP = get(gca,'Position');
% Pos = OP + [[TI(1) 0.5*TI(2)], ([-6*TI(1) -0*TI(2)] -TI(3:4)) ];  
% Pos = [TI(1:2) IP(3:4)] ; 
% set( gca,'Position',Pos);
% pos = [10 5 TI(3:4)];% figure position and size
% pos = [10 5 figsizeX figsizeY];% figure position and size
% 
% % set legend position to 'Best'
% set(legend_handle,'Location','Best')% legend font size
box on
end


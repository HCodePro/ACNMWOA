function save_pic(filename1)
%SAVE_PIC 此处显示有关此函数的摘要
%   此处显示详细说明


    a=findobj(gcf); % get the handles associated with the current figure
    allaxes=findall(a,'Type','axes');
    % alllines=findall(a,'Type','line');
    alltext=findall(a,'Type','text');
    set(allaxes,'FontName','Times','LineWidth',1,'FontSize',14,'FontWeight','bold');
    % set(alllines,'Linewidth',1);
    set(alltext,'FontName','Times','FontSize',14,'FontWeight','bold')
    %
    set(gcf, 'PaperUnits', 'inches');
    % set(gcf, 'PaperUnits', 'centimeters');
    krare=3.5;
    % x_width=krare*1.618 ;
    x_width=krare*5/3 ;
    %  x_width=3*1;
    y_width=krare*4/3;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
    % set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 4])
    set(gca, ...
        'Box'         , 'on'     , ...
        'TickDir'     , 'in'     , ...
        'TickLength'  , [.02 .02] , ...
        'XMinorTick'  , 'on'      , ...
        'YMinorTick'  , 'on'      , ...
        'YGrid'       , 'off'      , ...
        'XGrid'       , 'off'      , ...
        'XColor'      , [.3 .3 .3], ...
        'YColor'      , [.3 .3 .3], ...
        'LineWidth'   , 1         );
    axis tight
%     grid on
    box on
    
    
%     filename1=  [[pathName '\'],[currentFilename(1:end-4) '_hist']];
%     saveas(gcf,[filename,'-',Function_name,'-carve'],'fig')
    saveas(gcf, filename1,'fig')
%     [[pathName '\'],[currentFilename(1:end-4) '_hist' '.png']]

    print(filename1,'-dtiff', '-r300'); %<-Save as PNG with 300 DPI



end


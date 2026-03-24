function save_pic_seg(image,filename,gBest,level)
    set (gcf,'position',[500,500,800,200] );
    plot(imhist(image),'LineWidth', 1);
    set(gca,'XTick',0:100:300)
    set(gca,'YTick',0:round(max(imhist(image))/3):max(imhist(image)))
    axis([0,300,0,max(imhist(image))]);
    hold on;
    for c=1:level
        plot([gBest(level+c) gBest(level+c)], get(gca, 'YLim'), '-r', 'LineWidth', 1);
    end
    hold off;
%     filename=[[pathName '\'],[currentFilename(1:end-sign) '_seg']];             
    a=findobj(gcf); % get the handles associated with the current figure
    allaxes=findall(a,'Type','axes');
    % alllines=findall(a,'Type','line');
    alltext=findall(a,'Type','text');
    set(allaxes,'FontName','Times','LineWidth',1,'FontSize',14,'FontWeight','bold');
    % set(alllines,'Linewidth',1);
    set(alltext,'FontName','Times','FontSize',14,'FontWeight','bold')
    %
    krare=3.5;
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
    saveas(gcf, filename,'fig')
    print(filename,'-dtiff', '-r300'); %<-Save as PNG with 300 DPI
end


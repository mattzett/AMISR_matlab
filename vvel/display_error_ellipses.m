% DISPLAY ERROR ELLIPSES

figure(100);

plot_idx = 1;
for y = Ny:-1:1,
    for x = 1:Nx,
        
        subplot(Ny,Nx,plot_idx),
        
        error_ellipse(squeeze(Pvest(y,x,1:2,1:2)),'conf',0.394);
        xlim([-200,200]);
        ylim([-150,150]);
        
        set(gca,'XTick',[-200,-100,0,100,200]);
        set(gca,'XTickLabel',{});
        
        set(gca,'YTick',[-150,-75,0,75,150]);
        set(gca,'YTickLabel',{});
        
        axis square;
        
        % set(gca,'XTick',[],'YTick',[]);
        
        if y == 1,
            xlabel('v_{pe} (m/s)');
            set(gca,'XTickLabel',[-200,-100,0,100,200]);
        end;
        
        if x == 1,
            ylabel('v_{pn} (m/s)');
            set(gca,'YTickLabel',[-150,-75,0,75,150]);
        end;
        
        plot_idx = plot_idx + 1;
    end;
end;

set(gcf,'PaperUnits','inches');
papersize = [5, 8];
set(gcf,'PaperSize',papersize);
set(gcf,'PaperPosition',[0 0 papersize]);
outfile = fullfile(outputdir,'error_ellipses');
print('-f100','-depsc2','-r300',outfile);

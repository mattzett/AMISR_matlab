% DISPLAY ERROR ELLIPSES

plot_idx = 1;
for y = Ny:-1:1,
    for x = 1:Nx,
        
        subplot(Ny,Nx,plot_idx),
        
        error_ellipse(squeeze(Pvest(y,x,1:3,1:3)),'conf',0.394);
        xlim([-200,200]);
        ylim([-150,150]);
        
        set(gca,'XTick',[-200,-100,0,100,200]);
        set(gca,'XTickLabel',{});
        
        set(gca,'YTick',[-150,-75,0,75,150]);
        set(gca,'YTickLabel',{});
        
        set(gca,'ZTick',[-150,-50,50,150]);
        set(gca,'ZTickLabel',{});
        
        axis equal;
        
        grid on;
        view(36,10);
        
        % set(gca,'XTick',[],'YTick',[]);
        
        if y == 1,
            xlabel('v_{pe} (m/s)');
            ylabel('v_{pn} (m/s)');
            set(gca,'XTickLabel',[-200,-100,0,100,200]);
            set(gca,'YTickLabel',[-150,-75,0,75,150]);
        end;
        
        if x == 1,
            zlabel('v_{ap} (m/s)');
            set(gca,'ZTickLabel',[-150,-50,50,150]);
        end;
        
        plot_idx = plot_idx + 1;
    end;
end;
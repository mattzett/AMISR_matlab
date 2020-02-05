% display_Ti

pos_im = get(gca,'Position');

a = 15;
contour(xp/1e3, yp/1e3,...
        filter2(1/a^2*ones(a),Ti_I,'same'),...
        logspace(log10(cmin),log10(cmax),10),...
        'LineWidth',1);
caxis([cmin,cmax]);
colormap(hot(20));
cbar_handle = colorbar;
set(get(cbar_handle,'ylabel'),'String','Temperature (K)','FontSize',6);

% Resize + reposition colorbar
pos_cbar = get(cbar_handle,'Position');
pos_cbar(3) = 0.85*pos_cbar(3);                                 % "width"
pos_cbar(1:2) = pos_im(1:2) + [pos_im(3),0] + [pos_cbar(3),0];  % "x" & "y"

set(cbar_handle,'Position',pos_cbar);

% set(cbar_handle,'OuterPosition',[1 1 sqrt(1/2) 1].*get(cbar_handle,'OuterPosition'));

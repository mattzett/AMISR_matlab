% display_dvlos_overlay

dotsize = min(3 + dvlos(:)/10,...
              30);

Ncolors = 64;
vlos_rounded = linspace(-1000,1000,Ncolors);
vlos_color = zeros(size(vlos(:)));
for pos_idx = 1:length(vlos(:)),
    vlos_color(pos_idx) = vlos_rounded( abs(vlos_rounded - vlos(pos_idx)) == min( abs(vlos_rounded - vlos(pos_idx)) ) );
end;

scatter(xr(:)/1e3,yr(:)/1e3,dotsize,vlos_color,'Linewidth',1);

h_cbar = colorbar;
set(h_cbar,'CLim',[-1000,1000],'YLim',[-1000,1000]);

colormap(bluewhitered(Ncolors));


% % This second version did plots dvlos as size and vlos as color, but the
% % colormap is weird.  Turns out MATLAB has a built-in function for scatter
% % plots, and that is the uncommented version above.
% 
% [dvlos_sorted,dvlos_order] = sort(dvlos);
% 
% cmap_vlos = jet(2*numel(vlos));
% vlos_interp = linspace(-500,500,2*numel(vlos)); 
% 
% hold on;
% 
% for pos_idx = 1:numel(xr),
%     
%   this_color = cmap_vlos( abs(vlos(pos_idx) - vlos_interp) == min(abs(vlos(pos_idx) - vlos_interp)), : );
%   this_size  = min( 3 + dvlos(pos_idx)/200, ...
%                     30 );
%   
%   plot(xr(pos_idx)/1e3,yr(pos_idx)/1e3,'o',...
%        'MarkerEdgeColor','none',...
%        'MarkerFaceColor',this_color,... % 'MarkerFaceColor',cmap_dvlos(dvlos_idx,:),...
%        'MarkerSize', this_size );
%   
% %   fprintf('%3d, (%f,%f,%f) ... %f\n',pos_idx,cmap_dvlos(dvlos_idx,:),dvlos(dvlos_idx));
% end;
% 
% % set(gca,'Colormap',cmap_dvlos);
% % set(gca,'CLim',[10, 1000]);
% % colorbar



% % This original version encoded dvlos as color, and |vlos| as marker size
% 
% [dvlos_sorted,dvlos_order] = sort(dvlos);
% 
% cmap_dvlos = jet(2*numel(dvlos));
% dvlos_interp = linspace(10,1000,2*numel(dvlos)); 
% 
% hold on;
% 
% for dvlos_idx = 1:numel(xr),         % Though we plot (xr,yr) in order of increasing dvlos, ...
%   pos_idx = dvlos_order(dvlos_idx);  % ... xr and yr are in a different order, so we have to index into them appropriately.
%   
%   this_color = cmap_dvlos(abs(dvlos(dvlos_idx) - dvlos_interp) == min(abs(dvlos(dvlos_idx) - dvlos_interp)),:);
%   this_size = 3 + abs(vlos(dvlos_idx))/150;
%   plot(xr(pos_idx)/1e3,yr(pos_idx)/1e3,'o',...
%        'MarkerEdgeColor','none',...
%        'MarkerFaceColor',this_color,... % 'MarkerFaceColor',cmap_dvlos(dvlos_idx,:),...
%        'MarkerSize', this_size );
%   
% %   fprintf('%3d, (%f,%f,%f) ... %f\n',pos_idx,cmap_dvlos(dvlos_idx,:),dvlos(dvlos_idx));
% end;
% 
% % set(gca,'Colormap',cmap_dvlos);
% % set(gca,'CLim',[10, 1000]);
% % colorbar

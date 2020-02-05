% Arrange los velocities into x and y axes for display
vlos_display = zeros(Ny,Nx);
for k = 1:Nr,
    vlos_display(PART==k) = vlos(k);
end;

figure(51),
  imagesc(x/1e3,y/1e3,vlos_display);
  colormap(jet(64));
  h=colorbar; set(get(h,'ylabel'),'String','Line-of-sight velocity (m/s)');
  axis xy;
  hold on;
  plot(0,0,'ko','MarkerSize',16,'LineWidth',3,'MarkerFaceColor','g');
  plot(xr/1e3,yr/1e3,'rx','MarkerSize',14,'LineWidth',3);
  quiver(X/1e3,Y/1e3,qscaling*v(:,:,1),qscaling*v(:,:,2),0,'c');
  hold off;
  
  xlabel('Ground distance east (km)');
  ylabel('Ground distance north (km)');
  title('Line-of-sight measurements');
  
% if output,
% 
%     FILENAME  = sprintf('%s/truth_%dx%d_Nr%d.png',outputdir,Nx,Ny,Nr);
% 
%     if ~exist(outputdir,'file'),
%         mkdir(outputdir);
%     end;
% 
%     print('-f51','-dpng','-r90',FILENAME);
% 
% end;


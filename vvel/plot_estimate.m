figure(52),

% Show the phantom vector velocity field
quiver(X/1e3,             Y/1e3,...
       v(:,:,1)*qscaling, v(:,:,2)*qscaling,...
       0,'c');
hold on;
quiver(-120,-50,1e3*qscaling,0,0,'c');
text(-130,-45,'1 km/s');

% Show the beam positions & radar position (at origin)
% plot(0,0,'ko','MarkerSize',10,'LineWidth',2,'MarkerFaceColor','g');
plot(xr/1e3,yr/1e3,'rx','MarkerSize',8,'LineWidth',2);

% Show the shear boundary
    plot(X(1,:)/1e3,vy/vx*(X(1,:)-xoffset)/1e3+yoffset/1e3,'k','LineWidth',2);
    
% Show the reconstruction pixel boundaries
for i = 1:size(VERTICES,1),
    xb = [xr(VERTICES(i,:)), xr(VERTICES(i,1))];
    yb = [yr(VERTICES(i,:)), yr(VERTICES(i,1))];
    zb = [zr(VERTICES(i,:)), zr(VERTICES(i,1))];
    plot3(xb/1e3,yb/1e3,zb/1e3,'k:')
end;

% Show the reconstructed velocity field (+ a reference vector)
quiver([X2/1e3;    -95],          [Y2/1e3; -50],...
       [vest(:,1);1000]*qscaling, [vest(:,2);0]*qscaling,...
       0,'b');
text(-95,-45,'1 km/s');

hold off;

axis equal;
axis([X(1),X(end),Y(1),Y(end)]/1e3);
xlabel('Ground distance east (km)');
ylabel('Ground distance north (km)');
title('Cyan: truth, Blue: estimate, Red: radar beam centers')

if output,
    if Bayesian,
        outputdir1 = sprintf('%s/%s/n%d_bayesian',outputdir,reconstruction_type,Nstd);
    else
        outputdir1 = sprintf('%s/%s/n%d',outputdir,reconstruction_type,Nstd);
    end;
    FILENAME = sprintf('%s/estimate_%03d.png',outputdir1,t);

    if ~exist(outputdir1,'file'),
        mkdir(outputdir1);
    end;

    print('-f52','-dpng','-r90',FILENAME);

end;


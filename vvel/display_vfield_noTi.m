% display_vfield_noTi
% For displaying vfield over optical without Ti.
% 
% % Compute the coordinates just at HEIGHT (where Ti was sampled)
% xr = rr.*cos(el) .* sin(az);
% yr = rr.*cos(el) .* cos(az);
% zr = rr.*sin(el);
% 
% % Rotated coordinates
% xyzgmag_atHEIGHT = [xr' yr' zr'] * Rgmag';
% xgmag_atHEIGHT = xyzgmag_atHEIGHT(:,1);
% ygmag_atHEIGHT = xyzgmag_atHEIGHT(:,2);
% zgmag_atHEIGHT = xyzgmag_atHEIGHT(:,3);
% 
% % 
% x3 = linspace(min(xgmag_atHEIGHT), max(xgmag_atHEIGHT), Nx+2);
% y3 = linspace(min(ygmag_atHEIGHT), max(ygmag_atHEIGHT), Ny+2);
% 
% xc = x3(2:end-1)';
% yc = y3(2:end-1)';

[XC,YC] = meshgrid(x2,y2);
temp = [XC(:), YC(:), HEIGHT*ones(numel(XC),1)] * Rgmag;
xc = temp(:,1);
yc = temp(:,2);
vx = vest(:,:,1);
vy = vest(:,:,2);

% Show the reconstructed velocity field
quiver(xc/1e3,     yc/1e3, ...
       vx(:)*qfac, vy(:)*qfac, ...
       0,...
       'Color','w',...
       'LineWidth',0.5); %,...
%       'MaxHeadSize',0.1);

% Quiver legend
quiver(100, -40, 1e3*qfac, 0, 0, 'Color', 'w', 'LineWidth', 0.5);
text(100,-20,'1 km/s','color','w','FontSize',6)
% quiver(-140, -40, 1e3*qfac, 0, 0, 'Color', 'w', 'LineWidth', 0.5);
% text(-140,-20,'1 km/s','color','w','FontSize',6)


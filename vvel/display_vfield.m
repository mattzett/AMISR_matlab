% display_vfield

% Compute the coordinates just at HEIGHT (where Ti was sampled)
xr = rr.*cos(el) .* sin(az);
yr = rr.*cos(el) .* cos(az);
zr = rr.*sin(el);

% Rotated coordinates
xyzgmag_atHEIGHT = [xr' yr' zr'] * Rgmag';
xgmag_atHEIGHT = xyzgmag_atHEIGHT(:,1);
ygmag_atHEIGHT = xyzgmag_atHEIGHT(:,2);
zgmag_atHEIGHT = xyzgmag_atHEIGHT(:,3);

% 
x3 = linspace(min(xgmag_atHEIGHT), max(xgmag_atHEIGHT), Nx+2);
y3 = linspace(min(ygmag_atHEIGHT), max(ygmag_atHEIGHT), Ny+2);

xc = x3(2:end-1)';
yc = y3(2:end-1)';
[XC,YC] = meshgrid(xc,yc);
temp = [XC(:), YC(:), HEIGHT*ones(numel(XC),1)] * Rgmag;
xc = temp(:,1);
yc = temp(:,2);
vx = vest(:,:,1);
vy = vest(:,:,2);

% Show the reconstructed velocity field
quiver(xc/1e3,     yc/1e3, ...
       vx(:)*qfac, vy(:)*qfac, ...
       0,...
       'Color','c',...
       'LineWidth',2,...
       'MaxHeadSize',0.1);

% Quiver legend
quiver(-90, -40, 1e3*qfac, 0, 0, 'Color', 'c', 'LineWidth', 2);
text(-90,-30,'1 km/s','color','c')


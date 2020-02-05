% scatterplot_vi_Ti
%
% After running test04_bayesian, plot a scatter plot of vi versus Ti at
% zenith.

test04_bayesian

% Find the time index corresponding to 0900 - 1000 UT.
t1 = find( mtime(1,:) >= datenum(2008,03,26,11,30,00),1,'first' );
t2 = find( mtime(2,:) <= datenum(2008,03,26,12,30,00),1,'last' );

% What to do next:
% Find all the beams that lie within the pixels closest to the zenith at
% <HEIGHT> km.  Assign the corresponding Ti values to column 1 of an array,
% with column 2 being the corresponding |v|.
% Repeat for all t = t1..t2.

kx = cos(el) .* sin(az);
ky = cos(el) .* cos(az);
kz = sin(el);
% direction_vectors = [kx ky kz];

xr = rr.*kx;
yr = rr.*ky;
zr = rr.*kz;

Rgmag = [cos(dec),         -sin(dec),          0;
         sin(dip)*sin(dec), cos(dec)*sin(dip), cos(dip);
        -cos(dip)*sin(dec),-cos(dec)*cos(dip), sin(dip)];
direction_vectors = [kx' ky' kz'] * Rgmag';
xyzgmag = [xr' yr' zr'] * Rgmag';
xgmag = xyzgmag(:,1);
ygmag = xyzgmag(:,2);
zgmag = xyzgmag(:,3);

x2 = linspace(min(xgmag), max(xgmag), Nx+2);
y2 = linspace(min(ygmag), max(ygmag), Ny+2);

scatter_values = [];
for t = t1:t2,
    % Grab data from the pixels nearest zenith.
    for x = 2:3, % Pixel indices in longitude
        for y = 2:3, % Pixel indices in latitude

            MASK = (xgmag >= x2(x) & xgmag <= x2(x+2) ...
                  & ygmag >= y2(y) & ygmag <= y2(y+2));
            
            this_v  = sqrt(VX(x,y,t)^2 + VY(x,y,t)^2 + VZ(x,y,t)^2);
            this_Ti = Ti_all(MASK,t);

            new_scat = [repmat(this_v,length(this_Ti),1), this_Ti(:)];
            scatter_values = [scatter_values; new_scat];
            
        end; % for y
    end; % for x
    
end; %for t

% scatter_values = [];
% for t = t1:t2,
%     % Pixel (2,1)
%     x = 2; y = 1;
%     Ti_idx = [13,25,8];
%     this_v = sqrt(VX(x,y,t)^2 + VY(x,y,t)^2 + VZ(x,y,t)^2);
%     this_Ti = Ti_all(Ti_idx,t);
%     
%     new_scat = [repmat(this_v,3,1), this_Ti(:)];
%     scatter_values = [scatter_values; new_scat];
%     
%     % Pixel (3,1)
%     x = 3; y = 1;
%     Ti_idx = [13,8,26];
%     
%     this_v = sqrt(VX(x,y,t)^2 + VY(x,y,t)^2 + VZ(x,y,t)^2);
%     this_Ti = Ti_all(Ti_idx,t);
%     
%     new_scat = [repmat(this_v,3,1), this_Ti(:)];
%     scatter_values = [scatter_values; new_scat];
%     
%     % Pixel (2,2)
%     x = 2; y = 2;
%     Ti_idx = [25,8];
%     this_v = sqrt(VX(x,y,t)^2 + VY(x,y,t)^2 + VZ(x,y,t)^2);
%     this_Ti = Ti_all(Ti_idx,t);
%     
%     new_scat = [repmat(this_v,2,1), this_Ti(:)];
%     scatter_values = [scatter_values; new_scat];    
%     
%     % Pixel (3,2)
%     x = 3; y = 2;
%     Ti_idx = [8,26];
%     this_v = sqrt(VX(x,y,t)^2 + VY(x,y,t)^2 + VZ(x,y,t)^2);
%     this_Ti = Ti_all(Ti_idx,t);
%     
%     new_scat = [repmat(this_v,2,1), this_Ti(:)];
%     scatter_values = [scatter_values; new_scat];
%     
% end; %t

figure(6);
plot(scatter_values(:,1), scatter_values(:,2),'bo');
xlabel('Ion Velocity (m/s)');
ylabel('Ion Temperature (K)');
axis([0,2500,0,2500]);

%%
figure(5);
plot(xgmag,ygmag,'bx')
vline(x2,'k--')
hline(y2,'k--')
text(xgmag,ygmag,reshape(sprintf('%02d',1:26),2,26)')
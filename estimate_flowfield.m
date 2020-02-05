function estimate_flowfield(datadir,filelab,itstart,itfin)

%This script computes the ISR flow field according to the method presented 
% in Semeter et al, 2010 and Butler et al, 2010:
%
% http://onlinelibrary.wiley.com/doi/10.1029/2009JA014931/full
% http://onlinelibrary.wiley.com/doi/10.1029/2010RS004364/full


%PROCESSED LONGPULSE DATA
load([datadir,'/data_mat/',filelab,'_rawdata.mat']); % This is redundant but useful when working the code in segments
load([datadir,'/data_mat/',filelab,'_fieldgrid.mat']);


%SIZE OF DATA SET
[lz,lb,lt]=size(isvi);


%LINE-OF-SIGHT DIRECTIONS FOR EACH BEAM
azrad=az*pi/180; elrad=el*pi/180;
dec=22*pi/180; dip=77.5*pi/180;


% Compute unit vectors for each beam.
el1 = repmat(elrad,size(Rz,1),1);            % Replicate az & el, making each the same size as Alt (or Range)
az1 = repmat(azrad,size(Rz,1),1);
%range_idx = find( Rz>150 & Rz<400 ); % Find just the points within these altitudes.
%el2 = el1(range_idx);                      % Extract only those az & el
%az2 = az1(range_idx);


% % Unit vectors of the beam at each sample point.
% kx = cos(el2) .* sin(az2);
% ky = cos(el2) .* cos(az2);
% kz = sin(el2);


% % Cartesian coordinates of each sample point.
% xr = Rx(range_idx);
% yr = Ry(range_idx);
% zr = Rz(range_idx);


% % Rotation matrix from geo -> gmag
% Rgmag = [cos(dec),         -sin(dec),          0;
%          sin(dip)*sin(dec), cos(dec)*sin(dip), cos(dip);
%         -cos(dip)*sin(dec),-cos(dec)*cos(dip), sin(dip)];

    
% % Rotate direction vectors from geo -> gmag
% direction_vectors = [kx ky kz] * Rgmag';


% % Rotate coordinates to geomagnetic
% xyzgmag = [xr,yr,zr] * Rgmag';
% xgmag = xyzgmag(:,1);
% ygmag = xyzgmag(:,2);
% zgmag = xyzgmag(:,3);


%CALL THOMAS'S ESTIMATION CODE
%Nx=4; Ny=4;   %butler default
%alpha=3;      %regularization parameter, butler default
%Nx=3; Ny=3;
%alpha=5;
%Nx=4; Ny=4;
Nx=7; Ny=7;
alpha=5;


%REMOVE NEGATIVE DENSITIES AND RELATED DATA
Neg = find(isne<1E-100);
isne(Neg)=NaN;
isti(Neg)=NaN;
isvi(Neg)=NaN;


% vest=zeros(Ny,Nx,3,lt);
vest=zeros(Ny,Nx,3,itfin-itstart);
vest_geog=vest;

for it=1:itfin-itstart
    %LET USER KNOW WHAT IS HAPPENING
    if mod(it,25)==0
       fprintf('VVEL_RASTER.M --> Processing scan %d of %d.  \n',itstart+it-1,itfin); 
    end
    
    
    %ESTIMATE VELOCITY
    vlos=isvi(:,:,it);
    dvlos=isdvi(:,:,it);
    
    
    %LOTS OF NANS IN SOME DATA SETS
    inds=find(isnan(vlos));
    vlos(inds)=0;    %for lack of better idea of how to keep this from crashing on the svd
    dvlos(inds)=1;   %ditto
    
    
    %PICK THE ALTITUDES TO USE IN FITS
    SNRnow=SNR(:,:,it);
    vthreshold=300;
    SNRthreshold=0.1;
    range_idx = find( Rz>150 & Rz<400 & dvlos<vthreshold & SNRnow>SNRthreshold);   %this filters on altitude as well as velocity error and signal to noise
    el2 = el1(range_idx);   % Extract only those az & el
    az2 = az1(range_idx);
    xr = Rx(range_idx);
    yr = Ry(range_idx);
    zr = Rz(range_idx);
    
    kx = cos(el2) .* sin(az2);
    ky = cos(el2) .* cos(az2);
    kz = sin(el2);
    

    % Rotation matrix from geo -> gmag
    Rgmag = [cos(dec),         -sin(dec),          0;
             sin(dip)*sin(dec), cos(dec)*sin(dip), cos(dip);
            -cos(dip)*sin(dec),-cos(dec)*cos(dip), sin(dip)];
    

    % Rotate direction vectors from geo -> gmag
    direction_vectors = [kx ky kz] * Rgmag';
    

    % Rotate coordinates to geomagnetic
    xyzgmag = [xr,yr,zr] * Rgmag';
    xgmag = xyzgmag(:,1);
    ygmag = xyzgmag(:,2);
    zgmag = xyzgmag(:,3);
    

    %APPLY THE BUTLER ALGORITHM
    [vest(:,:,:,it),~] = vfield_holistic_incompressible(vlos(range_idx),dvlos(range_idx),xgmag,ygmag,direction_vectors,Nx,Ny,alpha);
    
    
    %ROTATE VECTORS TO GEOGRAPHIC COORDS.
    vest_geog_tmp = zeros(Ny,Nx,3);
    for ix = 1:Nx
        for iy = 1:Ny
            vest_geog_tmp(iy,ix,:) = Rgmag'*squeeze(vest(iy,ix,:,it));
        end
    end
    vest_geog(:,:,:,it)=vest_geog_tmp;
end


%DEFINE A GRID FOR ESTIMATED VELOCITY
zref=300;    
[~,iz]=min(abs(Rz(:,1)-zref));
xyzgmag_atHEIGHT = [Rx(iz,:)' Ry(iz,:)' Rz(iz,:)'] * Rgmag';
xgmag_atHEIGHT = xyzgmag_atHEIGHT(:,1);
ygmag_atHEIGHT = xyzgmag_atHEIGHT(:,2);
xvm=linspace(min(xgmag_atHEIGHT),max(xgmag_atHEIGHT),Nx+2);
yvm=linspace(min(ygmag_atHEIGHT),max(ygmag_atHEIGHT),Ny+2);
xvm=xvm(2:end-1)';
yvm=yvm(2:end-1)';
[Xvm,Yvm] = meshgrid(xvm,yvm);
temp=[Xvm(:), Yvm(:), mean(Rz(iz,:))*ones(numel(Xvm),1)]*Rgmag;
Xvg=temp(:,1);
Yvg=temp(:,2);
Zvg=temp(:,3);


%SAVE PROCESSED DATA TO OUTPUT FILE
filelab2=datestr(exp_date(1,:),'ddmmmyyyy');
save([datadir,'/data_mat/',filelab2,'_',num2str(tint),'min','_vvelscans.mat'], ...
    'vest','vest_geog','Nx','Ny','range_idx','Xv*','Yv*','Zv*','xv*','yv*','Rgmag');



%MAKE A KEOGRAM BY AVERAGING THE TWO CENTER CUTS
Xvg=reshape(Xvg,[Ny Nx]);
Yvg=reshape(Yvg,[Ny Nx]);
i1=floor(Nx/2);
i2=i1+1;
Xvgcut=1/2*(Xvg(:,i1)+Xvg(:,i2));
Yvgcut=1/2*(Yvg(:,i1)+Yvg(:,i2));
Veastcuts=squeeze(1/2*(vest_geog(:,i1,1,:)+vest_geog(:,i2,1,:)));
Vnorthcuts=squeeze(1/2*(vest_geog(:,i1,2,:)+vest_geog(:,i2,2,:)));
t=datenum(exp_date);


figure; 

subplot(211);
imagesc(t(itstart:itfin),Yvgcut,Veastcuts);
axis xy;
datetick;
caxis([-3e3 3e3]);
xlabel('UT');
ylabel('dist N (km)');
c=colorbar;
ylabel(c,'eastward drift (m/s)');

subplot(212);
imagesc(t(itstart:itfin),Yvgcut,Vnorthcuts);
axis xy;
datetick;
caxis([-2e3 2e3]);
xlabel('UT');
ylabel('dist N (km)');
c=colorbar;
ylabel(c,'northward drift (m/s)');

print([datadir,'/plot_imgfiles/',filelab,'/drift_keo.png'],'-dpng','-r300') 

end

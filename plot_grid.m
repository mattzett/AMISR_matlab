function plot_grid(datadir,filelab,payloadflag)

%LOAD DATA
% addpath ./geom;
% load ./datasets/data_mat/18Jan2017_0.6min_rawdata.mat;

load([datadir,'/data_mat/',filelab,'_rawdata.mat']);


lp=length(az);

%MAGNETIC ZENITH
zenaz=205.7;
zenel=77.5;
[zentheta,zenphi]=azel2thetaphi(zenaz,zenel);

%CONVERT ANGLES
[radtheta,radphi]=azel2thetaphi(az,el);

%BEAM E-REGION INTERSECTION IN O
thetaB=radtheta;
phiB=radphi;
zB=isz(2,:); zB=zB(:);
rB=zB.*sec(thetaB);
[xB,yB,zBtmp]=spher2cart(rB,thetaB,phiB);

%BEAM E-REGION INTERSECTION IN O''
thetaBpp=zentheta;
phiBpp=zenphi;
zBpp=zB;
rBpp=zBpp.*sec(thetaBpp);
[xBpp,yBpp,zBpptmp]=spher2cart(rBpp,thetaBpp,phiBpp);

%FIELD LINE FOOTPOINT IN O'
xCp=0-xBpp;
yCp=0-yBpp;
zCp=0-zBpp;

%RADAR LOCATION IN O'
xAp=0-xB;
yAp=0-yB;
zAp=0-zB;

%FIELD LINE FOOTPOINT IN O
xC=xCp-xAp;
yC=yCp-yAp;
zC=zCp-zAp;

%SLOPE OF FIELD LINES
mx=xB-xC;
my=yB-yC;
mz=zB-zC;

%FIELD LINE GRID
zline=100:10:400;
lz=length(zline);
Bx=zeros(lz,lp); By=Bx; Bz=Bx;
for ip=1:lp
    t=(zline-zC(ip))/mz(ip);
    
    Bx(:,ip)=mx(ip)*t+xC(ip);
    By(:,ip)=my(ip)*t+yC(ip);
    Bz(:,ip)=mz(ip)*t+zC(ip);
end



%RADAR GRID
Rx=zeros(size(isz)); Ry=Rx; Rz=Rx;
for ip=1:lp
    range=isz(:,ip)*sec(thetaB(ip));
    
    Rx(:,ip)=range.*sin(thetaB(ip)).*cos(phiB(ip));
    Ry(:,ip)=range.*sin(thetaB(ip)).*sin(phiB(ip));
    Rz(:,ip)=range.*cos(thetaB(ip));
end



%ADD TRAJECTORY
if payloadflag==1   % ISINGLASS A, 2/22/2017, 10:14:00-10:24:52 UT (main payload), 10:14:00-10:23:55 UT (sub payload)
   
    %PFISR RADAR COORDS
    lat_radar = 65.1192;
    long_radar = -147.4303;

    %READ IN ISINGLASS TRAJECTORY - Main Payload
    trajdata=xlsread('./ISINGLASS/A/36303MainSubCombinedGPS.xlsx','A13:L572');  % grabs main payload trajectory data from t=0 through downleg till 90km, can go lower if needed
%     time=trajdata(:,1); % time since launch in s
    rlat=trajdata(:,8); % main payload geomagnetic (or geographic?) latitude
    rlon=trajdata(:,9); % main payload geomagnetic (or geographic?) longitude
    ralt=trajdata(:,10); % main payload altitude from ground in km
    
    [SRmainx,SRmainy,SRmainz]=Geodetic2ENU(rlat,rlon,ralt,lat_radar,long_radar);
    
    %READ IN ISINGLASS TRAJECTORY - Sub Payload
    trajdata2=xlsread('./ISINGLASS/A/36303MainSubCombinedGPS.xlsx','N13:Q572');  % grabs main payload trajectory data from t=0 through downleg till 90km, can go lower if needed
%     time2=trajdata2(:,1); % time since launch in s
    r2lat=trajdata2(:,3); % main payload geomagnetic (or geographic?) latitude
    r2lon=trajdata2(:,2); % main payload geomagnetic (or geographic?) longitude
    r2alt=trajdata2(:,4)/1e3; % main payload altitude from ground conterted to km
    
    [SRsubx,SRsuby,SRsubz]=Geodetic2ENU(r2lat,r2lon,r2alt,lat_radar,long_radar);

 
    
elseif payloadflag==2   % ISINGLASS B, 3/02/2017, ???-??? UT (main payload), ???-??? UT (sub payload)
   
    %PFISR RADAR COORDS
    lat_radar = 65.1192;
    long_radar = -147.4303;

    %READ IN ISINGLASS TRAJECTORY - Main Payload
    trajdata=xlsread('./ISINGLASS/B/36304_MainTrajectory.xlsx','A12:D12738');
%     time=trajdata(:,1); % time since launch in s
    rlat=trajdata(:,2); % main payload geomagnetic (or geographic?) latitude
    rlon=trajdata(:,3); % main payload geomagnetic (or geographic?) longitude
    ralt=trajdata(:,4); % main payload altitude from ground in km
    
    [SRmainx,SRmainy,SRmainz]=Geodetic2ENU(rlat,rlon,ralt,lat_radar,long_radar);
    
    %READ IN ISINGLASS TRAJECTORY - Sub Payload
    trajdata2=xlsread('./ISINGLASS/B/36304SubPayloadGPS.xls');  
%     time2=trajdata2(:,1); % time since launch in s
    r2lat=trajdata2(:,3); % main payload geomagnetic (or geographic?) latitude
    r2lon=trajdata2(:,2); % main payload geomagnetic (or geographic?) longitude
    r2alt=trajdata2(:,4)/1e3; % main payload altitude from ground conterted to km
    
    [SRsubx,SRsuby,SRsubz]=Geodetic2ENU(r2lat,r2lon,r2alt,lat_radar,long_radar);
    
    
    
elseif payloadflag==3  % VISIONS
    
    %PFISR RADAR COORDS
    lat_radar = 65.1192;
    long_radar = -147.4303;

    %READ IN VISIONS TRAJECTORY
    load('./VISIONS/posdat_ephem_alt_lat_long_t.mat')
    [SRx,SRy,SRz]=Geodetic2ENU(lat,lon,alt./1e3,lat_radar,long_radar);

end


%PLOTS
figure;
hold on;
for ip=1:lp
    mlines=plot3(Bx(:,ip),By(:,ip),Bz(:,ip),'k-');
end

for ip=1:lp
    rlines=plot3(Rx(:,ip),Ry(:,ip),Rz(:,ip),'b--');
    %    text(Rx(end,ip),Ry(end,ip),Rz(end,ip),num2str(ip))
    text(Rx(end,ip),Ry(end,ip),Rz(end,ip),{['az ',num2str(az(ip)),', el ',num2str(el(ip))],['BeamCode: ',num2str(beamcodes(ip))]})
end

pfisr=plot3(0,0,0,'bo','MarkerSize',15);

if payloadflag==1
    main=plot3(SRmainx,SRmainy,SRmainz,'go');
    sub=plot3(SRsubx,SRsuby,SRsubz,'r*');
    legend([mlines,rlines,main,sub,pfisr],{'Magnetic Field Lines','Radar Beams','Main Payload','Sub Payload','PFISR'})
    legend boxoff;
    
elseif payloadflag==2
    main=plot3(SRmainx,SRmainy,SRmainz,'go');
    sub=plot3(SRsubx,SRsuby,SRsubz,'r*');
    legend([mlines,rlines,main,sub,pfisr],{'Magnetic Field Lines','Radar Beams','Main Payload','Sub Payload','PFISR'})
    legend boxoff;
    
elseif payloadflag==3
    main=plot3(SRx,SRy,SRz,'go');
    legend([mlines,rlines,main,pfisr],{'Magnetic Field Lines','Radar Beams','VISIONS','PFISR'})
    legend boxoff;

end

axis tight;
grid on;
xlabel('dist. E of PFISR (km)');
ylabel('dist. N of PFISR (km)');
zlabel('altitude (km)')
% legend([mlines,rlines,main,sub,pfisr],{'Magnetic Field Lines','Radar Beams','Main Payload','Sub Payload','PFISR'})
% legend boxoff;

mkdir([datadir,'/plot_imgfiles/',filelab]);
print([datadir,'/plot_imgfiles/',filelab,'/ISRgrid.png'],'-dpng','-r300');

end



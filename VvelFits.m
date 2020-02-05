%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%Vvel FITS%%%%%%%%%%%%%%
%%%%%%TO BE COMPARED TO Butler%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PATH AND FILENAME
clear;
path='./datasets/data_hdf5/';

%DOWNLOADED DATA SET
% filename='20170222.002_lp_1min-fitcal-vvelsLat-60sec.h5';
%filename='20170302.002_lp_2min-fitcal_vvelsLat-60sec_ISINGLASS_B.h5';
filename='20150318.001_lp_1min-vvelsLat-180sec.h5';


%EXTRACT VMAG, EMAG, MAGNETICLATITUDE
vmag=cast(hdf5read([path, filename],'VectorVels/Vmag'),'double');
vdir=cast(hdf5read([path, filename],'VectorVels/Vdir'),'double');
vest=cast(hdf5read([path, filename],'VectorVels/Vest'),'double');
emag=cast(hdf5read([path, filename],'VectorVels/Emag'),'double');
maglat=cast(hdf5read([path, filename],'VectorVels/MagneticLatitude'),'double');
mlat=mean(maglat,1);    %why are we doing this???


%EXTRACT TIME STAMP
dtime=cast(hdf5read([path, filename],'Time/dtime'),'double');
UT=dtime(1,:);   %use beginning of the integration period as time referece.
% month = cast(hdf5read([path, filename],'/Time/Month'),'double');
% day = cast(hdf5read([path, filename],'/Time/Day'),'double');
% year = cast(hdf5read([path, filename],'/Time/Year'),'double');


% %SET UP A DATE STRUCTURE
os=ones(size(UT(:)));
zs=zeros(size(UT(:)));
exp_date=[2017*os, 02*os, 22*os, UT(:), zs, zs];
t=datenum(exp_date);
%exp_date=datevec(datenum(exp_date));


%SIDE BY SIDE VMAG AND EMAG
figure(2);

subplot(2,1,1)
imagesc(t,mlat,squeeze(vest(2,:,:)));
axis xy;
datetick;
%xlim([736748+9/24 736748+16/24]); % base year + decimal hours because of datetick
title('Vvels')
xlabel('UT')
ylabel('mlat.')
clb=colorbar;
caxis ([-1500 1500]);
ylabel(clb,'Mag. east veloc.');

subplot(2,1,2);
t=datenum(exp_date);
imagesc(t,mlat,squeeze(vest(1,:,:)));
axis xy;
datetick;
%xlim([736748+9/24 736748+16/24]); % base year + decimal hours because of datetick
xlabel('UT')
ylabel('mlat.')
clb=colorbar;
caxis ([-1000 1000]);
ylabel(clb,'Mag. north veloc.');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%BUTLER FITS%%%%%%%%%%%%%%
%%%%%%TO BE COMPARED TO VVELS%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

%PROCESSED LONGPULSE DATA
load('./datasets/data_mat/22Feb2017_2.1min_vvelscans.mat');
load('./datasets/data_mat/22Feb2017_2.1min_fieldgrid.mat'); % needed for exp_date information
% load('./datasets/data_mat/02Mar2017_2.1min_vvelscans.mat');
% load(['./datasets/data_mat/02Mar2017_2.1min_fieldgrid.mat']);

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

figure(1); 

subplot(211);
imagesc(t,Yvgcut,Veastcuts);
axis xy;
datetick;
xlim([736748+9/24 736748+16/24]); % base year + decimal hours because of datetick
title('Butler')
caxis([-1.5e3 1.5e3]);
xlabel('UT');
ylabel('dist N (km)');
c=colorbar;
ylabel(c,'eastward drift (m/s)');

subplot(212);
imagesc(t,Yvgcut,Vnorthcuts);
axis xy;
datetick;
xlim([736748+9/24 736748+16/24]); % base year + decimal hours because of datetick
caxis ([-1000 1000]);
xlabel('UT');
ylabel('dist N (km)');
c=colorbar;
ylabel(c,'northward drift (m/s)');
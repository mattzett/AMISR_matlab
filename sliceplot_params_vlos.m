function sliceplot_params_vlos(Rx,Ry,Rz,data1,data2,data3,data4,title1)

lp=50;


xmin=min(Rx(:)); xmax=max(Rx(:));
ymin=min(Ry(:)); ymax=max(Ry(:));
zmin=min(Rz(:)); zmax=max(Rz(:));
% x=linspace(-50,250,lp);
% y=linspace(-50,400,lp);
% z=linspace(100,400,lp);
x=linspace(-50,150,lp);
y=linspace(-50,250,lp);
z=linspace(100,400,lp);
[X,Y,Z]=meshgrid(x,y,z);


datFN=TriScatteredInterp(Rx(:),Ry(:),Rz(:),data1(:));
datacube1=datFN(X(:),Y(:),Z(:));
datacube1=reshape(datacube1,[lp,lp,lp]);

datFN=TriScatteredInterp(Rx(:),Ry(:),Rz(:),data2(:));
datacube2=datFN(X(:),Y(:),Z(:));
datacube2=reshape(datacube2,[lp,lp,lp]);

datFN=TriScatteredInterp(Rx(:),Ry(:),Rz(:),data3(:));
datacube3=datFN(X(:),Y(:),Z(:));
datacube3=reshape(datacube3,[lp,lp,lp]);

datFN=TriScatteredInterp(Rx(:),Ry(:),Rz(:),data4(:));
datacube4=datFN(X(:),Y(:),Z(:));
datacube4=reshape(datacube4,[lp,lp,lp]);

sx=[];
sy=[];
dz=zmax-zmin;
dz=dz/10;
sz=linspace(zmin,zmax,4);
%sz=[120,150,200,250]; %Selected altitudes for plotting, updated 2017/03/31
sz=[120,150,250,325];

set(gcf,'PaperPosition',[0 0 12 4.5]);
az=-62.5;
%% el=24;
%el=18;
el=15;

nlims=[10.6 11.6];
Tlims=[0 3e3];

subplot(6,3,[1:3]);
title(title1,'FontSize',18);
axis off;

here=find(datacube1<=0);
datacube1(here)=1;

% subplot(221);
subplot(6,4,[5:4:17]);
slice(X,Y,Z,log10(datacube1),sx,sy,sz);
view([az,el]);
shading flat;
axis tight;
caxis(nlims);
xlabel('E [km]');
ylabel('N [km]');
zlabel('alt. [km]');
h=colorbar;
ylabel(h,sprintf('log_{10} n_e [m^{-3}]'))

% subplot(222);
subplot(6,4,[6:4:18]);
slice(X,Y,Z,datacube2,sx,sy,sz);
view([az,el]);
shading flat;
axis tight;
caxis(Tlims);
xlabel('E [km]');
ylabel('N [km]');
%zlabel('alt. [km]');
h=colorbar;
ylabel(h,'T_i [K]')

% subplot(223);
subplot(6,4,[7:4:19]);
slice(X,Y,Z,datacube3,sx,sy,sz);
view([az,el]);
shading flat;
axis tight;
caxis(Tlims);
xlabel('E [km]');
ylabel('N [km]');
%zlabel('alt. [km]');
h=colorbar;
ylabel(h,'T_e [K]')

subplot(6,4,[8:4:20]);
slice(X,Y,Z,datacube4,sx,sy,sz);
view([az,el]);
shading flat;
axis tight;
caxis([-500 500]);
xlabel('E [km]');
ylabel('N [km]');
%zlabel('alt. [km]');
h=colorbar;
ylabel(h,'v_{los} [m/s]')

end

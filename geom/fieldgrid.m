function [Bx,By,Bz]=fieldgrid(Rx,Ry,Rz,zentheta,zenphi,zline)

%DEFINE A 3D FIELD LINE GRID CORRESPONDING TO SET OF 2D OBSERVATION POINTS

%BEAM E-REGION INTERSECTION IN O
xB=Rx(:); yB=Ry(:); zB=Rz(:);


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


%FIELD LINE GRID INDEXED BY ALTITUDE AND POSITION NUMBER
lz=length(zline);
lp=length(Rx);
Bx=zeros(lz,lp); By=Bx; Bz=Bx;
for ip=1:lp
    t=(zline-zC(ip))/mz(ip);
    
    Bx(:,ip)=mx(ip)*t+xC(ip);
    By(:,ip)=my(ip)*t+yC(ip);
    Bz(:,ip)=mz(ip)*t+zC(ip);
end

end
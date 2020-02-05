function grid_ISR(datadir,filelab)

%This script computes positions of each beam in a local ENU coordinate
%system.  This system is used by the velocity fitting routine.  
% 
% Development notes:
%    1) This should also be redone as a function that accepts the isr beam 
%    data and returns the radar
%    coordinates (R*) and also the field line grid (B* - not currently used
%    but could be helpful in the future).  The end part of the program
%    interpolates onto a field line grid, which we don't need for now, but
%    keep the commented code in case it is useful later on.

load([datadir,'/data_mat/',filelab,'_rawdata.mat']); % This is redundant but useful when working the code in segments

% tic
% fprintf('ISR2FIELDGRID.M --> Calculating field lines\n');
lp=length(az);

%MAGNETIC ZENITH
zenaz=205.7;
zenel=77.5;
[zentheta,zenphi]=azel2thetaphi(zenaz,zenel);


%CONVERT ANGLES
[radtheta,radphi]=azel2thetaphi(az,el);


%WHICH BEAM INDEX CORRESPONDS TO MAGNETIC ZENITH
[mintheta,itheta]=min(abs(radtheta-zentheta));
[minphi,iphi]=min(abs(radphi-zenphi));
if(mintheta<0.01 & minphi<0.01 & itheta==iphi)
  ibupB=itheta;
else
  ibupB=0;
end


%RADAR GRID
lr=size(isz,1);
Rx=zeros(lr,lp); Ry=Rx; Rz=Rx;
for ip=1:lp
    range=isz(:,ip)*sec(radtheta(ip));
    
    Rx(:,ip)=range.*sin(radtheta(ip)).*cos(radphi(ip));
    Ry(:,ip)=range.*sin(radtheta(ip)).*sin(radphi(ip));
    Rz(:,ip)=range.*cos(radtheta(ip));
end


%ESTABLISH A REASONABLE FIELD LINE GRID
zline=100:10:400;
[~,iz]=min(abs(Rz(:,1)-150));
[Bx,By,Bz]=fieldgrid(Rx(iz,:)',Ry(iz,:)',Rz(iz,:)',zentheta,zenphi,zline);


%ASPECT ANGLE BETWEEN BEAM AND FIELD
aspect_ang=zeros(lp,1);
for ip=1:lp
    Rpos=[sin(radtheta(ip))*cos(radphi(ip)); sin(radtheta(ip))*sin(radphi(ip)); cos(radtheta(ip))];
    Bpos=[sin(zentheta)*cos(zenphi); sin(zentheta)*sin(zenphi); cos(zentheta)];
    
    aspect_ang(ip)=acos(Rpos'*Bpos);
end
ang=repmat(aspect_ang',lr,1);


%SAVE THE GRID INFORMATION
filelab=datestr(exp_date(1,:),'ddmmmyyyy');
save([datadir,'/data_mat/',filelab,'_',num2str(tint),'min','_fieldgrid.mat'],'exp_date','az','el','R*','B*','tint','zen*');


% rmpath ./geom

%{
%INTERPOLATE ASPECT ANGLE ONTO FIELD LINE GRID (ONLY AN APPROXIMATE)
lfl=size(Bx,2);
lfr=size(Bx,1);
angFN=TriScatteredInterp(Rx(:),Ry(:),Rz(:),ang(:));
angp=angFN(Bx(:),By(:),Bz(:));
angp=reshape(angp,lfr,lfl);
inan=find(isnan(angp));
angp(inan)=0;            %kludge
isangB=angp;


%LOOP OVER DATASET, APPLYING LINEAR INTERPOLATION OVER IRREGULAR GRID
lt=length(exp_date);
istiB=zeros(lfr,lfl,lt); isteB=istiB; isviB=istiB;
isneB=istiB; isdtiB=istiB; isdteB=istiB; isdviB=istiB;
isdneB=istiB; ispB=istiB;
for it=1:lt
    if(mod(it,10) == 0)
        fprintf('ISR2FIELDGRID.M --> Interpolating scan %d \n',it);
    end
    
    
    %GRAB SCAN FROM THIS PARTICULAR TIME INTERVAL
    Ti=squeeze(isti(:,:,it));
    Te=squeeze(iste(:,:,it));
    vi=squeeze(isvi(:,:,it));
    ne=squeeze(isne(:,:,it));
    p=squeeze(isp(:,:,it));
    dTi=squeeze(isdti(:,:,it));
    dTe=squeeze(isdte(:,:,it));
    dvi=squeeze(isdvi(:,:,it));
    dne=squeeze(isdne(:,:,it));
    
    
    %CALCULATE INTERPOLANT FUNCTIONS FOR EACH PARAMETER
    TiFN=TriScatteredInterp(Rx(:),Ry(:),Rz(:),Ti(:));
    TeFN=TriScatteredInterp(Rx(:),Ry(:),Rz(:),Te(:));
    viFN=TriScatteredInterp(Rx(:),Ry(:),Rz(:),vi(:));
    neFN=TriScatteredInterp(Rx(:),Ry(:),Rz(:),ne(:));
    pFN=TriScatteredInterp(Rx(:),Ry(:),Rz(:),p(:));
    dTiFN=TriScatteredInterp(Rx(:),Ry(:),Rz(:),dTi(:));
    dTeFN=TriScatteredInterp(Rx(:),Ry(:),Rz(:),dTe(:));
    dviFN=TriScatteredInterp(Rx(:),Ry(:),Rz(:),dvi(:));
    dneFN=TriScatteredInterp(Rx(:),Ry(:),Rz(:),dne(:));
    
    
    %APPLY INTERPOLANT FOR FIELD LINE GRID
    Tip=TiFN(Bx(:),By(:),Bz(:));
    Tep=TeFN(Bx(:),By(:),Bz(:));
    vip=viFN(Bx(:),By(:),Bz(:));
    nep=neFN(Bx(:),By(:),Bz(:));
    pp=pFN(Bx(:),By(:),Bz(:));
    dTip=dTiFN(Bx(:),By(:),Bz(:));
    dTep=dTeFN(Bx(:),By(:),Bz(:));
    dvip=dviFN(Bx(:),By(:),Bz(:));
    dnep=dneFN(Bx(:),By(:),Bz(:));
    
    
    %PROPERLY SHAPE ARRAYS
    Tip=reshape(Tip,lfr,lfl);
    Tep=reshape(Tep,lfr,lfl);
    vip=reshape(vip,lfr,lfl);
    nep=reshape(nep,lfr,lfl);
    pp=reshape(pp,lfr,lfl);
    dTip=reshape(dTip,lfr,lfl);
    dTep=reshape(dTep,lfr,lfl);
    dvip=reshape(dvip,lfr,lfl);
    dnep=reshape(dnep,lfr,lfl);    
            
    
    %TREAT UP B BEAM (IF IT EXISTS) AS A SPECIAL CASE B/C TRIANGULATION
    %TENDS TO SCREW THIS UP OTHERWISE.
    if(ibupB ~= 0)
        Tip(:,ibupB)=max(interpolate(Ti(:,ibupB),Rz(:,ibupB),Bz(:,ibupB),'lin','lin'),0);
        Tep(:,ibupB)=max(interpolate(Te(:,ibupB),Rz(:,ibupB),Bz(:,ibupB),'lin','lin'),0);
        vip(:,ibupB)=interpolate(vi(:,ibupB),Rz(:,ibupB),Bz(:,ibupB),'lin','lin');
        nep(:,ibupB)=max(interpolate(ne(:,ibupB),Rz(:,ibupB),Bz(:,ibupB),'lin','lin'),0);
        pp(:,ibupB)=min(max(interpolate(p(:,ibupB),Rz(:,ibupB),Bz(:,ibupB),'lin','lin'),0),1);
        dTip(:,ibupB)=max(interpolate(dTi(:,ibupB),Rz(:,ibupB),Bz(:,ibupB),'lin','lin'),0);
        dTep(:,ibupB)=max(interpolate(dTe(:,ibupB),Rz(:,ibupB),Bz(:,ibupB),'lin','lin'),0);
        dvip(:,ibupB)=max(interpolate(dvi(:,ibupB),Rz(:,ibupB),Bz(:,ibupB),'lin','lin'),0);
        dnep(:,ibupB)=max(interpolate(dne(:,ibupB),Rz(:,ibupB),Bz(:,ibupB),'lin','lin'),0);
    end
    
    
    %CORRECT EXTRAPOLATED PARTS OF EACH FIELD LINE PROFILE.  OR PARTS WHOSE
    %INTERPOLATION PRODUCED A NEGATIVE RESULT FOR TI,TE,NE.
    for ifl=1:lfl
        inan=find(isnan(Tip(:,ifl)) | Tip(:,ifl)<0 | Tep(:,ifl)<0 | nep(:,ifl)<0);
        
        if length(inan) ~= 0
            if length(inan)>=lfr-1    %data is total crap
                Tip(:,ifl)=zeros(size(Tip(:,ifl)));
                Tep(:,ifl)=zeros(size(Tip(:,ifl)));
                vip(:,ifl)=zeros(size(Tip(:,ifl)));
                nep(:,ifl)=zeros(size(Tip(:,ifl)));
                pp(:,ifl)=zeros(size(Tip(:,ifl)));
                dTip(:,ifl)=ones(size(Tip(:,ifl)));
                dTep(:,ifl)=ones(size(Tip(:,ifl)));
                dvip(:,ifl)=ones(size(Tip(:,ifl)));
                dnep(:,ifl)=ones(size(Tip(:,ifl)));
            else
                idef=setdiff(1:lfr,inan);
                zfl=Bz(inan,ifl);
                ztmp=Bz(idef,ifl);
                
                Tiptmp=Tip(idef,ifl);
                Teptmp=Tep(idef,ifl);
                viptmp=vip(idef,ifl);
                neptmp=nep(idef,ifl);
                pptmp=pp(idef,ifl);
                dTiptmp=dTip(idef,ifl);
                dTeptmp=dTep(idef,ifl);
                dviptmp=dvip(idef,ifl);
                dneptmp=dnep(idef,ifl);
                
                Tipcorr=max(interpolate(Tiptmp,ztmp,zfl,'lin','lin'),0);
                Tepcorr=max(interpolate(Teptmp,ztmp,zfl,'lin','lin'),0);
                vipcorr=interpolate(viptmp,ztmp,zfl,'lin','lin');
                nepcorr=max(interpolate(neptmp,ztmp,zfl,'lin','lin'),0);
                ppcorr=min(max(interpolate(pptmp,ztmp,zfl,'lin','lin'),0),1);
                dTipcorr=max(interpolate(dTiptmp,ztmp,zfl,'lin','lin'),0);
                dTepcorr=max(interpolate(dTeptmp,ztmp,zfl,'lin','lin'),0);
                dvipcorr=max(interpolate(dviptmp,ztmp,zfl,'lin','lin'),0);
                dnepcorr=max(interpolate(dneptmp,ztmp,zfl,'lin','lin'),0);
                
                Tip(inan,ifl)=Tipcorr;
                Tep(inan,ifl)=Tepcorr;
                vip(inan,ifl)=vipcorr;
                nep(inan,ifl)=nepcorr;
                pp(inan,ifl)=ppcorr;
                dTip(inan,ifl)=dTipcorr;
                dTep(inan,ifl)=dTepcorr;
                dvip(inan,ifl)=dvipcorr;
                dnep(inan,ifl)=dnepcorr;
            end
        end
    end
    
    
    %STORE COMPLETED FIELD LINE PROFILES
    istiB(:,:,it)=Tip;
    isteB(:,:,it)=Tep;
    isviB(:,:,it)=vip;
    isneB(:,:,it)=nep;
    ispB(:,:,it)=pp;
    isdtiB(:,:,it)=dTip;
    isdteB(:,:,it)=dTep;
    isdviB(:,:,it)=dvip;
    isdneB(:,:,it)=dnep;
end


%CREATE A NEW DATA SET USING FIELD LINE GRID
filelab=datestr(exp_date(1,:),'ddmmmyyyy');
save(['./datasets/data_mat/',filelab,'_',num2str(tint),'min','_fieldgrid.mat'], ...
    'is*B','exp_date','az','el','R*','B*','tint','zen*');


rmpath ./geom ./optim_fns;
toc
%}

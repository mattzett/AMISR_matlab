function fileout=hdf5_extract(filename,saveflag)
%This script reads in the amisr hdf5 files that are posted on SRIs website:
% http://amisr.com/amisr/links/data-access/  .  Note that, depending on when
% the data were fitted, the exact parameters in the "Fits" array may
% differ and as far as I can tell the only real way to know what is there
% is to look at the metadata in the hdf5 file.

%PATH AND FILENAME
path='./datasets/data_hdf5/';


%LOOK INSIDE THE h5 FILE
% h5disp([path,filename])


%AZ,EL,ALT FOR EXPERIMENT
BeamCodes=cast(hdf5read([path, filename],'BeamCodes'),'double');
beamcodes=BeamCodes(1,:);
az=BeamCodes(2,:);
el=BeamCodes(3,:);
z=cast(hdf5read([path, filename],'/Geomag/Altitude'),'double')./1e3;
beamalt=cast(hdf5read([path, filename],'/Geomag/Altitude'),'double');
beamlat=cast(hdf5read([path, filename],'/Geomag/Latitude'),'double');
beamlon=cast(hdf5read([path, filename],'/Geomag/Longitude'),'double');


%READ FITTED PARAMETERS AND TIME.  DIMENSIONS:
% param  ion  range  beam  time
%   4     3    17    26   698,   these are not likely to be the same for every experiment
Ne_in = cast(hdf5read([path, filename],'/FittedParams/Ne'),'double');
dtime = cast(hdf5read([path, filename],'/Time/dtime'),'double');
month = cast(hdf5read([path, filename],'/Time/Month'),'double');
day = cast(hdf5read([path, filename],'/Time/Day'),'double');
year = cast(hdf5read([path, filename],'/Time/Year'),'double');
Fits = cast(hdf5read([path, filename],'/FittedParams/Fits'),'double');
Errors = cast(hdf5read([path, filename],'/FittedParams/Errors'),'double');
Range = cast(hdf5read([path, filename],'/FittedParams/Range'),'double');


%TRUNCATE RECORDS TO SAME SIZE
inan=[];
lb=size(BeamCodes,2);
for ib=1:lb
    inannew=find(isnan(z(:,ib)));
    if length(inannew) > length(inan)
        inan=inannew;
    end
end
izs=1:min(inan)-1;


%SORT INTO MATRICES
isvi=squeeze(Fits(4,1,izs,:,:));
isti=squeeze(Fits(2,1,izs,:,:));
iste=squeeze(Fits(2,6,izs,:,:));
isne=squeeze(Ne_in(izs,:,:));
isp=squeeze(Fits(1,1,izs,:,:));

isdvi=squeeze(Errors(4,1,izs,:,:));
isdti=squeeze(Errors(2,1,izs,:,:));
isdte=squeeze(Errors(2,6,izs,:,:));
isdne=squeeze(Errors(4,1,izs,:,:));

isz=z(izs,:);
UT=dtime(1,:);


%ALSO GRAB THE SNR FROM THE NEPOWER HDF5GROUP (IT HAS ITS OWN RANGE GRID)
SNRfull=cast(hdf5read([path,filename],'/NeFromPower/SNR'),'double');
RangeSNR=cast(hdf5read([path,filename],'/NeFromPower/Range'),'double');


%INTERPOLATE THE SNR TO EACH BEAM RANGE
[lr,lb,lt]=size(isne);
for it=1:lt
   for ib=1:lb
      SNR(:,ib,it)=interp1(RangeSNR(:),SNRfull(:,ib,it),Range(:,ib));
   end
end
SNR=SNR(izs,:,:);


%SET UP A DATE STRUCTURE
os=ones(size(UT(:)));
zs=zeros(size(UT(:)));
exp_date=[year(1,:)', month(1,:)', day(1,:)', UT(:), zs, zs];
exp_date=datevec(datenum(exp_date));


%INTEGRATION TIME FOR THESE DATA
tint=60*(UT(2)-UT(1));
tint=round(10*tint)/10;


%SAVE TO A .MAT FILE FOR FURTHER PROCESSING
filelab=datestr(exp_date(1,:),'ddmmmyyyy');
if saveflag==1
    save(['./datasets/data_mat/',filelab,'_',num2str(tint),'min','_rawdata.mat'],'is*','exp_date','az','el','tint','beam*','SNR');
end
fileout=[filelab,'_',num2str(tint),'min'];
end





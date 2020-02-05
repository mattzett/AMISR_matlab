function strlab=datelab(UT,dmy)


UTsec=round(UT*3600);
if UTsec==86400
    UTstr='00000';
    [tmp,dmynew]=dateinc(UT,dmy,1);
else
    UTstr=num2str(round(UT*3600));
    lstr=numel(UTstr);
    while lstr<5
        UTstr=['0',UTstr];
        lstr=lstr+1;
    end
    dmynew=dmy;
end

yearstr=num2str(dmynew(3));
lstr=numel(yearstr);
while lstr<4
    yearstr=['0',yearstr];
    lstr=lstr+1;
end

monthstr=num2str(dmynew(2));
lstr=numel(monthstr);
while lstr<2
    monthstr=['0',monthstr];
    lstr=lstr+1;
end

daystr=num2str(dmynew(1));
lstr=numel(daystr);
while lstr<2
    daystr=['0',daystr];
    lstr=lstr+1;
end

strlab=[yearstr,monthstr,daystr,'_',UTstr];

end
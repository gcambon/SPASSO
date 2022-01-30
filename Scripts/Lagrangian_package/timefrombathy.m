function [touched,lons,lats,lonf,latf,touchedlon,touchedlat]=timefrombathy(day0,lons0,lats0,bathylvl,numdays,product,bathy_product)
%timeplateau=fromKplateau(day0,lons0,lats0,bathylvl,numdays,product)
%Gives the first day at which each particle has touched for the first time the bathymentric level bathylvl
%-1 if never touched.

[lonf,latf,lons,lats]=aviso_lonlatadvect(day0,lons0,lats0,numdays,product);

if(nargin==6)
	bathy_product='hr';
	end


Nstep=round(abs(numdays)/(60*60*24)*8)+2;

lons=lons(1:8:end,:);%precision of 1 day
lats=lats(1:8:end,:);
sz=size(lons);
touched=zeros(size(lons0));

trajdepth=oceandepth(lons,lats,bathy_product);
trajdepthM=reshape(trajdepth,sz);

for ct=1:length(lons0(:))
        
        [tmp,touch]=max(trajdepthM(:,ct)>bathylvl);
        if(tmp==0),touch=-1;,end
        touched(ct)=touch;
        end


touchedpos=touched;
nottouched=find(touchedpos(:)<0); %before was <2
touchedpos(nottouched)=1;
touchedlat=zeros(size(touched));
touchedlon=zeros(size(touched));

for ct=1:length(touchedpos(:))
touchedlat(ct)=lats(touchedpos(ct),ct);
touchedlon(ct)=lons(touchedpos(ct),ct);
end

touchedlat(nottouched)=NaN;
touchedlon(nottouched)=NaN;


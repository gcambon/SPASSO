function [Udegs,Vdegs,lonw,latw,Ucms,Vcms,u,v,lon,lat]=getncUV_aviso_cmems_product_med(path_to_product,nday,wrk_dir)
%[Udegs,Vdegs,lonw,latw,Ucms,Vcms,u,v,lon,lat]=getncUV_aviso_product(path_to_product,nday) reads any Aviso u,v file.
%nday is in the format: [2010 12 27]. u,v,lo,lat are the variables as in the netcdf file. Udegs, Vdegs are the currents in deg/s.
%If the product covers the longitudinally the entire globe, two latitudinal bands at 360 are repeated at O deg, in order to have trajectories over a sphere.
%Should be campatible with Matlab 

lvars=localvars(wrk_dir);

RT=6371e5;

%From e.g. [2010 12 27] to 20101227
ndaystr=sprintf("%d%02d%02d",nday(1),nday(2),nday(3));

%fnames=glob(sprintf("%s/*_allsat_phy_l4_%s_*.nc.gz",lvars.datadir,ndaystr));%STEPH
fnames=glob(sprintf("%s/*_allsat_phy_l4_%s_*.nc",lvars.datadir,ndaystr));

%If fnames is empty, try to find the files in a directory tree without the year in the path (this is the convenction of NRT files, for instance)
%if(isempty(fnames))
%fnames=glob(sprintf("%s/products/ftp.aviso.oceanobs.com/%s/*_%s_*.nc.gz",lvars.datadir,path_to_product,ndaystr));
%endif

 if(isempty(fnames))
                disp("No Aviso file found.");
                Udegs=[];,Vdegs=[];,lonw=[];,latw=[];,Ucms=[];,Vcms=[];,u=[];,v=[];,lon=[];,lat=[];
		return; 
                else
                fname=fnames{end};     
                endif
fname
ncfiletmp=sprintf("%s/%s_%08d_unzipped",lvars.datadir,ndaystr,rand(1)*10000000);

%cmd=sprintf("gunzip -c %s >%s",fname,ncfiletmp);
cmd=sprintf("cp %s %s",fname,ncfiletmp);%STEPH
system(cmd); %STEPH

lon=ncread(ncfiletmp,"longitude");
lat=ncread(ncfiletmp,"latitude");
u=ncread(ncfiletmp,"ugos");
v=ncread(ncfiletmp,"vgos");

% Convert 2014 aviso variables format
u=squeeze(u);
v=squeeze(v);
Ucms=u*100;
Vcms=v*100; %conversion m/s to cm/s

% Create lon lat matrices
nlon=length(lon);
nlat=length(lat);

lon=reshape(lon,nlon,1);
lat=reshape(lat,nlat,1);


latw=lat;
%FInding whether lon goes all around the globe
lonw=lon;

% if dataset is for longitude [-180;180], use this script like that
%% if dataset is for longitude [0;360], umcomment all lines wich are commented with the mention '% steph uncommented'


%fw=find(lonw>180); %AR 03/04/2019 commented % steph uncommented
%lonw(fw)=lonw(fw)-360; %AR 03/04/2019 commented % steph uncommented
%[lonw,lonsrti]=sort(lonw); %AR 03/04/2019 commented % steph uncommented
%lonw=lonw([(nlon-1) nlon 1:nlon 1]); %AR 03/04/2019 commented % steph uncommented
%lonw(1:2)=lonw(1:2)-360; %AR 03/04/2019 commented % steph uncommented
%lonw(end)=lonw(end)+360; %AR 03/04/2019 commented % steph uncommented



%if(sum(diff(diff(lonw)))<1e-9) %If the raccorded lonw is growing with constant step, lon goes around all the globe %AR 03/04/2019 commented % steph uncommented

%Ucms=Ucms(lonsrti,:);%AR 03/04/2019 commented % steph uncommented
%Vcms=Vcms(lonsrti,:);%AR 03/04/2019 commented % steph uncommented
%Ucms=Ucms([(nlon-1) nlon 1:nlon 1],:);%AR 03/04/2019 commented % steph uncommented
%Vcms=Vcms([(nlon-1) nlon 1:nlon 1],:);%AR 03/04/2019 commented % steph uncommented

%else %If the raccorded lonw has a break, then lon was regional: hence resetting lonw to lon %AR 03/04/2019 commented % steph uncommented
%lonw=lon;%AR 03/04/2019 commented % steph uncommented
%end %AR 03/04/2019 commented % steph uncommented

if(max(lonw(:))>360)
	lonw=lonw-360;
	end

% Velocities in deg/s
[glon,glat]=meshgrid(lonw,latw);
Udegs=Ucms./(RT.*cos(glat'./180*pi)).*180/pi;
Vdegs=Vcms*180/pi/RT;

system(sprintf("rm %s",ncfiletmp));




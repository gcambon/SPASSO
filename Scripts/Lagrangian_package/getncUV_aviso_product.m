function [Udegs,Vdegs,lonw,latw,Ucms,Vcms,u,v,lon,lat]=getncUV_aviso_product(path_to_product,nday,wrk_dir)
%[Udegs,Vdegs,lonw,latw,Ucms,Vcms,u,v,lon,lat]=getncUV_aviso_product(path_to_product,nday) reads any Aviso u,v file.
%nday is in the format: [2010 12 27]. u,v,lo,lat are the variables as in the netcdf file. Udegs, Vdegs are the currents in deg/s.
%If the product covers the longitudinally the entire globe, two latitudinal bands at 360 are repeated at O deg, in order to have trajectories over a sphere.
%Should be campatible with Matlab 

lvars=localvars(wrk_dir);



RT=6371e5;

%From e.g. [2010 12 27] to 20101227
disp(path_to_product);

ndaystr=sprintf("%d%02d%02d",nday(1),nday(2),nday(3));
fnames=glob(sprintf("%s/%s/*_%s_*.nc.gz",lvars.datadir,path_to_product,ndaystr));


%disp(sprintf("%s/%s/*_%s_*.nc.gz",lvars.datadir,path_to_product,ndaystr));
%disp(fnames)

 if(isempty(fnames))
                disp("No Aviso file found.");
                Udegs=[];,Vdegs=[];,lonw=[];,latw=[];,Ucms=[];,Vcms=[];,u=[];,v=[];,lon=[];,lat=[];
		return; 
                else
                fname=fnames{end};     
                endif
ncfiletmp=sprintf("%s/%s_%08d_unzipped",lvars.datadir,ndaystr,rand(1)*10000000);

cmd=sprintf("gunzip -c %s >%s",fname,ncfiletmp);
system(cmd);

lon=ncread(ncfiletmp,"lon");
lat=ncread(ncfiletmp,"lat");
u=ncread(ncfiletmp,"u");
v=ncread(ncfiletmp,"v");
% Convert 2014 aviso variables format
u=squeeze(u);
v=squeeze(v);
Ucms=u*100;
Vcms=v*100; %conversion m/s to cm/s

Ucms(isnan(Ucms))=0;
Vcms(isnan(Vcms))=0;

% Create lon lat matrices
nlon=length(lon);
nlat=length(lat);

lon=reshape(lon,nlon,1);
lat=reshape(lat,nlat,1);


latw=lat;

%FInding whether lon goes all around the globe

lonw=lon;

fw=find(lonw>180);
lonw(fw)=lonw(fw)-360;
[lonw,lonsrti]=sort(lonw);
lonw=lonw([(nlon-1) nlon 1:nlon 1]);
lonw(1:2)=lonw(1:2)-360;
lonw(end)=lonw(end)+360;



if(sum(diff(diff(lonw)))<1e-9) %If the raccorded lonw is growing with constant step, lon goes around all the globe 

Ucms=Ucms(lonsrti,:);
Vcms=Vcms(lonsrti,:);
Ucms=Ucms([(nlon-1) nlon 1:nlon 1],:);
Vcms=Vcms([(nlon-1) nlon 1:nlon 1],:);

else %If the raccorded lonw has a break, then lon was regional: hence resetting lonw to lon
lonw=lon;
end

if(max(lonw(:))>360)
	lonw=lonw-360;
	end


% Velocities in deg/s
[glon,glat]=meshgrid(lonw,latw);
Udegs=Ucms./(RT.*cos(glat'./180*pi)).*180/pi;
Vdegs=Vcms*180/pi/RT;

system(sprintf("rm %s",ncfiletmp));




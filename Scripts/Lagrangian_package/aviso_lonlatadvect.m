function [lonf,latf,lons,lats]=aviso_lonlatadvect(dayv,lon,lat,numdays,product)
% Compute longitude and latitude advection
day0=datenum(dayv)-datenum([1950 1 1]);


dayf=day0+numdays;
aviso_load(min([day0,dayf]),max([day0,dayf]),product);

pts=zeros(length(lon(:))*2,1);
numpts=length(lon(:));
pts(1:2:end)=lon(:);
pts(2:2:end)=lat(:);

tspan=([day0 dayf]'*60*60*24);
Nstep=round(abs(diff(tspan))/(60*60*24)*8)+2;
trj=RK4(tspan,pts,Nstep);
[lons,lats]=trj2pos(trj,Nstep,numpts);
latf=trj(2*Nstep:Nstep*2:end);
lonf=trj(Nstep:Nstep*2:end);

%out=FSLEn_v(tspan,pts,3,delta0,Nstep,delta);


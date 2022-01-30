function [lon0,lat0,lonf,latf,sz]=aviso_bigsqadvect(dayv,lonv,latv,delta0,numdays,product)

day0=datenum(dayv)-datenum([1950 1 1]);


dayf=day0+numdays;
aviso_load(min([day0,dayf]),max([day0,dayf]),product);

xvc=[lonv(1):delta0:lonv(end)];
yvc=[latv(1):delta0:latv(end)];

[xg,yg]=meshgrid(xvc,yvc);

sz=size(xg);
lonf=zeros(size(xg));
latf=zeros(size(xg));
pts=zeros(sz(2)*2,1);

for ctsz=1:sz(1);


	pts(1:2:end)=xg(ctsz,:);
	pts(2:2:end)=yg(ctsz,:);

	tspan=([day0 dayf]'*60*60*24);
	Nstep=round(abs(diff(tspan))/(60*60*24)*8)+2;
	trj=RK4(tspan,pts,Nstep);

	latf(ctsz,:)=trj(2*Nstep:Nstep*2:end);
	lonf(ctsz,:)=trj(Nstep:Nstep*2:end);
end

lon0=xg(:);
lat0=yg(:);
%out=FSLEn_v(tspan,pts,3,delta0,Nstep,delta);


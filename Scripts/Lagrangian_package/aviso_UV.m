function [U,V,xvc,yvc,Ucms,Vcms]=aviso_UV(dayc,lonv,latv,delta0,product)
% Get AVISO uv, reshape matrix and convert u and v into cm/s
RT=6371e5;

deltad0=8;
deltadf=8;

daycc=datenum(dayc)-datenum([1950 1 1]);

if(~strcmp(product,'none'))

	switch product
		case {'nrt_global' 'rt_global' 'rt_med' 'nrt_med','nrt_global_2010','nrt_kerguelen','nrt_global_2014'}
			deltad0=1;
			deltadf=0;
	end

	day0=datenum(dayc)-datenum([1950 1 1])-deltad0;
	dayf=datenum(dayc)-datenum([1950 1 1])+deltadf;

	aviso_load(floor(day0),ceil(dayf),product);	
end

% lonv=([0 360]);
% latv=([-90 90]);

xvc=[lonv(1):delta0:lonv(end)];
yvc=[latv(1):delta0:latv(end)];

[xg,yg]=meshgrid(xvc,yvc);
pts=zeros(length(xg(:))*2,1);

pts(1:2:end)=xg(:);
pts(2:2:end)=yg(:);

t=daycc*24*60*60;

out=UVext_v(t,pts);
U=reshape(out(1,:),length(yvc),length(xvc));
V=reshape(out(2,:),length(yvc),length(xvc));

Ucms=U*pi/180.*(RT.*cos(yg./180*pi));
Vcms=V/180*pi*RT;

function [ow,dUdx,dUdy,dVdx,dVdy,U,V]=owtraj(tv,lonv,latv,ds)

sz=size(latv);

ow=zeros(size(latv));
dUdx=zeros(size(latv));
dUdy=zeros(size(latv));
dVdx=zeros(size(latv));
dVdy=zeros(size(latv));
U=zeros(size(latv));
V=zeros(size(latv));


RT=6371e5;
for ct=1:length(latv(:))

	conv=[pi/180.*(RT.*cos(latv(ct)./180*pi)); 1/180*pi*RT];
	out=xyp(tv(ct),lonv(ct),latv(ct));
	outx=xyp(tv(ct),lonv(ct)+ds,latv(ct));
	outy=xyp(tv(ct),lonv(ct),latv(ct)+ds);
	outmx=xyp(tv(ct),lonv(ct)-ds,latv(ct));
	outmy=xyp(tv(ct),lonv(ct),latv(ct)-ds);

	vel=out.*conv;
	velx=outx.*conv;
	vely=outy.*conv;
	velmx=outmx.*conv;
	velmy=outmy.*conv;

	dsx=2*ds*conv(1);
	dsy=2*ds*conv(2);

	dUdx(ct)=(velx(1)-velmx(1))/dsx;
	dUdy(ct)=(vely(1)-velmy(1))/dsy;
	dVdx(ct)=(velx(2)-velmx(2))/dsx;
	dVdy(ct)=(vely(2)-velmy(2))/dsy;
	U(ct)=vel(1);
	V(ct)=vel(2);

end


sn=dUdx-dVdy;
ss=dUdy+dVdx;
vor=-dUdy+dVdx;
ow=sn.^2+ss.^2-vor.^2;


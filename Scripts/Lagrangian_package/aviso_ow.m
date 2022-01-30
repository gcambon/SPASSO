function [ow,dUdx,dUdy,dVdx,dVdy,Ucms,Vcms]=aviso_ow(dayc,lonv,latv,ds,product)
% Compute Okubo-Weiss parameter

R=6371e5;
% Convert velocities to cm/sec
% (for Lagrangian analysis velocities are in deg/sec)
[U,V,lonw,latw,Ucms,Vcms]=aviso_UV(dayc,lonv,latv,ds,product);

[glon,glat]=meshgrid(lonw,latw);

%dx=(glon(2:end,:)-glon(1:(end-1),:))/180*pi*R.*cos(glat(2:end,:)*pi/180);
%dy=(glat(:,2:end)-glat(:,1:(end-1)))/180*pi*R;


[dUdx,dUdy]=gradient(Ucms);
[dVdx,dVdy]=gradient(Vcms);

[Dx,tmp]=gradient(glon/180*pi*R.*cos(glat*pi/180));
[tmp,Dy]=gradient(glat/180*pi*R);

dUdx=dUdx./Dx;
dVdx=dVdx./Dx;
dUdy=dUdy./Dy;
dVdy=dVdy./Dy;



%[dudx,dudy]=gradient(U',lonw,latw);
%[dvdx,dvdy]=gradient(V',lonw,latw);

sn=dUdx-dVdy;
ss=dUdy+dVdx;
vor=-dUdy+dVdx;

ow=sn.^2+ss.^2-vor.^2;

%figure,pcolor(glon,glat,min(max(ow,-1e-10),1e-10));,colormap('jet');

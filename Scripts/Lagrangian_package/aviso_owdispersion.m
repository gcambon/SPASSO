function posexit=aviso_owdispersion(day0,lonv,latv,delta0,numdays,tollerance,product)


dayf=datevec(datenum(day0)+numdays);
daycnes=datenum([1950 1 1]);
ds=1/6;

xvc=[lonv(1):delta0:lonv(end)];
yvc=[latv(1):delta0:latv(end)];
[xg,yg]=meshgrid(xvc,yvc);
pts=zeros(length(xg(:))*2,1);
pts(1:2:end)=xg(:);
pts(2:2:end)=yg(:);

aviso_load(min([datenum(dayf) datenum(day0)])-daycnes,max([datenum(dayf) datenum(day0)])-daycnes,product);

day0n=datenum(day0)-daycnes;
dayfn=datenum(day0)-daycnes+numdays;

tspan=([day0n dayfn]'*60*60*24);
Nstep=round(abs(diff(tspan))/(60*60*24)*8)+2;



trj=RK4(tspan,pts,Nstep);

%[ssth,sstf,lonsstf,latsstf,sst]=amsre_adv(day0,lonv,latv,delta0,-3,0,0,'none');

[lons,lats,M]=trj2pos(trj,Nstep,length(xg(:)));
lons=lons(1:8:end,:);%precision of 1 day
lats=lats(1:8:end,:);

sz=size(lons);

Nstepe=ceil(Nstep/8); %precision of 1 day
tall=(((1:Nstepe)-1)'./(Nstepe-1)*(dayfn-day0n)+day0n)*ones(1,sz(2));

[ow,dUdx,dUdy,dVdx,dVdy,U,V]=owtraj(tall(:)*60*60*24,lons(:),lats(:),ds);

owM=reshape(ow,sz);

posexit=zeros(size(xg));

for ct=1:length(xg(:))

	[tmp,exit1]=max(sign(owM(tollerance:end,ct)));
	if(tmp<0)
		exit1=sz(1);
	end
	posexit(ct)=exit1;
	
end


%figure,imagesc(lonv,latv, posexit),axis xy

%owm=aviso_ow(day0,lonv,latv,ds,'none');
%figure,imagesc(lonv,latv, owm),axis xy
%figure,imagesc(lonv,latv,reshape(owM(1,:),size(xg))),axis xy



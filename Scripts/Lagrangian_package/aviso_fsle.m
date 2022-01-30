function [lambda,xvc,yvc,out]=aviso_fsle(dayv,lonv,latv,delta0,delta,product,numdays)
% Compute FSLEs using AVISO uv product
% Compute FSLEs using:
%
% - dayv is the day when the analysis is started
% - lonv and latv are the longitude latitude limits of the domain
% - delta0 is the resolution of the FSLE field & 
%                                 the initial distance of separation
% - delta is the final distance of separation
% - product is the velocity product to be used in the analysis:
%                       if product='none' then no new velocity field is loaded
%                       and the one already in memory is used
% - days is the total number of days when to compute trajectories

% Set total days if not defined
if(nargin<7)
	numdays=200;
end
% Compute initial and final day
% (day0>dayf => backward FSLE !!! )
dayf=datenum(dayv)-datenum([1950 1 1]);
day0=dayf-numdays;

% load velocity field
% (nothing loaded if product='none')
aviso_load(floor(day0),ceil(dayf),product);	
% lon lat matrices to deploy particles
xvc=[lonv(1):delta0:lonv(end)];
yvc=[latv(1):delta0:latv(end)];
[xg,yg]=meshgrid(xvc,yvc);
% Create vector of particle position to compute FSLE
pts=zeros(length(xg(:))*2,1);
pts(1:2:end)=xg(:);
pts(2:2:end)=yg(:);
% time parameters for FSLE analysis
tspan=round([dayf day0]'*60*60*24);
Nstep=abs(diff(tspan))/(60*60*24)*4;

unfreezetime();

% Compute FSLE
% (functions defined in lamta.dev/lamta_all.cc)
out=FSLEext_v(tspan,pts,delta0,Nstep,delta);
%out=FSLEn_v(tspan,pts,3,delta0,Nstep,delta);

lambda=reshape(out(1,:),length(yvc),length(xvc));

function fname=all_lagrangian_diags(day0,product,out_path,ctday,numdays)
% Launches all the functions to compute the lagrangian diagnostics
% Modified from:
% $DATA/KEOPS2/octave/keops/analysis/realtime_cron/all_lagrangiandiags_keops2.m
%
	
% Load domain limits
% (domain_limits.m is a hard link to domain_limits.m in Mfiles directory)
%%AD: THE LINK IS NOT MORE NECESSARY, path added in MAIN_Lagrangian
domain_limits
lonv=Lon;
latv=Lat;
clear Lon Lat
% Load local vars
%lvars=localvars;

%--- Set parameters for Lagrangian analyses
delta0=0.01;
%delta0=.1; % Initial particle separation [degree]
delta=0.2;
%delta=.6;  % Final particle separation [degree]
ds=1/6; % Not used; What to do with this???
%AD: MOVED TO MAIN_Lagrangian Maximum number of days to compute particle trajectories
%numdays=120;
%numdays=10;


% parameter for OW dispersion analysis
%tollerance=10; CHANGED FOR TESTING, SINCE WE HAVE ONLY 3 FILES - ADP: remove next two lines and comment this one later
tollerance=1;
%numdays=1;

% lon lat matrices for Lagrangian analysis
lonvv=lonv(1):delta0:lonv(end);
latvv=latv(1):delta0:latv(end);
[long,latg]=meshgrid(lonvv,latvv);
sz=size(long);

%--- Begin Lagrangian analysis
 disp(['%   aviso_fsle'])
 [lambda,xvc,yvc,out]=aviso_fsle(day0,lonv,latv,delta0,delta,'none',numdays);
 theta=reshape(out(3,:),length(yvc),length(xvc));
 disp(' ')
 
% %---------------------------------------------------
 disp(['%   owdispersion 1/16=grid mesh '])
 OWdispersion_b=aviso_owdispersion(day0,lonv,latv,1/16,-numdays/2,tollerance,'none');

% Interpolate OW parameter to FSLE grid
 lonvo=lonv(1):1/16:lonv(end);
 latvo=latv(1):1/16:latv(end);
 [glono,glato]=meshgrid(lonvo,latvo);
 OWdispersion_bi=interp2(glono,glato,OWdispersion_b,long,latg,'cubic');

% disp(' ')
% %---------------------------------------------------
 
 disp(['%   aviso_ow '])
 owm=aviso_ow(day0,lonv,latv,1/16,'none');
 disp(' ')
 % Interpolate OW parameter to FSLE grid
 lonvo=lonv(1):1/16:lonv(end);
 latvo=latv(1):1/16:latv(end);
 [glono,glato]=meshgrid(lonvo,latvo);
 owmi=interp2(glono,glato,owm,long,latg,'cubic');

[U,V,xvc,yvc,Ucms,Vcms]=aviso_UV(day0,lonv,latv,delta0,'none');

 % Replace with MODIS SST?
 
 %disp(['          amsre_adv'])
 %[ssth,sstf,lonsstf,latsstf,sst]=amsre_adv(day0,lonv,latv,.1,-3,0,0,'none');
 
 disp(['%   aviso_bigsquadvec1'])
 [lon0,lat0,lonf15,latf15,sz15]=aviso_bigsqadvect(day0,lonv,latv,delta0,-15,'none');
 disp(' ')
  
%  disp(['%   aviso_bigsquadvec2'])
%  [lon0,lat0,lonf45,latf45,sz]=aviso_bigsqadvect(day0,lonv,latv,delta0,-45,'none');
  disp(' ')
  
% disp(['%   aviso_bigsquadvec3'])
% [lon0,lat0,lonf90,latf90,sz]=aviso_bigsqadvect(day0,lonv,latv,delta0,-90,'none');
 disp(' ')
  
  lon0=reshape(lon0,sz);
  lat0=reshape(lat0,sz);
  lonf=reshape(lonf15,sz);
  latf=reshape(latf15,sz);

%---------------------------------------------------
% possegday0b=find((OWdispersion_b(:)>20));
%  
%  disp(['%   aviso_lonlatadvect1'])
%  [lonfowadv15,latfowadv15,lonsmowadv15,latsmowadv15]=aviso_lonlatadvect(day0,long(possegday0b),latg(possegday0b),-15,'none');
%  disp(' ')
%  
%  disp(['%   aviso_lonlatadvect2'])
%  [lonfowadv60,latfowadv60,lonsmowadv60,latsmowadv60]=aviso_lonlatadvect(day0,long(possegday0b),latg(possegday0b),-60,'none');
%  disp(' ')
%---------------------------------------------------

% Replace with time from STRASSE domain?

%%% REMOVED FOR OSCAHR
%disp(['          timefrombathy'])
%touched=timefrombathy(day0,long,latg,-700,-numdays,'none');
%disp(['          '])

%--- Save workspace;
%    (loaded afterwards to plot figures)
todayv=datevec(now);
%fname=sprintf('%s/all_lagr_diag_%s_%d%02d%02d_%d%02d%02d.mat',out_path,product,day0(1),day0(2),day0(3),todayv(1),todayv(2),todayv(3));
fname=sprintf('%s/%d%02d%02d_%s_d%s',out_path,day0(1),day0(2),day0(3),product,num2str(ctday-1));

save('-mat-binary',[fname,'.mat']);
save('-mat-binary',[fname,'_lambda_only.mat'],'lonvv','latvv','lambda');
%end


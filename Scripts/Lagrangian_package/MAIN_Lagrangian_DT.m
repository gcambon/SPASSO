% Make a Lagrangian analysis using AVISO uv products
%	Plot: -FSLE
%	      -Okubo-Weiss parameter
%	      -Lon and Lat advection
%	      -Velocity fields
%==============================================================================
% Modified version of:
% /mnt/netapp-nencioli/DATA/KEOPS2/octave/keops/analysis/realtime_cron/an01_kerguelen_products.m
%
% With respect to KEOPS this file is directly called from bash as opposed to
% keops_realtime_cron.bsh;
% So I added the arguments below to replace:
% - cd to lagrangian mfile directory
% - call to init.m that: - added a bunch of directories to the path
%                        - called loclavars which
%                          defined lvars.datadir and lvars.lamtadir
%                          (lvars.lamtadir is never used afterwrds)
%
% DATADIR and LAMTADIR are declared global because they are afterwards used
% inside localvars.m calle by multiple scripts;
%
%2015/03/15: AD&FdO: added regional for med sea
%2015/03/15: AD: modifed to run in operational mode at delayed time
%AD:02/02/2015:inserted TEST/OPERATIONAL MODE
%AD:03/02/2015:used numdays from all_lagrangian_diags also for the number of aviso files to download
%==============================================================================

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' BEGIN OF ALL LAGRANGIAN DIAG')
disp(' ')

%-----------------------------------------------------------
% Cell array with arguments passed from bash script
args=argv;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(args) %%% TEST MODE %%
%%%
main_path='/home/doglioli/SPASSO/'
lagrang_path=[main_path,'Scripts/Lagrangian_package/'];
global DATADIR=[main_path,'Data/'];
global LAMTADIR=[lagrang_path,'lamta.dev/'];
global out_path=[main_path,'Cruises/PEACETIME/Wrk/'];
dayprod=datenum(2016,10,31);
products='ftp.aviso.oceanobs.com/regional-mediterranean/near-real-time/grids/madt/all-sat-merged/uv/';   % Added ADP 31/10/2016
%%
else %%% OPERATIONAL  MODE %%%
%%% 
 lagrang_path=[args{1},'/'];
 global DATADIR=[args{2},'/'];
 global LAMTADIR=[args{3},'/'];
 global out_path=[args{4},'/Wrk_DT/'];

 %dayprod=now
% if isempty(args{5})
%   dayprod=now
% else
   dayprod=datenum(args{5},'yyyymmdd')
   datestr(dayprod)
% end
 products=[args{6}];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To replace init.m
cd(lagrang_path)
lvars=localvars(out_path);
addpath(lvars.lamtadir)
addpath(args{4}) %AD: to charge domain_limits.m

%products={'nrt_med_2014'}; 
% products={'nrt_global_2014'}; % Only global products for OUTPACE
% products={'nrt_global_2010'}; % Only global products for OUTPACE

%days=[1 7]; % To process d0 and d6 data
days=[1]; % To process d0 data
% (add 4 to the array to process also d3)

%Maximum number of days to compute particle trajectories
%and for the number of dayly aviso files to download
numdays=40;


for ctproduct=1:size(products,1)

	product=products(ctproduct,:)
	datestr(dayprod-numdays)
	%disp(product)
	err=aviso_load_product(datevec(dayprod-numdays),datevec(floor(dayprod)),product,out_path);

	if(~err)
		for ctday=days
			%day0=datevec(dayprod-ctday);
			day0=datevec(dayprod);
			disp('#===============================================')
			disp(['# Processing ',product,' for d',num2str(ctday-1)])
			% Changed from all_lagrangiandiags_keops2
			fname=all_lagrangian_diags(day0,product,out_path,ctday,numdays);
                        %all_lagrangian_diags
			disp('%-----------------------------------------------')

			 disp(['# Figures for product ',product,' for d',num2str(ctday-1)])
			 %all_lagrangian_figs(fname,product,out_path,ctday);
			 disp('#===============================================')
			 disp(' ')
		end
	end

end

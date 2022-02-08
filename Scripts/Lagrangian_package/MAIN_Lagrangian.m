% Make a Lagrangian analysis using CMEMS uv products
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
%==============================================================================

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' BEGIN OF ALL LAGRANGIAN DIAG')
disp(' ')
%start
%-----------------------------------------------------------
% Cell array with arguments passed from bash script
argv='';
args=argv;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(args) %%% TEST MODE %%
    %%%
    %main_path='/home/doglioli/SPASSO/'
    %lagrang_path=[main_path,'Scripts/Lagrangian_package/'];
    main_path='/home/gcambon/HCONFIGS_SPASSO/DEMO2/'
    lagrang_path=[main_path,'/Lagrangian_package/'];
    global DATADIR=[main_path,'/Data/ALTI/nrt.cmems-du.eu/Core/SEALEVEL_EUR_PHY_L4_NRT_OBSERVATIONS_008_060/dataset-duacs-nrt-europe-merged-allsat-phy-l4/'];
    global LAMTADIR=[lagrang_path,'lamta.dev/'];
    global out_path=[main_path,'Wrk/'];
    dayprod=datenum(2022,02,04);
    %products='nrt_global_allsat';
    products='nrt_cmems';
    %%
else %%% OPERATIONAL  MODE %%%
    %%% 
    lagrang_path=[args{1},'/'];
    global DATADIR=[args{2},'/'];
    global LAMTADIR=[args{3},'/'];
    global out_path=[args{4},'/Wrk/'];
    dayprod=datenum(args{5},'yyyymmdd')
    datestr(dayprod)
    products=[args{6}];
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To replace init.m
cd(lagrang_path)
lvars=localvars(out_path);
addpath(lvars.lamtadir)
%addpath(args{4})

days=[1]; % To process d0 data

% Maximum number of days to compute particle trajectories
% and for the number of dayly CMEMS files to download
numdays=2; 
%numdays=60; 


for ctproduct=1:size(products,1)
	product=products(ctproduct,:)
	datestr(dayprod-numdays)	
	disp(product)

	err=aviso_load(datevec(dayprod-numdays),datevec(floor(dayprod)),product,out_path,out_path);
        
	if(~err)
	        for ctday=days
		    day0=datevec(dayprod);
			disp('#===============================================')
			disp(['# Processing ',product,' for d',num2str(ctday-1)])
			% Changed from all_lagrangiandiags_keops2
			fname=all_lagrangian_diags(day0,product,out_path,ctday,numdays);
                        % all_lagrangian_diags
			disp('%-----------------------------------------------')
                        
			 disp(['# Figures for product ',product,' for d',num2str(ctday-1)])
			 % all_lagrangian_figs(fname,product,out_path,ctday);
			 disp('#===============================================')
			 disp(' ')
		end
	end
        
end

function lvars=localvars(args)
% Previously defined in MAIN_Lagrangian.m
global DATADIR
global LAMTADIR
global WRKDIR


%lvars.datadir='/home/SPASSO/Data';
lvars.datadir=DATADIR; %LR: get datadir defined in SPASSO
lvars.lamtadir='/home/SPASSO/Scripts/Lagrangian_package/lamta.dev/';
lvars.wrkdir=sprintf(args);




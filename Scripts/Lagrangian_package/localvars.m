function lvars=localvars(args)
% Previously defined in MAIN_Lagrangian.m
global DATADIR
global LAMTADIR
global WRKDIR

lvars.datadir=DATADIR; 
lvars.lamtadir=LAMTADIR
lvars.wrkdir=sprintf(args);




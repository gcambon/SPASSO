function depth=oceandepth(lons,lats)

lvars=localvars;

bathyfile=sprintf('%s/bathy/Bathy180_30.mat',lvars.datadir);
load(bathyfile);
[xg,yg]=meshgrid(x1,y1);
depth=interp2(xg,yg,z1',lons,lats,'linear');





function err=aviso_load_cmems_product(day0,dayf,path_to_product,wrk_dir)

err=0;

if(length(day0)>1)
	day0=datenum(day0);
	end
if(length(dayf)>1)
	dayf=datenum(dayf);
	end

numdays=(dayf-day0)+1;

select(4); 

ctr=0;
for ctday=day0:1:dayf
dayvec=datevec(ctday);
dayvec=dayvec(1:3);
[Udegs,Vdegs,lonw,latw]=getncUV_aviso_cmems_product(path_to_product,dayvec,wrk_dir);
if(isempty(Udegs))
		disp(sprintf('Error in finding Aviso velocity files (day %d)',ctday));
		err=1;
		return;
	endif

if(ctr==0) %For the first frame to load, allocate the array and set the geometry
	sz=size(Udegs);
	
	%%NB!!!!!!!!!
	%For backward compatibility, in the line below we have to maintain the days with 0 at 1/1/1950,
	% otherwise all the functions aviso_<diagnostic> do not work.
	%This is a workaround which should be changed in a future version. 
	%%!!!!!!!!!!!

	set_par([lonw(1) lonw(end) sz(1) latw(1) latw(end) sz(2) (day0-datenum([1950 1 1]))*60*60*24 (dayf-datenum([1950 1 1]))*60*60*24 numdays]');

	% field_geometry(disttype,datatype,gridtype) defined in lamta.dev/lamta_all.cc
	% disttype=1 (Euclidean), 2 (sphere).
	% datatype=1 (deg./sec.), 2 (cm/sec.)
	% gridtype=0 (flat), 1 (sphere, regular), 2 (sphere Mercator)
	% field_geometry(2,1,2);
	field_geometry(2,1,1); % Now (2014) all grids for AVISO netcdf files are Cartesian instead of Mercator (ALL products) 
	endif

LUT_frame_fill(Udegs,Vdegs,ctr);
ctr=ctr+1;
endfor


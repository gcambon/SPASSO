function err=aviso_load(day0,dayf,product,path_to_product,wrk_dir)
% The function choses different functions to read AVISO velocity fields depending on the product specified;
%
% day0 and dayf define the day range of data to be read;
%2015/03/15: AD&FdO: added regional for med sea
%if(length(day0)>1)
%	day0=datenum(day0)-datenum([1950 1 1]);
%end
%if(length(dayf)>1)
%	dayf=datenum(dayf)-datenum([1950 1 1]);
%end

err=0;
switch product
	case 'upd_global_egm'
		aviso_upd_egm_load(floor(day0),ceil(dayf));
	case 'upd_global'
		aviso_upd_load(floor(day0),ceil(dayf));
	case 'upd_global_2010'
		aviso_upd_2010_load(floor(day0),ceil(dayf));
	case 'upd_global_largescale'
		aviso_upd_largescale_load(floor(day0),ceil(dayf));
	case 'upd_global_largescale_mdt'
		aviso_upd_largescale_mdt_load(floor(day0),ceil(dayf));
	case 'upd_global_sa'
		aviso_sa_upd_load(floor(day0),ceil(dayf));
	case 'upd_global_sa_tpj1'
		aviso_dt_ref_global_tpj1_msla(floor(day0),ceil(dayf));
	case 'upd_global_ssh_tpj1'
		aviso_dt_ref_global_tpj1_ssh(floor(day0),ceil(dayf));
	case 'global_mdt'
		aviso_mdt_load(floor(dayf-8),ceil(dayf));
	case 'nrt_global'
		disp('nrt_global obsolete, using nrt_global_2010 instead.');
		err=aviso_nrt_2010_load(floor(day0),ceil(dayf));
	case 'nrt_global_2010'
		err=aviso_nrt_2010_load(floor(day0),ceil(dayf));
        case 'nrt_global_2014'
                err=aviso_nrt_2014_load(floor(day0),ceil(dayf));
	case 'upd_med'
		avisomed_upd_load(floor(day0),ceil(dayf));
	case 'upd_med_sa'
		avisomed_sa_upd_load(floor(day0),ceil(dayf));
	case 'nrt_med'
		err=avisomed_nrt_load(floor(day0),ceil(dayf));
	case 'nrt_med_2014'
		err=avisomed_nrt_2014_load(floor(day0),ceil(dayf));
	case 'rt_med'
		err=avisomed_rt_load(floor(day0),ceil(dayf));	
	case 'nrt_med_sa'
		err=avisomed_nrt_load(floor(day0),ceil(dayf),'ssa',0);
	case 'rt_med_sa'
		err=avisomed_nrt_load(floor(day0),ceil(dayf),'ssa',1);
	case 'rt_global'
		disp('rt_global obsolete, using nrt_global_2010 instead.');
		err=aviso_nrt_2010_load(floor(day0),ceil(dayf));	
	case 'mdt_global'
		err=aviso_mdt_load(floor(day0),ceil(dayf));
	case 'upd_kerguelen'
		aviso_upd_kerguelen_load(floor(day0),ceil(dayf));
	case 'nrt_kerguelen'
		err=aviso_nrt_kerguelen_load(floor(day0),ceil(dayf));
	case 'nrt_kerguelen_ekman'
		err=aviso_nrt_kerguelen_load(floor(day0),ceil(dayf),1);
	case 'dt_product'
		err=aviso_load_product(floor(day0),ceil(dayf),path_to_product,wrk_dir);
	case 'nrt_product'
		err=aviso_load_product(floor(day0),ceil(dayf),path_to_product,wrk_dir);
	case 'dt_coralsea'
		err=aviso_load_coralsea_product(floor(day0),ceil(dayf));
	case 'dt_cmems'
		err=aviso_load_cmems_product(floor(day0),ceil(dayf),path_to_product,wrk_dir);
        case 'nrt_cmems'
                err=aviso_load_cmems_product(floor(day0),ceil(dayf),path_to_product,wrk_dir);
        case 'nrt_cmems_med'
                err=aviso_load_cmems_product_med(floor(day0),ceil(dayf),path_to_product,wrk_dir); 
	case 'none'
		disp('product=none: Not loading any new AVISO velocity');
	otherwise
		error('Invalid aviso product');
end


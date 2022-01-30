function stations(mycolor)
%%% ADD CRUISE STATIONS
SD_stat=load('../../Cruises/PEACETIME/station_coord.txt');

lon_stat_SD=SD_stat(:,2);
lat_stat_SD=SD_stat(:,1);

plot(lon_stat_SD,lat_stat_SD,'-o','Markersize',6,'Color',mycolor,'linewidth',4);


return


%- Stations positions -% (degres east)
%lat_stat = -19;
%lon_stat_A = 155;
%lon_stat_B = 176.32;
%lon_stat_C = 200; %(°E)
%SD_stat=load('SD_stations_coord.txt');
%lon_stat_SD=SD_stat(2:18,2);
%lat_stat_SD=SD_stat(2:18,1);
%
%- Plot stations -%
%plot(lon_stat_A,lat_stat,'+','Color',mycolor,'Markersize',10,'linewidth',5);
%text(lon_stat_A,lat_stat+1.2,'A','Color',mycolor);
%plot(lon_stat_B,lat_stat,'+','Color',mycolor,'Markersize',10,'linewidth',5);
%text(lon_stat_B,lat_stat+1.2,'B','Color',mycolor);
%plot(lon_stat_C,lat_stat,'+','Color',mycolor,'Markersize',10,'linewidth',5);
%text(lon_stat_C,lat_stat+1.2,'C','Color',mycolor);
%plot(lon_stat_SD,lat_stat_SD,'+','Markersize',3,'Color',mycolor,'linewidth',4);
%text(lon_stat_SD(1:8),lat_stat_SD(1:8)+0.9,num2str([2:9]'),'Color',mycolor);
%text(lon_stat_SD(9:end)-0.5,lat_stat_SD(9:end)+0.9,num2str([11:19]'),'Color',mycolor);

%- Plot territorial waters -%

%load ZEE_NewCaledonia.txt
%plot([ZEE_NewCaledonia(:,2)],[ZEE_NewCaledonia(:,1)],'-','Color',mycolor,'Linewidth',2)
%text(156,-16,'NOUVELLE','Color',mycolor)%,'FontWeight','bold')
%text(156,-17,'CALEDONIE','Color',mycolor)%,'FontWeight','bold')

%- Plot VANUATU territorial water limits -%
%load ZEE_Vanuatu.txt
%plot([ZEE_Vanuatu(:,2)],[ZEE_Vanuatu(:,1)],'-','Color',mycolor,'Linewidth',2)
%plot([163.34 170.54 173.26 173.55 172.08 171.61 167.78 164.15 163.34],[-14.77 -21.64 -19.89 -18.53 -16.5 -14.96 -12.32 -12.65  -14.77],'-','Color',mycolor,'Linewidth',2)
%plot([174.5 174.5],[-19.7 -18.3],'-','Color',mycolor,'Linewidth',5)
%text(167,-17,'VANUATU','Color',mycolor)%,'FontWeight','bold')

%- Plot TONGA territorial water limits -%
%plot([189 189],[-19.7 -18.3],'-','Color',mycolor,'Linewidth',5)
%plot([194 194],[-19.7 -18.3],'-','Color',mycolor,'Linewidth',5)
%text(188.5,-21,'TONGA','Color',mycolor)%,'FontWeight','bold')

%- Plot ILES COOK territorial water limits -%
%plot([199 199],[-19.7 -18.3],'-','Color',mycolor,'Linewidth',5)
%plot([209 209],[-19.7 -18.3],'-','Color',mycolor,'Linewidth',5)
%text(201,-21,'ILES','Color',mycolor)%,'FontWeight','bold')
%text(201,-22,'COOK','Color',mycolor)%,'FontWeight','bold')

%- Plot NIUE territorial water limits -%
%plot([194 194],[-19.7 -18.3],'-','Color',mycolor,'Linewidth',5)
%plot([199 199],[-19.7 -18.3],'-','Color',mycolor,'Linewidth',5)
%text(194.5,-17,'NIUE','Color',mycolor)%,'FontWeight','bold')

%- Plot FIDJI territorial water limits -%
%plot([174.5 174.5],[-19.7 -18.3],'-','Color',mycolor,'Linewidth',5)
%plot([189 189],[-19.7 -18.3],'-','Color',mycolor,'Linewidth',5)
%text(180.5,-21,'FIDJI','Color',mycolor)%,'FontWeight','bold')


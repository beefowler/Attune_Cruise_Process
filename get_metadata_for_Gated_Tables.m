
%% Manually choose cruise to process
basepath = '\\sosiknas1\Lab_data\Attune\cruise_data\20201013_EN657\preserved';
restpath =  'https://nes-lter-data.whoi.edu/api/ctd/en657/';

cruisename = 'EN657';
elogpath = '\\sosiknas1\Lab_data\LTER\20201013_EN657\eLog\R2R_ELOG_EN657_FINAL_EVENTLOG_20201018_134037.csv'; %set to '' if there are no discrete underway samples 
uw_fullname = 'https://nes-lter-data.whoi.edu/api/underway/en657.csv';

%% some file structure setup
fpath = [basepath filesep 'FCS' filesep];
outpath = [basepath filesep 'outputs' filesep];
classpath = [outpath 'class' filesep];

if ~exist(outpath, 'dir')
    mkdir(outpath)
end
if ~exist(classpath, 'dir')
    mkdir(classpath)
end

T = load([outpath '\FCSlist.mat']);
T = T.FCSList; 

G = load([outpath '\Gated_Table.mat']);
gated_table = G.gated_table; 


if startsWith(restpath, 'https')
bottledata = webread([restpath 'bottles.csv']); 
metadata = webread([restpath 'metadata.csv']); 

castlist = unique(gated_table.Cast);

for f = 1:height(gated_table)

    if gated_table.Cast(f) == 0 
        continue
    end

    %first get cast metdata
    gated_table.Latitude(f) = metadata.latitude(metadata.cast == gated_table.Cast(f));
    gated_table.Longitude(f) = metadata.longitude(metadata.cast == gated_table.Cast(f));
    gated_table.date_sampled(f) = metadata.date(metadata.cast == gated_table.Cast(f));
    gated_table.nearest_station(f) = metadata.nearest_station(metadata.cast == gated_table.Cast(f));

    %now add depths for each niskin number
    gated_table.salinity(f) = bottledata.sal00(bottledata.cast == gated_table.Cast(f) & bottledata.niskin == gated_table.Niskin(f));
    gated_table.potemp090c(f) = bottledata.potemp090c(bottledata.cast == gated_table.Cast(f) & bottledata.niskin == gated_table.Niskin(f));
    gated_table.depth_m(f) = bottledata.depsm(bottledata.cast == gated_table.Cast(f) & bottledata.niskin == gated_table.Niskin(f));

    
end
    
    gated_table.Longitude(gated_table.Latitude == 0) = NaN;
    gated_table.salinity(gated_table.Latitude == 0) = NaN;
    gated_table.potemp090c(gated_table.Latitude == 0) = NaN;
    gated_table.depth_m(gated_table.Latitude == 0) = NaN;
    gated_table.Latitude(gated_table.Latitude == 0) = NaN;

else %use mat file for SPIROPA Cruises not the same format >:(

load(restpath)
castlist = unique(BTL.cast);
FilesToADD_w_meta = table(); 

for c = 1:length(castlist)

end


end


save([outpath '\Gated_Table.mat'], 'gated_table', '-append'); 

%% now underways

if ~isempty(elogpath)
load([outpath '\Gated_Table.mat'])

A = readtable(elogpath); 
temp = A((contains(A.Comment, 'FCM')), :) ;

%first pick times that match 
ind = find(gated_table.Cast ==0); 
temp_g = gated_table(ind, :) ; 
for s = 1:height(temp_g)
    disp(temp)
    pause
    UWnum = input([temp_g.fcslist{s} '-Which elog row should we use for this fcs file?']);
    stupiddate = temp.dateTime8601(UWnum); 
    stupiddate = replace(stupiddate, 'T', ' ');
        stupiddate = replace(stupiddate, '0000', '00:00');
    gated_table.date_sampled{ind(s)} = stupiddate;
    datenums(s) = datenum(temp.GPS_Time(UWnum)); 
end

%%
%then go get underway data 

udatenums = unique(datenums); %only need one underway match per UW sample

if strncmp ('http', uw_fullname, 4) %case for API 
    uw = webread(uw_fullname);
    dt = datetime(uw.date, 'InputFormat', 'yyyy-MM-dd HH:mm:ss+00:00');
    uw_mdate = datenum(dt);
else % works for R/V Endeavor cruise files
    fid = fopen(uw_fullname);
    t = fgetl(fid); fclose(fid);
    numHdrLines = str2num(regexprep(t, '#DataStartLine:',''));
    opts = delimitedTextImportOptions('VariableNamesLine', numHdrLines, 'DataLines', numHdrLines+1);
    uw = readtable(uw_fullname, opts);    
    dt = datetime(uw.DateTime_ISO8601, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSS''Z');
    uw_mdate = datenum(dt);
% new case needed for Armstrong--easiest to use API if possible!
end

tdiff = NaN(length(udatenums),1);
match_ind = tdiff;
for ii = 1:length(tdiff)
    [tdiff(ii), match_ind(ii)] = min(abs(udatenums(ii)-uw_mdate));

    %now pick variables to save 

    %env data parsing, my favorite! -- modified from 
    %get_crusie_voldists_fromEDItable2.m

    %Latitude & Longitude 
    if startsWith(cruisename, 'EN') 
        Exist_Column = strcmp('gps_furuno_latitude', uw.Properties.VariableNames); %use gps-furuno if quality flag is not 0
        Exist_Column2 = strcmp('gps_garmin741_latitude', uw.Properties.VariableNames); %use gps-furuno if quality flag is not 0

        if Exist_Column(Exist_Column ==1)
            lat = uw.gps_furuno_latitude;
            lon = uw.gps_furuno_longitude;
            lat(uw.gps_furuno_quality == 0) = NaN;
            lon(uw.gps_furuno_quality == 0) = NaN;
        elseif Exist_Column2(Exist_Column2 ==1) %super early cruises dont have gps-furuno
            lat = uw.gps_garmin741_latitude;
            lon = uw.gps_garmin741_longitude;
            lat(uw.gps_garmin741_quality == 0) = NaN;
            lon(uw.gps_garmin741_quality == 0) = NaN;
        end
       
    elseif startsWith(cruisename, 'AR')
        Exist_Column4 = strcmp('dec_lat',uw.Properties.VariableNames); 
        Exist_Column5 = strcmp('Dec_LAT', uw.Properties.VariableNames); 
        if Exist_Column4(Exist_Column4==1) 
            lat = uw.dec_lat; 
            lon = uw.dec_lon; 
        elseif Exist_Column5(Exist_Column5==1) 
            lat = uw.Dec_LAT; 
            lon = uw.Dec_LON;
        end
    else %other cruises
        Exist_Column1 = strcmp('gps_furuno_latitude', uw.Properties.VariableNames); %use gps-furuno if quality flag is not 0
        Exist_Column2 = strcmp('dec_lat',uw.Properties.VariableNames); 
        Exist_Column3 = strcmp('Dec_LAT', uw.Properties.VariableNames); 
        Exist_Column4 = strcmp('lat_flr', uw.Properties.VariableNames); 
        Exist_Column5 = strcmp('lat_SAMOS', uw.Properties.VariableNames);
        Exist_Column6 = strcmp('lat_tsg',uw.Properties.VariableNames); 
        Exist_Column7 = strcmp('lat',uw.Properties.VariableNames); 
        if Exist_Column1(Exist_Column1==1) 
           lat = uw.gps_furuno_latitude;
           lon = uw.gps_furuno_longitude;
        elseif Exist_Column2(Exist_Column2==1) 
            lat = uw.dec_lat; 
            lon = uw.dec_lon; 
        elseif Exist_Column3(Exist_Column3==1) 
            lat = uw.Dec_LAT; 
            lon = uw.Dec_LON;
        elseif Exist_Column4(Exist_Column4==1) 
            lat = uw.lat_flr; 
            lon = uw.lon_flr; 
        elseif Exist_Column5(Exist_Column5==1) 
            lat = uw.lat_SAMOS; 
            lon = uw.lon_SAMOS;
        elseif Exist_Column6(Exist_Column6==1) 
            lat = uw.lat_tsg; 
            lon = uw.lon_tsg; 
        elseif Exist_Column7(Exist_Column7==1) 
            lat = uw.lat; 
            lon = uw.lon; 
        end
    end
   

    %Temperature & Salinity 
    
     if startsWith(cruisename, 'EN')

         if strcmp(cruisename, 'EN627') %SUPER ODDBALL We use tsg2_temperature
              temperature = uw.tsg2_temperature; 
         else
            temperature = uw.tsg1_sst;
         end

         if sum(strcmp(cruisename, {'EN627'; 'EN644'; 'EN668'})) 
            salinity = uw.tsg2_salinity; 
         else
             salinity = uw.tsg1_salinity; 
         end

     elseif startsWith(cruisename, 'AR')
          Exist_Column5 = strcmp('sbe45s', uw.Properties.VariableNames);
          Exist_Column6 = strcmp('SBE45S', uw.Properties.VariableNames);

          if Exist_Column5(Exist_Column5==1)
                salinity = uw.sbe45s; 
                temperature = uw.sbe48t; 
          elseif Exist_Column6(Exist_Column6==1)
                salinity = uw.SBE45S; 
                temperature = uw.SBE48T; 
          end
     else %other cruise vessels 
        Exist_Column = strcmp('sbe45S',uw.Properties.VariableNames);  
        Exist_Column2 = strcmp('tsg1_salinity', uw.Properties.VariableNames);  %sometimes we want tsg2, depends on cruise
        Exist_Column3 = strcmp('salinity', uw.Properties.VariableNames);
        Exist_Column4 = strcmp('s', uw.Properties.VariableNames);
        Exist_Column5 = strcmp('sbe45s', uw.Properties.VariableNames);
        Exist_Column6 = strcmp('SBE45S', uw.Properties.VariableNames);
        if Exist_Column(Exist_Column==1) 
            salinity = uw.sbe45S; 
            temperature = uw.sbe45T; 
        elseif Exist_Column2(Exist_Column2==1)
            salinity = uw.tsg1_salinity; 
            temperature = uw.tsg1_sst; %this is good. We don't want tsg1_temperature. Even if we use tsg2 for salinity, we want tsg1_sst for temp. 
        elseif Exist_Column3(Exist_Column3==1)
            salinity = uw.salinity; 
            temperature = uw.temperature; 
        elseif Exist_Column4(Exist_Column4==1)
            salinity = uw.s; 
            temperature = uw.t1; 
        elseif Exist_Column5(Exist_Column5==1)
            salinity = uw.sbe45s; 
            temperature = uw.sbe48t; 
        elseif Exist_Column6(Exist_Column6==1)
            salinity = uw.SBE45S; 
            temperature = uw.SBE48T; 
        end
     end
   
    
    gated_table.Latitude(ind(datenums == udatenums(ii))) = lat(match_ind(ii)); 
    gated_table.Longitude(ind(datenums == udatenums(ii))) = lon(match_ind(ii)); 
    gated_table.salinity(ind(datenums == udatenums(ii))) = salinity(match_ind(ii)); 
    gated_table.potemp090c(ind(datenums == udatenums(ii))) = temperature(match_ind(ii)); 

    gated_table.depth_m(ind(datenums == udatenums(ii))) = 5; 

end

save([outpath '\Gated_Table.mat'], 'gated_table', '-append'); 


end
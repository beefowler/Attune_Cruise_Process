% want to convert the table that has rows for each file into a table that
% has a single row for each niskin / cast

%the challenge here will be looking to the file which has the most trusted
%counts for each group. 

%this will also be the step that removes prochlorococcus counts from the
%bacteria counts when applicable 

%% Manually choose cruise to process
basepath = '\\sosiknas1\Lab_data\Attune\cruise_data\20210203_EN661\preserved\';



%% some file structure setup
fpath = [basepath filesep 'FCS' filesep];
outpath = [basepath filesep 'outputs' filesep];

G = load([outpath '\Gated_Table.mat']);
gated_table = G.gated_table; 


%first some counting of underway samples
    temp = gated_table(gated_table.Cast == 0, :);
    uw_list = unique(string(temp.date_sampled)) ;
    uw_num = zeros(height(gated_table), 1); 
    for i = 1:length(uw_list)
        uw_num(gated_table.Cast == 0 & strcmp(string(gated_table.date_sampled), uw_list(i))) = i; 
    end


[G, C, N, uw] = findgroups(gated_table.Cast, gated_table.Niskin, uw_num);
ia = 1:max(G); 

CNTable = table; 
CNTable.Cast = C; 
CNTable.Niskin = N; 

%preallocate columns so we can get nans rather than zeros
hetbaccol = nan(height(CNTable), 1); 
syncol = hetbaccol; 
procol = syncol; 
eukcol = syncol; 
lp_eukcol = syncol; 
hp_eukcol = syncol; 
vols = nan(height(CNTable), 6); 

for g = 1:max(G)
    

    temp = gated_table(G == g, :);

    CNTable.latitude(g) = temp.Latitude(1); 
    CNTable.longitude(g) = temp.Longitude(1);
    CNTable.nearest_station(g) = temp.nearest_station(1);
    CNTable.salinity(g) = temp.salinity(1);
    CNTable.potemp090c(g) = temp.potemp090c(1);
    CNTable.depth_m(g) = temp.depth_m(1);
    CNTable.date_sampled(g) = temp.date_sampled(1); 

 
    %pick file to count Syn
    ind = find(contains(temp.fcslist, 'phyto_PE'));
    use_euk_for_syn = 0; 
    if length(ind) == 1 %if only one fcs file of this type, use that. 
        CNTable.Synfile(g) = temp.fcslist(ind); 
        syncol(g) = temp.Syn_conc(ind); 
        vols(g, 2) = temp.median_volumes_syn(ind); 
    elseif ~isempty(ind) %if more than 1, get most recent
        timesince = []; 
        for f = 1:length(ind)
            time = dir(fullfile(fpath, temp.fcslist{ind(f)}));
            time = datetime(time.date);  
            timesince = [timesince datetime()-time];
        end
        [~,truind] = min(timesince);
        CNTable.Synfile(g) = temp.fcslist(ind(truind));
        syncol(g) = temp.Syn_conc(ind(truind)); 
        vols(g, 2) = temp.median_volumes_syn(ind(truind)); 

    elseif isempty(ind) %no special PE runs
        use_euk_for_syn = 1; 
    end


    %now pick file to count Euks
    ind = find(contains(temp.fcslist, 'phyto_CHL') & ~contains(temp.fcslist, 'pro', 'IgnoreCase', true));
    if length(ind) == 1 %if only one fcs file of this type, use that. 
         if use_euk_for_syn
          CNTable.Synfile(g) = temp.fcslist(ind); 
          syncol(g) = temp.Syn_conc(ind); 
          vols(g, 2) = temp.median_volumes_syn(ind); 

         end
          CNTable.Eukfile(g) = temp.fcslist(ind); 
          eukcol(g) = temp.Euk_conc(ind); 
          lp_eukcol(g) = temp.LowP_Euk_conc(ind); 
          hp_eukcol(g) = temp.HighP_Euk_conc(ind); 
          vols(g, 1) = temp.median_volumes_euk(ind); 
          vols(g, 5) = temp.median_volumes_low_pe_euk(ind); 
          vols(g, 6) = temp.median_volumes_high_pe_euk(ind); 
    elseif ~isempty(ind) %if more than 1, get most recent
        timesince = []; 
        for f = 1:length(ind)
            time = dir(fullfile(fpath, temp.fcslist{ind(f)}));
            time = datetime(time.date);  
            timesince = [timesince datetime()-time];
        end
        [~,truind] = min(timesince);
        if use_euk_for_syn
            CNTable.Synfile(g) = temp.fcslist(ind(truind)) ; 
            syncol(g) = temp.Syn_conc(ind(truind)); 
            vols(g, 2) = temp.median_volumes_syn(ind(truind)); 
        end
        CNTable.Eukfile(g) = temp.fcslist(ind(truind)); 
        eukcol(g) = temp.Euk_conc(ind(truind)); 
        lp_eukcol(g) = temp.LowP_Euk_conc(ind(truind)); 
        hp_eukcol(g) = temp.HighP_Euk_conc(ind(truind)); 
        vols(g, 1) = temp.median_volumes_euk(ind(truind)); 
        vols(g, 5) = temp.median_volumes_low_pe_euk(ind(truind)); 
        vols(g, 6) = temp.median_volumes_high_pe_euk(ind(truind)); 

    end

    %pick file to count Prochlorococcus
    ind = find(contains(temp.fcslist, 'pro', 'IgnoreCase', true));
    noprofile = 1; 
    if isempty(ind)
        CNTable.ProFile{g} = '';
        noprofile = 1; 
    elseif length(ind) == 1 %if only one fcs file of this type, use that. 
        CNTable.ProFile(g) =  temp.fcslist(ind);
        procol(g) = temp.Pro_conc(ind); 
        vols(g, 4) = temp.median_volumes_pro(ind); 
         if isnan(procol(g)) %if there was no pro gate, make conc zero 
            procol(g) = 0; 
        end
    else %if more than 1, get most recent
        timesince = []; 
        for f = 1:length(ind)
            time = dir(fullfile(fpath, temp.fcslist{ind(f)}));
            time = datetime(time.date); 
            timesince = [timesince datetime()-time];
        end
        [~,truind] = min(timesince);
        CNTable.ProFile(g) =  temp.fcslist(ind(truind));
        procol(g) = temp.Pro_conc(ind(truind)); 
        vols(g, 4) = temp.median_volumes_pro(ind(truind)); 
        if isnan(procol(g)) %if there was no pro gate, make conc zero 
            procol(g) = 0; 
        end
    end




    %pick file to count bacteria, and subtract prochlorococs
    %if a pro file exists and there is no pro gate in bacteria file, then
    %subtract concentration 
    ind = find(contains(temp.fcslist, 'hbac'));
    if length(ind) == 1 %if only one fcs file of this type, use that. 
        CNTable.BacteriaFile(g) =  temp.fcslist(ind);

       if ~noprofile && isnan(temp.Pro_conc(ind)) %pro run, but no pro gate in bacteria run 
           hetbaccol(g) = temp.HetBact_conc(ind) - procol(g); 
       else
           hetbaccol(g) = temp.HetBact_conc(ind); 
       end

       vols(g, 3) = temp.median_volumes_hetbact(ind); 

    elseif ~isempty(ind) %if more than 1, get most recent
        timesince = []; 
        for f = 1:length(ind)
            time = dir(fullfile(fpath, temp.fcslist{ind(f)}));
            time = datetime(time.date); 
            timesince = [timesince datetime()-time];
        end
        [~,truind] = min(timesince);
        CNTable.BacteriaFile(g) =  temp.fcslist(ind(truind));
        
        if ~noprofile && isnan(temp.Pro_conc(ind(truind))) %pro run, but no pro gate in bacteria run 
            hetbaccol(g) = temp.HetBact_conc(ind(truind)) - procol(g); 
        else
            hetbaccol(g) = temp.HetBact_conc(ind(truind)); 
        end
        vols(g, 3) = temp.median_volumes_hetbact(ind(truind)); 

    end


    end


CNTable.euk_per_ml = eukcol; 
CNTable.syn_per_ml = syncol; 
CNTable.pro_per_ml = procol; 
CNTable.hetbac_per_ml = hetbaccol; 

CNTable.low_pe_euk_per_ml = lp_eukcol; 
CNTable.high_pe_euk_per_ml = hp_eukcol; 


%reorder so it isn't so hard to see
CNTable = CNTable(:, [1:10 15 11 14 18 19 12 16 13 17]); 

   CNTable.median_volumes_euk = vols(:,1);
    CNTable.median_volumes_syn = vols(:,2);
    CNTable.median_volumes_hetbact = vols(:,3);
    CNTable.median_volumes_pro = vols(:,4);
    CNTable.median_volumes_low_pe_euk = vols(:,5);
    CNTable.median_volumes_high_pe_euk = vols(:,6);


save([outpath '\SummaryTable.mat'], 'CNTable')






% Working from Process_Preserved_Sampels_with_AWS
% want to actually classify individual particles in FCS files 
% this is the subjective part of the 


%% Manually choose cruise to process
basepath = '\\sosiknas1\Lab_data\Attune\cruise_data\20201013_EN657\preserved';

%% some file structure setup
fpath = [basepath filesep 'FCS' filesep];
outpath = [basepath filesep 'outputs' filesep];
classpath = [outpath 'class' filesep];
awspath = [basepath filesep 'aws\']; 

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



%% Here is where we decide which gates we are interested in and how we will find them

    % gates_of_interest = {'syn'; 'Syn'; 'Euks'; 'euks'; 'pro'; 'Pro'; 'bacteria'; 'Bacteria'}; 

    classnames = {'Euks == 1'; 'Syn == 2'; 'Bacteria == 3'; 'Pro == 4'; 'LowPE_Euks = 5'; 'HighPE_Euks - 6'};

    class = cell(1,height(gated_table)); 
    syn_conc = nan(height(gated_table),1); 
    euk_conc = syn_conc; 
    bact_conc = syn_conc; 
    pro_conc = syn_conc; 
    lp_euk_conc = syn_conc; 
    hp_euk_conc = syn_conc; 

for i = 1:height(gated_table) 
    filename = gated_table.fcslist{i};
    [fcsdat, fcshdr]  = fca_readfcs([fpath filename]);

    a = size(gated_table.gate_assignments{i}); 
 
    if max(a) == 0
        class{i} = []; 
        concentrations{i} = []; 
        continue 
    end

    %first go through and generate a single class assingment for each
    %particle
    gate_assign_i =gated_table.gate_assignments{i}; 
    class_i = zeros(1, max(size(gate_assign_i))); 
    gate_names = gated_table.gate_names{i}; 


    % look for euks
    if sum(strcmp(gate_names, 'euk'))
            gate_num_1 = find(strcmp(gate_names, 'euk')); 
            class_i(class_i == 0' & gate_assign_i(gate_num_1, :)) = 1; 
    end
    if sum(strcmp(gate_names, 'Euk'))
            gate_num_1 = find(strcmp(gate_names, 'Euk')); 
            class_i(class_i == 0' & gate_assign_i(gate_num_1, :)) = 1; 
    end

    if sum(strcmp(gate_names, 'LowP_Euk'))
            gate_num_5 = find(strcmp(gate_names, 'LowPE_Euk')); 
            class_i(class_i == 0' & gate_assign_i(gate_num_5, :)) = 5; 
    end
    if sum(strcmp(gate_names, 'HighP_Euk'))
            gate_num_6 = find(strcmp(gate_names, 'LowPE_Euk')); 
            class_i(class_i == 0' & gate_assign_i(gate_num_6, :)) = 6; 
    end

    %next syn 
    if sum(strcmp(gate_names, 'Syn'))
            gate_num_2 = find(strcmp(gate_names, 'Syn')); 
            class_i(class_i == 0' & gate_assign_i(gate_num_2, :)) = 2; 
    end
    if sum(strcmp(gate_names, 'syn'))
            gate_num_2 = find(strcmp(gate_names, 'syn')); 
            class_i(class_i == 0' & gate_assign_i(gate_num_2, :)) = 2; 
    end



    %next bacteria 
    if sum(strcmp(gate_names, 'Bacteria'))
            gate_num_b = find(strcmp(gate_names, 'Bacteria')); 
            class_i(class_i == 0' & gate_assign_i(gate_num_b, :)) = 3; 
    end
    if sum(strcmp(gate_names, 'bacteria'))
            gate_num_b = find(strcmp(gate_names, 'bacteria')); 
            class_i(class_i == 0' & gate_assign_i(gate_num_b, :)) = 3; 
    end



    %next pro 
    if sum(strcmp(gate_names, 'pro'))
            gate_num_4 = find(strcmp(gate_names, 'pro')); 
            class_i(class_i == 0' & gate_assign_i(gate_num_4, :)) = 4; 
    end
    if sum(strcmp(gate_names, 'Pro'))
            gate_num_4 = find(strcmp(gate_names, 'Pro')); 
            class_i(class_i == 0' & gate_assign_i(gate_num_43, :)) = 4; 
    end
    if sum(strcmp(gate_names, 'Euk_sm'))
            gate_num_4 = find(strcmp(gate_names, 'Euk_sm')); 
            class_i(class_i == 0' & gate_assign_i(gate_num_4, :)) = 4; 
    end



    %ok now save class assignments 

    class{i} = class_i;

    
    %onto calculating concentrations
    concent_i = nan(1, 6);

    if exist('gate_num_1')
        concent_i(1) = sum(class_i == 1)./gated_table.Vol_ml(i); 
    end

    if exist('gate_num_2')
        concent_i(2) = sum(class_i == 2)./gated_table.Vol_ml(i); 
    end

    if exist('gate_num_4')
        concent_i(4) = sum(class_i == 4)./gated_table.Vol_ml(i); 
    end

    if exist('gate_num_5')
        concent_i(4) = sum(class_i == 5)./gated_table.Vol_ml(i); 
    end
    if exist('gate_num_6')
        concent_i(4) = sum(class_i == 6)./gated_table.Vol_ml(i); 
    end

    %account for bacteria time gate
    g = gated_table.gate_logic{i}; 
    if exist('gate_num_b') & contains(g{gate_num_b}, 'time_to_include')
        p = gated_table.polygon_vars{i}; 
        emptyCells = cellfun(@isempty,p) ;
        [a,b] = find(emptyCells); 
        p{a,b} = ''; 
        [a,b] = find(cellfun(@(x) contains(x, 'Time'), p));
        pv = gated_table.polygon_vals{i} ; 
        vals = pv{b};
        t_val = vals(:,1); %this is the actual value range used for the Time variable
        
        
        matdate1 = [fcshdr.date ' ' fcshdr.starttime];
        matdate1 = datenum(matdate1, 'dd-mmm-yyyy HH:MM:SS');
        matdate2 = [fcshdr.date ' ' fcshdr.stoptime];
        matdate2 = datenum(matdate2, 'dd-mmm-yyyy HH:MM:SS');

        fulltime = 1000*etime(datevec(matdate2), datevec(matdate1)); %this should be 60 seconds * 1000 but just to make sure
        
        fraction_of_time = (t_val(2) - t_val(1))./fulltime;

        volume_considered = fraction_of_time.*gated_table.Vol_ml(i); 
        concent_i(3) = sum(class_i == 3)./volume_considered; 
    end

    
    disp(gated_table(i,:))
    disp(concent_i)

    %now save concentrations
    euk_conc(i) = concent_i(1);
    syn_conc(i) = concent_i(2);
    bact_conc(i) = concent_i(3) ;
    pro_conc(i) = concent_i(4);
    lp_euk_conc(i) = concent_i(5); 
    hp_euk_conc(i) = concent_i(6); 
   
        clear gate_num_1 gate_num_2 gate_num_4 gate_num_b gate_assign_i class_i gate_num_5 gate_num_6
end

gated_table.class = class'; 
gated_table.Euk_conc = euk_conc; 
gated_table.Low_P_Euk_conc = lp_euk_conc; 
gated_table.HighP_Euk_conc = hp_euk_conc; 

gated_table.Syn_conc = syn_conc; 
gated_table.HetBact_conc = bact_conc; 
gated_table.Pro_conc = pro_conc; 

save([outpath '\Gated_Table.mat'], 'gated_table', '-append'); 




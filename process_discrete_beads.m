
%process discrete bead runs and create a beadstat table. 
%This table will span all lab runs, so we will need to pick which files to use from
%timestamps. 


%This script assumes that all beadfiles are run with OD2 filter on GL1 as
%in more recent (2022) cruises. 


beadpath = '\\sosiknas1\Lab_data\Attune\cruise_data\beads\FCB_bead_mix_experiment_settings\between_cruises\'; 
beadfigpath = [beadpath filesep 'outputs\bead_plots'];

beadlist = dir(fullfile(beadpath, '*fcs')); 



%initialize bead table
beadstat = table;


%run through bead files in FCS directory
for b = 1:length(beadlist)
    clear temp_table
    disp(b)
    
    [fcsdat, fcshdr] = fca_readfcs([beadpath filesep beadlist(b).name]);
    
    ssca_chn = strmatch('SSC-A', {fcshdr.par.name});
    gl1a_chn = strmatch('GL1-A', {fcshdr.par.name});

    bead_ch_names{1} = 'GL1-A';
    bead_ch_names{2} = 'BL3-H';
    bead_ch_names{3} = 'GL3-H';

    sidescatter_ch = gl1a_chn;

    ch(1) = sidescatter_ch;
    ch(2) = strmatch('BL3-H', {fcshdr.par.name});  %chlorophyl channe;
    ch(3) = strmatch('GL3-H', {fcshdr.par.name});
  

    temp_table = table;
    temp_table.filename = string(beadlist(b).name);
    temp_table.time = datetime([fcshdr.date, ' ', fcshdr.starttime]);
    temp_table.parname = {fcshdr.par.name};
    temp_table.OD2_filter = 'GL1';
    temp_table.hv = cell2mat({fcshdr.par.hv});
    
    n_clust = 3;
             
        % cluster analysis
        chl_noise_cut = 200;
        c = -1*ones(size(fcsdat,1),1); % c is a vector of cluster indeces, -1 indicates noise
        tempi = find(fcsdat(:,ch(2))>chl_noise_cut & fcsdat(:,ch(2))<1e6); %consider particles above threshhold
        minpts = floor(length(tempi)/75); %this is a parameter for clustering algorithm
        ctemp = dbscan(real(log10(fcsdat(tempi,ch)+1)), .005, minpts, 'distance', 'squaredeuclidean');
        c(tempi) = ctemp;
        nc = max(c); %this is number of cluster indexes
        meanc = NaN(length(ch), nc); %get mean signal on each channel for each cluster
        for iii = 1:nc
            meanc(:,iii) = mean(fcsdat(c==iii,ch));
        end
        clear iii
        [~,j] = sort(meanc(3,:));
        if (nc > 1 & meanc(1:2,j(1))./meanc(1:2,j(2)) > .8 & meanc(1:2,j(1))./meanc(1:2,j(2)) <1.5 & meanc(3,j(1))./meanc(3,j(2)) > 0.3 & meanc(3,j(1))./meanc(3,j(2)) < 1.5)
            c(c==j(1)) = j(2);
            for iii = 1:nc
                meanc(:,iii) = mean(fcsdat(c==iii,ch));
            end
        end
        clear iii
        [~,c2] = min(meanc(3,:)); %smallest on GL3, 1 micron
        tt = find(meanc(1,:) < meanc(1,c2));
        [~,c1] = min(meanc(2,tt)); %smallest on SSC (GL1), 0.5 micron
        c1 = tt(c1);
        
        
        % second cluster the stuff offscale on chlorophyll
        tempi = find(fcsdat(:,ch(2))>=1e6);
        if length(tempi) ~= 0
            minpts = max([10 floor(length(tempi)/15)]); %dbscan parameters
            ctemp2 = dbscan(real(log10(fcsdat(tempi,ch([1,3]))+1)), .003, minpts, 'distance', 'squaredeuclidean');
            c(tempi(ctemp2>0)) = ctemp2(ctemp2>0)+max(ctemp);
            nc = max(ctemp2);
            meanc = NaN(length(ch), nc);
            for ii = 1:nc
                meanc(:,ii) = mean(fcsdat(c==ii+max(ctemp),ch));
            end
        end
        [~,c3] = min(meanc(1,:)); %smallest on SSC (GL1)
        c3 = c3+max(ctemp);
        
        
        class = zeros(length(c), 1);
        class(c==c1) = 1; %0.5 micron
        class(c==c2) = 2; %1 micron
        class(c==c3) = 3; %6 micron
        
        
            %% get center values for each bead cluster on OD2 channel and not 
                m05_noOd2 = centend(fcsdat, class, 1, ssca_chn);
                m1_noOd2 = centend(fcsdat, class, 2, ssca_chn);
                m6_noOd2 = centend(fcsdat, class, 3, ssca_chn);
            
                m05_od2 = centend(fcsdat, class, 1, gl1a_chn);
                m1_od2 = centend(fcsdat, class, 2, gl1a_chn);
                m6_od2 = centend(fcsdat, class, 3, gl1a_chn);
                
                hv_od2 = fcshdr.par(gl1a_chn).hv;
                hv_noOd2 = fcshdr.par(ssca_chn).hv;
                
                od2_chn = gl1a_chn;
                noOd2_chn = ssca_chn; 
                

 %%  plot bead clusters for visual verification
            figure(99)
            clf
            set(gcf,'Position', [100 125 1100 650])
            
            
            figure(99), clf
            colormap(lines(4));
            subplot(2,2,1)
            scatter(fcsdat(:,ch(1)), fcsdat(:,ch(2)), 1, class, 'filled')
            set(gca, 'XScale', 'log', 'YScale', 'log')
            xlim([10 1.2e6]); ylim([10 1.2e6]);
            xlabel(bead_ch_names{1}); ylabel(bead_ch_names{2})
            
            subplot(2,2,2)
            scatter(fcsdat(:,ch(1)), fcsdat(:,ch(3)), 1, class, 'filled')
            set(gca, 'XScale', 'log', 'YScale', 'log')
            xlim([10 1.2e6]); ylim([10 1.2e6]);
            xlabel(bead_ch_names{1}); ylabel(bead_ch_names{3})
            
 
            subplot(2,2,3); hold on;
            cmap = colormap(lines(4));
            bins = logspace(1, 6.1, 256);
            histogram(fcsdat(class==1, od2_chn), bins, 'FaceColor', cmap(2,:), 'EdgeColor', 'none')
            histogram(fcsdat(class==2, od2_chn), bins, 'FaceColor', cmap(3,:), 'EdgeColor', 'none')
            histogram(fcsdat(class==3, od2_chn), bins, 'FaceColor', cmap(4,:), 'EdgeColor', 'none')
            yl = ylim;
            xlim([10 1.2e6])
            plot([m05_od2 m05_od2], [yl(1) yl(2)], '-k', 'LineWidth', 1)
            plot([m1_od2 m1_od2], [yl(1) yl(2)], '-k', 'LineWidth', 1)
            plot([m6_od2 m6_od2], [yl(1) yl(2)], '-k', 'LineWidth', 1)
            set(gca, 'XScale', 'log')
            legend('0.5 \mum', '1 \mum', '6 \mum')
            xlabel('OD2 SSC')
            

            subplot(2,2,4); hold on;
            cmap = colormap(lines(4));
            bins = logspace(1, 6.1, 256);
            histogram(fcsdat(class==1, noOd2_chn), bins, 'FaceColor', cmap(2,:), 'EdgeColor', 'none')
            histogram(fcsdat(class==2, noOd2_chn), bins, 'FaceColor', cmap(3,:), 'EdgeColor', 'none')
            histogram(fcsdat(class==3, noOd2_chn), bins, 'FaceColor', cmap(4,:), 'EdgeColor', 'none')
            yl = ylim;
            xlim([10 1.2e6])
            plot([m05_noOd2 m05_noOd2], [yl(1) yl(2)], '-k', 'LineWidth', 1)
            plot([m1_noOd2 m1_noOd2], [yl(1) yl(2)], '-k', 'LineWidth', 1)
            plot([m6_noOd2 m6_noOd2], [yl(1) yl(2)], '-k', 'LineWidth', 1)
            set(gca, 'XScale', 'log')
            legend('0.5 \mum', '1 \mum', '6 \mum')
            xlabel('SSC No Filter')


            subplot(2,2,1), title(beadlist(b).name, 'interpreter', 'none')
            subplot(2,2,3), title([datestr(temp_table.time) '; ' bead_ch_names{1} ' hv = ' num2str(fcshdr.par(ch(1)).hv) ' (' fcshdr.par(ch(1)).name ')'])
            

            print(figure(99), fullfile(beadfigpath, regexprep(beadlist(b).name, '.fcs', '.png')), '-dpng')
            

        if m1_od2<(m6_od2) && m05_od2<(m1_od2)
            QC_flag = 1;
        else
            QC_Flag = 0; 


        end
        temp_table.QC_flag = logical(QC_flag);




    for iii = 1:n_clust
        nstr = num2str(iii);
        temp_table.(['mean' nstr]) = mean(fcsdat(class==iii,:));
        temp_table.(['std' nstr]) = std(fcsdat(class==iii,:));
        temp_table.(['median' nstr]) = median(fcsdat(class==iii,:));
        temp_table.(['number' nstr]) = sum(class==iii);
        temp_table.(['centend' nstr])= centend(fcsdat, class, iii, ch(1));
    end
    clear iii
    
        %Store center on SSC-H channel only
        temp_table.OD2centers = [m05_od2 m1_od2 m6_od2];
        temp_table.OD2_hv = hv_od2;
    
        temp_table.NoOD2centers = [m05_noOd2 m1_noOd2 m6_noOd2];
        temp_table.NoOD2_hv = hv_noOd2;


    %store beadstat statistics
    beadstat(b,:) = temp_table;

end

notes = ['beads processed ' string(datetime)];

% save beadstat
save([beadpath filesep 'outputs\beadstat.mat'],'beadstat', 'notes', 'bead_ch_names')

disp('Result file saved:')
disp([beadpath filesep 'outputs\beadstat.mat'])




function c = centend(fcsdat, class_vec, class, channel)

bins = logspace(1, 6.1, 1000);
[n, edg] = histcounts(fcsdat(class_vec==class, channel), bins);
edg = edg(1:end-1) + (edg(2) - edg(1))/2;
s = smooth(n', 5);
[~, ind] = max(s);
c = edg(ind);

end


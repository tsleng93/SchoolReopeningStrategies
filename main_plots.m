%This .m file produces Figure 1 and Figure 2 of the paper
%N.B. Violinplot-Matlab is required to generate Figure 1 b), c), and d)
%N.B breakyaxis is required to generate Figure 1 e)

    %Positive test profiles
    PCR_test_sym = readtable('PCR_Curve_summary.csv');
    PCR_test_sym = table2array(PCR_test_sym(:, 2:4));
    PCR_test_asym = csvread('PCR_Curve_asym.csv');

    lat_test_sym = readtable('lat_Curve_summary.csv');
    lat_test_sym = table2array(lat_test_sym(:, 2:4));
    lat_test_asym = csvread('lat_Curve_asym.csv');
       
    runs = 10000;    
    isolationparams = zeros(runs,11); 
    
    %generate parameters
    rng(1);
    for i = 1:runs
        isolationparams(i,1) = 3*rand+1;
        isolationparams(i,2) = 0.0013133;
        isolationparams(i,3) = 1;
        isolationparams(i,4) = 1;
        isolationparams(i,5) = 1;
        isolationparams(i,6) = 0.02;
        isolationparams(i,7) = 0.2;
        isolationparams(i,8) = 0;
        isolationparams(i,9) = 0.4*rand + 0.3;
        isolationparams(i,10) = 0.12 + 0.19*rand;
        isolationparams(i,11) = 0;
    end
    
    testingbaselineparams = isolationparams;
    testingbaselineparams(:,3) = 0;
    
    testinglowparams = testingbaselineparams;
    testinglowparams(:,4) = 2;
    testinghighparams = testingbaselineparams;
    testinghighparams(:,4) = 3;
    
    notestingparams = testingbaselineparams;
    notestingparams(:,4) = 4;
    
    weeklyparams = testingbaselineparams;
    weeklyparams(:,4) = 5;
    
    
    weeklyscparams = testingbaselineparams;
    weeklyscparams(:,4) = 7;

    %run simulatios
    
    
    %run simulations
    tic
    for i = 1:runs
        
        
       
        %isolation
        history = Interactingyeargroups(isolationparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
        [Prev_true, Prev_asym, Known_asym, pos_tests, Absences, Tot_Infected, Schooldays_missed, ~, ~,AsymsCaptured,AsymsTotal, ~, ~, Days_tested] = Modeloutputs(history);
        Prevs_1(i,:) = sum(Prev_true);
        Asym_tot_1(i,:) = sum(Prev_asym);
        Asym_known_1(i,:) = sum(Known_asym);
        tests_1(i,:) = sum(pos_tests);
        Absences_1(i,:) = sum(Absences);
        Tot_Infected_1(i) = sum(Tot_Infected);
        Schooldays_missed_1(i, :) = Schooldays_missed(:);
        Days_tested_1(i,:) = Days_tested;
        AsymsCaptured_1(i) = sum(AsymsCaptured);
        AsymsTotal_1(i) = sum(AsymsTotal);
        [Incidence, ~, Num_PCR, ~, School_Inc] = Moremodeloutputs(history);
        Incidence_1(i,:) = sum(Incidence);
        Num_PCR_1(i) = Num_PCR;
        School_Inc_1(i) = School_Inc;
        
        %baseline testing
        history = Interactingyeargroups(testingbaselineparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
        [Prev_true, Prev_asym, Known_asym, pos_tests, Absences, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, ~, ~, Days_tested, PresymsCaptured, day_caught_vec] = Modeloutputs(history); 
        Prevs_2(i,:) = sum(Prev_true);
        Asym_tot_2(i,:) = sum(Prev_asym);
        Asym_known_2(i,:) = sum(Known_asym);
        tests_2(i,:) = sum(pos_tests);
        Absences_2(i,:) = sum(Absences);
        Tot_Infected_2(i) = sum(Tot_Infected);
        Schooldays_missed_2(i,:) = Schooldays_missed(:);
        Days_tested_2(i,:) = Days_tested;
        AsymsCaptured_2(i) = sum(AsymsCaptured);
        AsymsTotal_2(i) = sum(AsymsTotal);
        PresymsCaptured_2(i) = sum(PresymsCaptured);
        SymsTotal_2(i) = Tot_Infected_2(i) - AsymsTotal_2(i);
        Days_caught_2(i, :) = day_caught_vec;
        [Incidence, Num_lat, Num_PCR, Num_conf, School_Inc] = Moremodeloutputs(history);
        Incidence_2(i,:) = sum(Incidence);
        Num_lat_2(i) = Num_lat;
        Num_PCR_2(i) = Num_PCR;
        Num_conf_2(i) = Num_conf;
        School_Inc_2(i) = School_Inc;      
        
        
        %only isolating positives
        history = Interactingyeargroups(notestingparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
        [Prev_true, Prev_asym, Known_asym, pos_tests, Absences, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, ~, ~, Days_tested] = Modeloutputs(history);
        Prevs_5(i,:) = sum(Prev_true);
        Asym_tot_5(i,:) = sum(Prev_asym);
        Asym_known_5(i,:) = sum(Known_asym);
        tests_5(i,:) = sum(pos_tests);
        Absences_5(i,:) = sum(Absences);
        Tot_Infected_5(i) = sum(Tot_Infected);
        Schooldays_missed_5(i,:) = Schooldays_missed(:);
        Days_tested_5(i,:) = Days_tested;
        AsymsCaptured_5(i) = sum(AsymsCaptured);
        AsymsTotal_5(i) = sum(AsymsTotal);
       [Incidence, Num_lat, Num_PCR, Num_conf, School_Inc] = Moremodeloutputs(history);
        Incidence_5(i,:) = sum(Incidence);
        Num_lat_5(i) = Num_lat;
        Num_PCR_5(i) = Num_PCR;
        Num_conf_5(i) = Num_conf;
        School_Inc_5(i) = School_Inc;
        
        
        
        %weekly tests
        history = Interactingyeargroups(weeklyparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
        [Prev_true, Prev_asym, Known_asym, pos_tests, Absences, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, ~, ~, Days_tested, PresymsCaptured] = Modeloutputs(history);
        Prevs_6(i,:) = sum(Prev_true);
        Asym_tot_6(i,:) = sum(Prev_asym);
        Asym_known_6(i,:) = sum(Known_asym);
        tests_6(i,:) = sum(pos_tests);
        Absences_6(i,:) = sum(Absences);
        Tot_Infected_6(i) = sum(Tot_Infected);
        Schooldays_missed_6(i,:) = Schooldays_missed(:);
        Days_tested_6(i,:) = Days_tested;
        AsymsCaptured_6(i) = sum(AsymsCaptured);
        AsymsTotal_6(i) = sum(AsymsTotal);
        PresymsCaptured_6(i) = sum(PresymsCaptured);
        SymsTotal_6(i) = Tot_Infected_6(i) - AsymsTotal_6(i);
        [Incidence, Num_lat, Num_PCR, Num_conf, School_Inc] = Moremodeloutputs(history);
        Incidence_6(i,:) = sum(Incidence);
        Num_lat_6(i) = Num_lat;
        Num_PCR_6(i) = Num_PCR;
        Num_conf_6(i) = Num_conf;
        School_Inc_6(i) = School_Inc;
        
        %weekly tests serial contact tracing
        history = Interactingyeargroups(weeklyscparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
        [Prev_true, Prev_asym, Known_asym, pos_tests, Absences, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, ~, ~, Days_tested, PresymsCaptured] = Modeloutputs(history);
        Prevs_8(i,:) = sum(Prev_true);
        Asym_tot_8(i,:) = sum(Prev_asym);
        Asym_known_8(i,:) = sum(Known_asym);
        tests_8(i,:) = sum(pos_tests);
        Absences_8(i,:) = sum(Absences);
        Tot_Infected_8(i) = sum(Tot_Infected);
        Schooldays_missed_8(i,:) = Schooldays_missed(:);
        Days_tested_8(i,:) = Days_tested;
        AsymsCaptured_8(i) = sum(AsymsCaptured);
        AsymsTotal_8(i) = sum(AsymsTotal);
        PresymsCaptured_8(i) = sum(PresymsCaptured);
        SymsTotal_8(i) = Tot_Infected_8(i) - AsymsTotal_8(i);
       [Incidence, Num_lat, Num_PCR, Num_conf, School_Inc] = Moremodeloutputs(history);
        Incidence_8(i,:) = sum(Incidence);
        Num_lat_8(i) = Num_lat;
        Num_PCR_8(i) = Num_PCR;
        Num_conf_8(i) = Num_conf;
        School_Inc_8(i) = School_Inc;
        
        
    end
    
    toc
    %}
        
    %Key stats 
    
    %simulations where mass testing has a higher number of cases
    sum(Tot_Infected_2 > Tot_Infected_1)
    % increase in cases compared to an isolation strategy
    (sum(Tot_Infected_2) - sum(Tot_Infected_1))/sum(Tot_Infected_1)
    %school days missed isolation
    mean(mean(Schooldays_missed_1))
    %school days missed mass testing
    mean(mean(Schooldays_missed_2))
    
    %asymptomatics missed
    100*(1 - sum(AsymsCaptured_2)/sum(AsymsTotal_2))
    
    
    
    %Updated figures
    
    %Figure 1A    
    p_1 = mean(Schooldays_missed_1');
    q_1 = 100*(Tot_Infected_1/1000);

    p_2 = mean(Schooldays_missed_2');
    q_2 = 100*(Tot_Infected_2/1000);
    
    p_3 = mean(Schooldays_missed_6');
    q_3 = 100*(Tot_Infected_6/1000);

    p_4 = mean(Schooldays_missed_8');
    q_4 = 100*(Tot_Infected_8/1000);

    figure;
    set(gcf, 'Position', [300, 300, 600, 400]);
    
    %Markers for legend
    plot(100, 100, '.', 'MarkerSize', 20, 'color', [1.00, 0.41, 0.16]); hold on
    plot(100, 100, '.', 'MarkerSize', 20, 'color', [0.30, 0.75, 0.93]);
    plot(100, 100, '.', 'MarkerSize', 20, 'color', [0.72, 0.27, 1.00]);
    plot(100, 100, '.', 'MarkerSize', 20, 'color', [0.39, 0.83, 0.07]);
    
    plot(q_1, p_1, '.', 'MarkerSize' ,5, 'color', [1, 0.41, 0.16], 'LineStyle', 'none');   
    plot(q_2,p_2, '.', 'MarkerSize', 5, 'color', [0.30, 0.75, 0.93], 'LineStyle', 'none');
    plot(q_3, p_3, '.', 'MarkerSize' ,5, 'color', [0.72, 0.27, 1.00], 'LineStyle', 'none');   
    plot(q_4,p_4, '.', 'MarkerSize', 5, 'color', [0.39, 0.83, 0.07], 'LineStyle', 'none');  

    xlabel('Total infected by end of half term (%)');
    ylabel('Mean school days missed per pupil');
    xlim([4 36]);
    ylim([0 18]);
    legend('Isolating year groups', 'serial contact testing', 'weekly mass testing', 'combined');
    set(gca, 'fontsize', 12);
    
    
    %Figure 1B
    figure;
    set(gcf, 'Position', [300, 300, 600, 400]);
    plot(100, 100, '.', 'MarkerSize', 20, 'color', [1.00, 0.41, 0.16]); hold on
    plot(100, 100, '.', 'MarkerSize', 20, 'color', [0.30, 0.75, 0.93]);
    plot(100, 100, '.', 'MarkerSize', 20, 'color', [0.72, 0.27, 1.00]);
    plot(100, 100, '.', 'MarkerSize', 20, 'color', [0.39, 0.83, 0.07]);
    plot(100, 100, '.', 'MarkerSize', 20, 'color', [0.65, 0.65, 0.65]);
    
    v = violinplot([100*(Tot_Infected_1/1000); 100*(Tot_Infected_2/1000); 100*(Tot_Infected_6/1000); 100*(Tot_Infected_8/1000); 100*(Tot_Infected_5/1000)]'); hold on
    v(1).ViolinColor = [1.00, 0.41, 0.16];
    v(2).ViolinColor = [0.30, 0.75, 0.93];
    v(3).ViolinColor = [0.72, 0.27, 1.00];
    v(4).ViolinColor = [0.39, 0.83, 0.07];
    v(5).ViolinColor = [0.65, 0.65, 0.65];
    xlim([0 6]);
    legend('Isolating year groups', 'serial contact testing', 'weekly mass testing', 'combined', 'no control'); %extra legend entries must be removed manually
    xticklabels({'(i)', '(ii)', '(iii)', '(iv)', '(v)'});
     set(gca, 'fontsize', 12);
     xlabel('Reopening stratagy');
     ylabel('Total infected by end of half term (%)');
     
     
    %Figure 1C
    figure; %violin curves must be removed manually
    set(gcf, 'Position', [300, 300, 600, 400]);   
    v = violinplot([100*(School_Inc_1/1000); 100*(School_Inc_2/1000); 100*(School_Inc_6/1000); 100*(School_Inc_8/1000); 100*(School_Inc_5/1000)]'); hold on
    v(1).ViolinColor = [1.00, 0.41, 0.16];
    v(2).ViolinColor = [0.30, 0.75, 0.93];
    v(3).ViolinColor = [0.72, 0.27, 1.00];
    v(4).ViolinColor = [0.39, 0.83, 0.07];
    v(5).ViolinColor = [0.65, 0.65, 0.65];
    xlim([0 6]);
    xticklabels({'(i)', '(ii)', '(iii)', '(iv)', '(v)'});
     set(gca, 'fontsize', 12);
     xlabel('Reopening stratagy');
     ylabel('In school incidence at end of term (%)');
     %}
     
    %Figure 1D 
    figure; 
    set(gcf, 'Position', [300, 300, 600, 400]);   
    v = violinplot([100*(AsymsCaptured_2./AsymsTotal_2); 100*(AsymsCaptured_6./AsymsTotal_6); 100*(AsymsCaptured_8)./AsymsTotal_8]');   
    v(1).ViolinColor = [0.30, 0.75, 0.93];
    v(2).ViolinColor = [0.72, 0.27, 1.00];
    v(3).ViolinColor = [0.39, 0.83, 0.07];
    xlim([-1 5]);
    xticklabels({ '(ii)', '(iii)', '(iv)'});
     set(gca, 'fontsize', 12);
     xlabel('Reopening stratagy');
     ylabel('% of asymptomatic cased identified');
     
     %Figure 1E %position of y-axis label, xticklabels, and yticklabels to
     %be adjusted manually, because of breakyaxis
     figure;
     set(gcf, 'Position', [300, 300, 600, 400]);   
     bar(1, mean(Schooldays_missed_1(:)), 'FaceColor', [1.00, 0.41, 0.16]); hold on
     bar(2, mean(Schooldays_missed_2(:)), 'FaceColor', [0.30, 0.75, 0.93]); hold on
     bar(3, mean(Schooldays_missed_6(:)), 'FaceColor', [0.72, 0.27, 1.00]); hold on
     bar(4, mean(Schooldays_missed_8(:)), 'FaceColor', [0.39, 0.83, 0.07]); hold on
     bar(5, mean(Schooldays_missed_5(:)), 'FaceColor', [0.65, 0.65, 0.65]); hold on
     xlim([0 6]);
     ylim([0 8]);
     xticklabels({'(i)', '(ii)', '(iii)', '(iv)', '(v)'});     
     xlabel('Reopening stratagy');
     ylabel('Mean school days missed per pupil');
     breakyaxis([0.6 7.6]);
     set(gca, 'fontsize', 12);
     
     
     %Figure 1F %xticklabels to be adjusted manually
     figure;
     set(gcf, 'Position', [300, 300, 600, 400]); 
     bar(2, mean(Num_lat_2/1000), 'FaceColor', [0.30, 0.75, 0.93]); hold on
     bar(3, mean(Num_lat_6/1000), 'FaceColor', [0.72, 0.27, 1.00]); hold on
     bar(4, mean(Num_lat_8/1000), 'FaceColor', [0.39, 0.83, 0.07]); hold on
     xticklabels({'(ii)', '(iii)', '(iv)'});     
     xlim([0 6]);
     xlabel('Reopening stratagy');
     ylabel('Mean number of LFTs per pupil');
     set(gca, 'fontsize', 12);
     
     
     %Figure 2A
     figure;
     set(gcf, 'Position', [300, 300, 600, 400]); 
     plot((1:56)/7 - 1, 100*mean(Prevs_5)/1000, 'linewidth', 1.5, 'color', [0.65, 0.65, 0.65]); hold on
     above =  quantile(100*Prevs_5/1000, 0.75);
     below =  quantile(100*Prevs_5/1000, 0.25);
     p = patch([(1:56)/7 - 1 fliplr((1:56)/7 - 1)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [0.65, 0.65, 0.65];
     plot((1:56)/7 - 1, 100*mean(Prevs_6)/1000, 'linewidth', 1.5, 'color', [0.72, 0.27, 1.00]); hold on
     above =  quantile(100*Prevs_6/1000, 0.75);
     below =  quantile(100*Prevs_6/1000, 0.25);
     p = patch([(1:56)/7 - 1 fliplr((1:56)/7 - 1)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [0.72, 0.27, 1.00];
     plot((1:56)/7 - 1, 100*mean(Prevs_2)/1000, 'linewidth', 1.5, 'color', [0.30, 0.75, 0.93]); hold on
     above =  quantile(100*Prevs_2/1000, 0.75);
     below =  quantile(100*Prevs_2/1000, 0.25);
     p = patch([(1:56)/7 - 1 fliplr((1:56)/7 - 1)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [0.30, 0.75, 0.93];
     plot((1:56)/7 - 1, 100*mean(Prevs_1)/1000, 'linewidth', 1.5, 'color', [1.00, 0.41, 0.16]); hold on
     above =  quantile(100*Prevs_1/1000, 0.75);
     below =  quantile(100*Prevs_1/1000, 0.25);
     p = patch([(1:56)/7 - 1 fliplr((1:56)/7 - 1)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [1.00, 0.41, 0.16];
     plot((1:56)/7 - 1, 100*mean(Prevs_8)/1000, 'linewidth', 1.5, 'color', [0.39, 0.83, 0.07]); hold on
     above =  quantile(100*Prevs_8/1000, 0.75);
     below =  quantile(100*Prevs_8/1000, 0.25);
     p = patch([(1:56)/7 - 1 fliplr((1:56)/7 - 1)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [0.39, 0.83, 0.07];
     xlim([0 7]);
     legend('no control', 'weekly mass testing', 'serial contact testing', 'Isolating year groups', 'combined');
     xlabel('Week');
     ylabel('Prevalence (%)');
     set(gca, 'fontsize', 12);
     
     %Figure 2B
     figure;
     
     plot((1:56)/7 - 1, 100*(sum(Asym_known_6)./sum(Asym_tot_6)), 'linewidth', 1.5, 'color', [0.72, 0.27, 1.00]); hold on
     above =  quantile(100*(Asym_known_6./Asym_tot_6), 0.75);
     below =  quantile(100*(Asym_known_6./Asym_tot_6), 0.25);
     p = patch([(1:56)/7 - 1 fliplr((1:56)/7 - 1)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [0.72, 0.27, 1.00];
     plot((1:56)/7 - 1, 100*(sum(Asym_known_2)./sum(Asym_tot_2)), 'linewidth', 1.5, 'color', [0.30, 0.75, 0.93]); hold on
     above =  quantile(100*(Asym_known_2./Asym_tot_2), 0.75);
     below =  quantile(100*(Asym_known_2./Asym_tot_2), 0.25);
     p = patch([(1:56)/7 - 1 fliplr((1:56)/7 - 1)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [0.30, 0.75, 0.93];
     plot((1:56)/7 - 1, 100*(sum(Asym_known_8)./sum(Asym_tot_8)), 'linewidth', 1.5, 'color', [0.39, 0.83, 0.07]); hold on
     above =  quantile(100*(Asym_known_8./Asym_tot_8), 0.75);
     below =  quantile(100*(Asym_known_8./Asym_tot_8), 0.25);
     p = patch([(1:56)/7 - 1 fliplr((1:56)/7 - 1)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [0.39, 0.83, 0.07];
     xlim([0 7]);
     xlabel('Week');
     ylabel('% of asymptomatic cases identified');
     set(gca, 'fontsize', 12);
     
     
     %Figure 2C
     figure;
     set(gcf, 'Position', [300, 300, 600, 400]); 
     plot((1:56)/7 - 1, 100*mean(Absences_5)/1000, 'linewidth', 1.5, 'color', [0.65, 0.65, 0.65]); hold on
     above =  quantile(100*Absences_5/1000, 0.75);
     below =  quantile(100*Absences_5/1000, 0.25);
     p = patch([(1:56)/7 - 1 fliplr((1:56)/7 - 1)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [0.65, 0.65, 0.65];
     
     plot((1:56)/7 - 1, 100*mean(Absences_6)/1000, 'linewidth', 1.5, 'color', [0.72, 0.27, 1.00]); hold on
     above =  quantile(100*Absences_6/1000, 0.75);
     below =  quantile(100*Absences_6/1000, 0.25);
     p = patch([(1:56)/7 - 1 fliplr((1:56)/7 - 1)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [0.72, 0.27, 1.00];
     plot((1:56)/7 - 1, 100*mean(Absences_2)/1000, 'linewidth', 1.5, 'color', [0.30, 0.75, 0.93]); hold on
     above =  quantile(100*Absences_2/1000, 0.75);
     below =  quantile(100*Absences_2/1000, 0.25);
     p = patch([(1:56)/7 - 1 fliplr((1:56)/7 - 1)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [0.30, 0.75, 0.93];
     plot((1:56)/7 - 1, 100*mean(Absences_1)/1000, 'linewidth', 1.5, 'color', [1.00, 0.41, 0.16]); hold on
     above =  quantile(100*Absences_1/1000, 0.75);
     below =  quantile(100*Absences_1/1000, 0.25);
     p = patch([(1:56)/7 - 1 fliplr((1:56)/7 - 1)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [1.00, 0.41, 0.16];
     plot((1:56)/7 - 1, 100*mean(Absences_8)/1000, 'linewidth', 1.5, 'color', [0.39, 0.83, 0.07]); hold on
     above =  quantile(100*Absences_8/1000, 0.75);
     below =  quantile(100*Absences_8/1000, 0.25);
     p = patch([(1:56)/7 - 1 fliplr((1:56)/7 - 1)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [0.39, 0.83, 0.07];
     xlim([0 7]);
     xlabel('Week');
     ylabel('% of pupils absent');
     set(gca, 'fontsize', 12);
     %}
     
     %Figure 2D
     figure;
     
     plot((1:56)/7 - 1, 100*mean(Days_tested_6)/5, 'linewidth', 1.5, 'color', [0.72, 0.27, 1.00]); hold on
     above =  quantile(100*Days_tested_6/5, 0.75);
     below =  quantile(100*Days_tested_6/5, 0.25);
     p = patch([(1:56)/7 - 1 fliplr((1:56)/7 - 1)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [0.72, 0.27, 1.00];
     plot((1:56)/7 - 1, 100*mean(Days_tested_2)/5, 'linewidth', 1.5, 'color', [0.30, 0.75, 0.93]); hold on
     above =  quantile(100*Days_tested_2/5, 0.75);
     below =  quantile(100*Days_tested_2/5, 0.25);
     p = patch([(1:56)/7 - 1 fliplr((1:56)/7 - 1)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [0.30, 0.75, 0.93];
     plot((1:56)/7 - 1, 100*mean(Days_tested_8)/5, 'linewidth', 1.5, 'color', [0.39, 0.83, 0.07]); hold on
     above =  quantile(100*Days_tested_8/5, 0.75);
     below =  quantile(100*Days_tested_8/5, 0.25);
     p = patch([(1:56)/7 - 1 fliplr((1:56)/7 - 1)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [0.39, 0.83, 0.07];
     xlim([0 7]);
     xlabel('Week');
     ylabel('% of pupils tested in school');
     set(gca, 'fontsize', 12);
     

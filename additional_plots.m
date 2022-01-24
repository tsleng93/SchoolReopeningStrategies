%This file generate Figures 3, S7-S10, %%To generate Figure S11 - change
%the parameter determining frequency dependent transmission for Figure S10


clear

%Positive test profiles
    PCR_test_sym = readtable('PCR_Curve_summary.csv');
    PCR_test_sym = table2array(PCR_test_sym(:, 2:4));
    PCR_test_asym = csvread('PCR_Curve_asym.csv');

    lat_test_sym = readtable('lat_Curve_summary.csv');
    lat_test_sym = table2array(lat_test_sym(:, 2:4));
    lat_test_asym = csvread('lat_Curve_asym.csv');
       

runs = 10; %For quick run
%runs = 2000; %Number of simulations for paper
    %set baseline parameters
    rng(1);
    for i = 1:runs

        isolationparams(i,1) = 1 +4*rand;       
        isolationparams(i,2) = 0.001026404983400;
        isolationparams(i,3) = 1;
        isolationparams(i,4) = 1;
        isolationparams(i,5) = 1;
        isolationparams(i,6) = 0.02;
        isolationparams(i,7) = 0.2;
        isolationparams(i,8) = 0;
        isolationparams(i,9) = 0.4*rand + 0.3;
        isolationparams(i,10) = 0.12 + 0.19*rand;
        isolationparams(i,11) = 0;
        isolationparams(i,12) = 200;
        isolationparams(i,13) = 5;
        isolationparams(i,14) = 0;
        isolationparams(i,15) = 0;

    end
  
%


%% Figure 3 - Impact of Uptake/participation %%

   tic
for j = 1:101
        
j
        
    isolationparams(:,3) = 0.01*(j-1);
        
    testingbaselineparams = isolationparams;
    %testingbaselineparams(:,3) = 0;
    
    testinglowparams = testingbaselineparams;
    testinglowparams(:,4) = 2;
    testinghighparams = testingbaselineparams;
    testinghighparams(:,4) = 3;

    
    biweeklyparams = testingbaselineparams;
    biweeklyparams(:,4) = 6;
    
    
    
    biweeklyscparams = testingbaselineparams;
    biweeklyscparams(:,4) = 8;
    
        
    notestingparams = testingbaselineparams;
    notestingparams(:,4) = 4;

    testisolparams = testingbaselineparams;
    testisolparams(:,4) = 14;
        
        parfor i = 1:runs
            
            
            %baseline testing
            
            
            
            history = Interactingyeargroups(testingbaselineparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
            [Prev_true, Prev_asym, Known_asym, pos_tests, Absences, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, ~, ~, Days_tested, PresymsCaptured, day_caught_vec] = Modeloutputs(history); 
            Tot_Infected_2(i,j) = sum(Tot_Infected);
            Schooldays_missed_2(i,j) = sum(Schooldays_missed(:));
            
            %weekly mass testing 
            history = Interactingyeargroups(biweeklyparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
            [Prev_true, Prev_asym, Known_asym, pos_tests, Absences, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, ~, ~, Days_tested, PresymsCaptured, day_caught_vec] = Modeloutputs(history); 
            Tot_Infected_6(i,j) = sum(Tot_Infected);
            Schooldays_missed_6(i,j) = sum(Schooldays_missed(:));
            
            
            %weekly mass testing + SCT
             history = Interactingyeargroups(biweeklyscparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
            [Prev_true, Prev_asym, Known_asym, pos_tests, Absences, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, ~, ~, Days_tested, PresymsCaptured, day_caught_vec] = Modeloutputs(history); 
            Tot_Infected_8(i,j) = sum(Tot_Infected);
            Schooldays_missed_8(i,j) = sum(Schooldays_missed(:));
            
            %biweekly mass testing + SCT
            history = Interactingyeargroups(testisolparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
            [Prev_true, Prev_asym, Known_asym, pos_tests, Absences, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, ~, ~, Days_tested, PresymsCaptured, day_caught_vec] = Modeloutputs(history); 
            Tot_Infected_9(i,j) = sum(Tot_Infected);
            Schooldays_missed_9(i,j) = sum(Schooldays_missed(:)); 
          
            %{
            history = Interactingyeargroups(notestingparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
            [Prev_true, Prev_asym, Known_asym, pos_tests, Absences, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, ~, ~, Days_tested, PresymsCaptured, day_caught_vec] = Modeloutputs(history); 
            Tot_Infected_5(i,j) = sum(Tot_Infected);
            Schooldays_missed_5(i,j) = sum(Schooldays_missed(:)); 
            %}
            
        end
        
        
end
toc
    

figure;
    for i = 1:101
    
    Tot_Infected_5(:,i) = Tot_Infected_6(:,end);
    Tot_Infected_1(:,i) = Tot_Infected_2(:,end);
    end

    
    xvals = 0:1:100;
    
    Tot_Infected_5  = Tot_Infected_5(:,end:-1:1);
    Tot_Infected_6  = Tot_Infected_6(:,end:-1:1);
    Tot_Infected_8  = Tot_Infected_8(:,end:-1:1);
    Tot_Infected_9  = Tot_Infected_9(:,end:-1:1);
    Tot_Infected_2  = Tot_Infected_2(:,end:-1:1);
    Tot_Infected_1  = Tot_Infected_1(:,end:-1:1);
    
    
     above =  quantile(100*Tot_Infected_1/1000, 0.75);
     below =  quantile(100*Tot_Infected_1/1000, 0.25);
     p = patch([xvals fliplr(xvals)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [1.00, 0.41, 0.16];
    
     above =  quantile(100*Tot_Infected_9/1000, 0.75);
     below =  quantile(100*Tot_Infected_9/1000, 0.25);
     p = patch([xvals fliplr(xvals)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [0.93, 0.69, 0.13];


        
     above =  quantile(100*Tot_Infected_2/1000, 0.75);
     below =  quantile(100*Tot_Infected_2/1000, 0.25);
     p = patch([xvals fliplr(xvals)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [0.30, 0.75, 0.93];
     
     



     above =  quantile(100*Tot_Infected_8/1000, 0.75);
     below =  quantile(100*Tot_Infected_8/1000, 0.25);
     p = patch([xvals fliplr(xvals)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [0.39, 0.83, 0.07];
    
    
     above =  quantile(100*Tot_Infected_6/1000, 0.75);
     below =  quantile(100*Tot_Infected_6/1000, 0.25);
     p = patch( [xvals fliplr(xvals)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [0.72, 0.27, 1.00];
    

     

    
    
     above =  quantile(100*Tot_Infected_5/1000, 0.75);
     below =  quantile(100*Tot_Infected_5/1000, 0.25);
     p = patch([xvals fliplr(xvals)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [0.65, 0.65, 0.65];
     
             
                       plot(xvals, 100*mean(Tot_Infected_1)/1000, 'linewidth', 1.5, 'color', [1.00, 0.41, 0.16]); hold on
                         plot(xvals, 100*mean(Tot_Infected_9)/1000, 'linewidth', 1.5, 'color', [0.93, 0.69, 0.13]); hold on

             plot(xvals, 100*mean(Tot_Infected_2)/1000, 'linewidth', 1.5, 'color', [0.30, 0.75, 0.93]); hold on
      plot(xvals, 100*mean(Tot_Infected_8)/1000, 'linewidth', 1.5, 'color', [0.39, 0.83, 0.07]); hold on
            plot(xvals, 100*mean(Tot_Infected_6)/1000, 'linewidth', 1.5, 'color', [0.72, 0.27, 1.00]); hold on
             plot(xvals, 100*mean(Tot_Infected_5)/1000, 'linewidth', 1.5, 'color', [0.65, 0.65, 0.65]); hold on

    
    
     xlim([0 100]);
     legend('isolating year groups', 'twice weekly mass tests + isolating year groups', 'serial contact testing', 'twice weekly mass tests + serial contact testing', 'twice weekly mass tests', 'no control');
     xlabel('% pupils participating in lateral flow testing ');
     ylabel('Total infected by end of half term (%)');
     set(gca, 'fontsize', 14);
     ylim([5 25]);
     
 figure;
 
     figure;
   
    for i = 1:101
    Schooldays_missed_5(:,i) = Schooldays_missed_6(:,101);
    Schooldays_missed_1(:,i) = Schooldays_missed_2(:,101);
    end
    
    Schooldays_missed_6 = Schooldays_missed_6(:,end:-1:1);
    Schooldays_missed_8 = Schooldays_missed_8(:, end:-1:1);
    Schooldays_missed_9 = Schooldays_missed_9(:,end:-1:1);
    Schooldays_missed_2 = Schooldays_missed_2(:,end:-1:1);
    
    xvals = 0:1:100;
    
    
     plot(xvals, mean(Schooldays_missed_5)/1000, 'linewidth', 1.5, 'color', [0.65, 0.65, 0.65]); hold on
     above =  quantile(Schooldays_missed_5/1000, 0.75);
     below =  quantile(Schooldays_missed_5/1000, 0.25);
     p = patch([xvals fliplr(xvals)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [0.65, 0.65, 0.65];
     

     plot(xvals, mean(Schooldays_missed_6)/1000, 'linewidth', 1.5, 'color', [0.72, 0.27, 1.00]); hold on
     above =  quantile(Schooldays_missed_6/1000, 0.75);
     below =  quantile(Schooldays_missed_6/1000, 0.25);
     p = patch( [xvals fliplr(xvals)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [0.72, 0.27, 1.00];
          plot(xvals, mean(Schooldays_missed_2)/1000, 'linewidth', 1.5, 'color', [0.30, 0.75, 0.93]); hold on
     above =  quantile(Schooldays_missed_2/1000, 0.75);
     below =  quantile(Schooldays_missed_2/1000, 0.25);
     p = patch([xvals fliplr(xvals)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [0.30, 0.75, 0.93];
     
     
     plot(xvals, mean(Schooldays_missed_1)/1000, 'linewidth', 1.5, 'color', [1.00, 0.41, 0.16]); hold on
     above =  quantile(Schooldays_missed_1/1000, 0.75);
     below =  quantile(Schooldays_missed_1/1000, 0.25);
     p = patch([xvals fliplr(xvals)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [1.00, 0.41, 0.16];


     plot(xvals, mean(Schooldays_missed_8)/1000, 'linewidth', 1.5, 'color', [0.39, 0.83, 0.07]); hold on
     above =  quantile(Schooldays_missed_8/1000, 0.75);
     below =  quantile(Schooldays_missed_8/1000, 0.25);
     p = patch([xvals fliplr(xvals)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [0.39, 0.83, 0.07];
     
     plot(xvals, mean(Schooldays_missed_9)/1000, 'linewidth', 1.5, 'color', [0.93, 0.69, 0.13]); hold on
     above =  quantile(Schooldays_missed_9/1000, 0.75);
     below =  quantile(Schooldays_missed_9/1000, 0.25);
     p = patch([xvals fliplr(xvals)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [0.93, 0.69, 0.13];
     xlim([0 100]);
     legend('no control', 'twice weekly mass tests', 'serial contact testing', 'isolating year groups', 'twice weekly mass tests + serial contact testing', 'twice weekly mass tests + isolating year groups');
     xlabel('% pupils participating in lateral flow testing ');
     ylabel('Mean school days missed per pupil');
     set(gca, 'fontsize', 14);

%}
     
   

%% Supplementary Figure S7 %% 
%Impact of infection on school days 

clearvars -except runs PCR_test_sym PCR_test_asym lat_test_sym lat_test_asym 
  for i = 1:runs

        isolationparams(i,1) = 1 +4*rand;       
        isolationparams(i,2) = 0.001026404983400;
        isolationparams(i,3) = 1;
        isolationparams(i,4) = 1;
        isolationparams(i,5) = 1;
        isolationparams(i,6) = 0.02;
        isolationparams(i,7) = 0.2;
        isolationparams(i,8) = 0;
        isolationparams(i,9) = 0.4*rand + 0.3;
        isolationparams(i,10) = 0.12 + 0.19*rand;
        isolationparams(i,11) = 0;
        isolationparams(i,12) = 200;
        isolationparams(i,13) = 5;
        isolationparams(i,14) = 0;
        isolationparams(i,15) = 0;

    end


for j = 1:2
    j
    
    isolationparams(:,3) = 1; 
    
    if j == 2
        isolationparams(:,14) = 1;
    end
    
    
    testingbaselineparams = isolationparams;
    testingbaselineparams(:,3) = 0;
    
    testinglowparams = testingbaselineparams;
    testinglowparams(:,4) = 2;
    testinghighparams = testingbaselineparams;
    testinghighparams(:,4) = 3;
    
    notestingparams = testingbaselineparams;
    notestingparams(:,4) = 4;
    
    biweeklyparams = testingbaselineparams;
    biweeklyparams(:,4) = 6;   
    biweeklyscparams = testingbaselineparams;
    biweeklyscparams(:,4) = 8;
    
    
    testisolparams = testingbaselineparams;
    testisolparams(:,4) = 14;
    
parfor i = 1:runs
                   %baseline testing
            

            
            history = Interactingyeargroups(testingbaselineparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
            [Prev_true, Prev_asym, Known_asym, pos_tests, Absences, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, ~, ~, Days_tested, PresymsCaptured, day_caught_vec] = Modeloutputs(history); 
            Tot_Infected_2(i,j) = sum(Tot_Infected);
            Schooldays_missed_2(i,j) = sum(Schooldays_missed(:));
            
            %weekly mass testing 
            history = Interactingyeargroups(biweeklyparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
            [Prev_true, Prev_asym, Known_asym, pos_tests, Absences, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, ~, ~, Days_tested, PresymsCaptured, day_caught_vec] = Modeloutputs(history); 
            Tot_Infected_6(i,j) = sum(Tot_Infected);
            Schooldays_missed_6(i,j) = sum(Schooldays_missed(:));
            
            
            %weekly mass testing + SCT
             history = Interactingyeargroups(biweeklyscparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
            [Prev_true, Prev_asym, Known_asym, pos_tests, Absences, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, ~, ~, Days_tested, PresymsCaptured, day_caught_vec] = Modeloutputs(history); 
            Tot_Infected_8(i,j) = sum(Tot_Infected);
            Schooldays_missed_8(i,j) = sum(Schooldays_missed(:));
            
            %biweekly mass testing + SCT
            history = Interactingyeargroups(testisolparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
            [Prev_true, Prev_asym, Known_asym, pos_tests, Absences, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, ~, ~, Days_tested, PresymsCaptured, day_caught_vec] = Modeloutputs(history); 
            Tot_Infected_9(i,j) = sum(Tot_Infected);
            Schooldays_missed_9(i,j) = sum(Schooldays_missed(:)); 
          
            history = Interactingyeargroups(notestingparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
            [Prev_true, Prev_asym, Known_asym, pos_tests, Absences, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, ~, ~, Days_tested, PresymsCaptured, day_caught_vec] = Modeloutputs(history); 
            Tot_Infected_5(i,j) = sum(Tot_Infected);
            Schooldays_missed_5(i,j) = sum(Schooldays_missed(:)); 

        
end
 
end


figure;
A = [100*(Tot_Infected_9/1000), NaN*ones(runs,1), 100*(Tot_Infected_2/1000), NaN*ones(runs,1), 100*(Tot_Infected_8/1000), NaN*ones(runs,1), 100*(Tot_Infected_6/1000)]';
    plot(100, 100, '.', 'MarkerSize', 20, 'color', [1.00, 1.00, 0.00]); hold on
    plot(100, 100, '.', 'MarkerSize', 20, 'color', [0.64, 0.08, 0.18]);

v = violinplot(A'); hold on


xticks([1.5, 4.5, 7.5, 10.5]);
 xticklabels({'mass testing + isol', 'serial contact testing', 'mass testing + SCT', 'mass testing'});

for i = 1:size(A,1)
   if mod(i,3) == 1
           v(i).ViolinColor = [1.0, 1.0, 0.0];
   elseif mod(i,3) == 2
       v(i).ViolinColor = [0.64, 0.08, 0.18];
   end
end

set(gcf, 'Position', [300, 300, 600, 450]);
set(gca, 'fontsize', 14);
xlabel('Reopening strategy');
ylabel('Total infected by end of half term (%)');
legend('no risk of transmitting on +ve test day', 'risk of transmitting on +ve test day');
     
ylim([0 40]);

%}

clearvars -except runs PCR_test_sym PCR_test_asym lat_test_sym lat_test_asym 
%% Figure S8 - Impact of levels of immunity %%
  for i = 1:runs

        isolationparams(i,1) = 1 +4*rand;       
        isolationparams(i,2) = 0.001026404983400;
        isolationparams(i,3) = 1;
        isolationparams(i,4) = 1;
        isolationparams(i,5) = 1;
        isolationparams(i,6) = 0.02;
        isolationparams(i,7) = 0.2;
        isolationparams(i,8) = 0;
        isolationparams(i,9) = 0.4*rand + 0.3;
        isolationparams(i,10) = 0.12 + 0.19*rand;
        isolationparams(i,11) = 0;
        isolationparams(i,12) = 200;
        isolationparams(i,13) = 5;
        isolationparams(i,14) = 0;
        isolationparams(i,15) = 0;

    end

for j = 1:91
   
    %parameters for isolation
    
    isolationparams(:,7) = (j-1)*0.01;
    
    %parameters for serial contact testing
    testingbaselineparams = isolationparams;
    testingbaselineparams(:,3) = 0;
    
    %parameters for twice weekly mass testing
    biweeklyparams = testingbaselineparams;
    biweeklyparams(:,4) = 6;
    
    %parameters for combined
    biweeklyscparams = testingbaselineparams;
    biweeklyscparams(:,4) = 8;
            
    testisolparams = testingbaselineparams;
    testisolparams(:,4) = 14;
        
    notestingparams = testingbaselineparams;
    notestingparams(:,4) = 4;
    
   tic
   parfor i = 1:runs
        %isolation
        history = Interactingyeargroups(isolationparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
       [~, ~, ~, ~, ~, Tot_Infected,Schooldays_missed , Infected_within_school, Infected_during_term, AsymsCaptured, AsymsTotal, Peak_infect, Tot_Isolating] =   Modeloutputs(history);
       Tot_Infected_1(j,i) = sum(Tot_Infected);
       AsymsCaptured_1(j,i) = sum(AsymsCaptured);
       AsymsTotal_1(j,i) = sum(AsymsTotal);
       Peak_infect_1(j,i) = Peak_infect;
       Tot_Isolating_1(j,i) = mean(Tot_Isolating(:));
       Infected_during_term_1(j,i) = sum(Infected_during_term);
       Infected_within_school_1(j,i) = sum(Infected_within_school);
       Schooldays_missed_1(j,i) = mean(Schooldays_missed(:));
       
        %baseline testing
        history = Interactingyeargroups(testingbaselineparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
       [~, ~, ~, ~, ~, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, Peak_infect, Tot_Isolating] =   Modeloutputs(history);
       Tot_Infected_2(j,i) = sum(Tot_Infected);
       AsymsCaptured_2(j,i) = sum(AsymsCaptured);
       AsymsTotal_2(j,i) = sum(AsymsTotal);
       Peak_infect_2(j,i) = Peak_infect;
       Tot_Isolating_2(j,i) = mean(Tot_Isolating(:));
       Schooldays_missed_2(j,i) = mean(Schooldays_missed(:));
       [Incidence, Num_lat, Num_PCR, Num_conf, School_Inc, Num_lat_t] = Moremodeloutputs(history);
       Num_lat_2(j,i) = Num_lat;


       
       
       %weekly mass testing
        history = Interactingyeargroups(biweeklyparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
       [~, ~, ~, ~, ~, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, Peak_infect, Tot_Isolating] =   Modeloutputs(history);
       Tot_Infected_3(j,i) = sum(Tot_Infected);
       AsymsCaptured_3(j,i) = sum(AsymsCaptured);
       AsymsTotal_3(j,i) = sum(AsymsTotal);
       Peak_infect_3(j,i) = Peak_infect;
       Tot_Isolating_3(j,i) = mean(Tot_Isolating(:));
       Schooldays_missed_3(j,i) = mean(Schooldays_missed(:));
       [Incidence, Num_lat, Num_PCR, Num_conf, School_Inc, Num_lat_t] = Moremodeloutputs(history);
       Num_lat_3(j,i) = Num_lat;

        %weekly mass testing and serial contact testing
        history = Interactingyeargroups(biweeklyscparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
       [~, ~, ~, ~, ~, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, Peak_infect, Tot_Isolating] =   Modeloutputs(history);
       Tot_Infected_4(j,i) = sum(Tot_Infected);
       AsymsCaptured_4(j,i) = sum(AsymsCaptured);
       AsymsTotal_4(j,i) = sum(AsymsTotal);
       Peak_infect_4(j,i) = Peak_infect;
       Tot_Isolating_4(j,i) = mean(Tot_Isolating(:));
       Schooldays_missed_4(j,i) = mean(Schooldays_missed(:));
       [Incidence, Num_lat, Num_PCR, Num_conf, School_Inc, Num_lat_t] = Moremodeloutputs(history);
       Num_lat_4(j,i) = Num_lat;
    
       
        % testing and isol
       history = Interactingyeargroups(testisolparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
       [~, ~, ~, ~, ~, Tot_Infected,Schooldays_missed , Infected_within_school, Infected_during_term, AsymsCaptured, AsymsTotal, Peak_infect, Tot_Isolating] =   Modeloutputs(history);
       Tot_Infected_5(j,i) = sum(Tot_Infected);
       AsymsCaptured_5(j,i) = sum(AsymsCaptured);
       AsymsTotal_5(j,i) = sum(AsymsTotal);
       Peak_infect_5(j,i) = Peak_infect;
       Tot_Isolating_5(j,i) = mean(Tot_Isolating(:));
       Infected_during_term_5(j,i) = sum(Infected_during_term);
       Infected_within_school_5(j,i) = sum(Infected_within_school);
       Schooldays_missed_5(j,i) = mean(Schooldays_missed(:));
       [Incidence, Num_lat, Num_PCR, Num_conf, School_Inc, Num_lat_t] = Moremodeloutputs(history);
       Num_lat_5(j,i) = Num_lat;       
     
       %no testing
       history = Interactingyeargroups(notestingparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
       [~, ~, ~, ~, ~, Tot_Infected,Schooldays_missed , Infected_within_school, Infected_during_term, AsymsCaptured, AsymsTotal, Peak_infect, Tot_Isolating] =   Modeloutputs(history);
       Tot_Infected_6(j,i) = sum(Tot_Infected);
       AsymsCaptured_6(j,i) = sum(AsymsCaptured);
       AsymsTotal_6(j,i) = sum(AsymsTotal);
       Peak_infect_6(j,i) = Peak_infect;
       Tot_Isolating_6(j,i) = mean(Tot_Isolating(:));
       Infected_during_term_6(j,i) = sum(Infected_during_term);
       Infected_within_school_6(j,i) = sum(Infected_within_school);
       Schooldays_missed_6(j,i) = mean(Schooldays_missed(:));
       
   end
   toc
   
   j
end


B = 0:1:90;
 


figure; 
set(gcf,'Position',[300,300,600,400] );
p = mean(Tot_Infected_1');
q = mean(Tot_Infected_2');
r = mean(Tot_Infected_3');
s = mean(Tot_Infected_4');
t = mean(Tot_Infected_5');
u = mean(Tot_Infected_6');


plot(B, u/10, 'linewidth', 1.5, 'color', [0.5, 0.5, 0.5]); hold on
above =  quantile(100*Tot_Infected_6'/1000, 0.75);
below =  quantile(100*Tot_Infected_6'/1000, 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.5, 0.5, 0.5];

plot(B, q/10, 'linewidth', 1.5, 'color', [0.30, 0.75, 0.93]); hold on
above =  quantile(100*Tot_Infected_2'/1000, 0.75);
below =  quantile(100*Tot_Infected_2'/1000, 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.30, 0.75, 0.93];

p = mean(Tot_Infected_1');
plot(B, p/10, 'linewidth', 1.5, 'color', [1, 0.41, 0.16]); hold on
above =  quantile(100*Tot_Infected_1'/1000, 0.75);
below =  quantile(100*Tot_Infected_1'/1000, 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [1.00, 0.41, 0.16];

 
plot(B, r/10, 'linewidth', 1.5, 'color', [0.72, 0.27, 1.00]); hold on
above =  quantile(100*Tot_Infected_3'/1000, 0.75);
below =  quantile(100*Tot_Infected_3'/1000, 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.72, 0.27, 1.00]; 




plot(B, s/10, 'linewidth', 1.5, 'color', [0.39, 0.83, 0.07]); hold on
above =  quantile(100*Tot_Infected_4'/1000, 0.75);
below =  quantile(100*Tot_Infected_4'/1000, 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.39, 0.83, 0.07];

plot(B, t/10, 'linewidth', 1.5, 'color', [0.93, 0.69, 0.13]); hold on
above =  quantile(100*Tot_Infected_5'/1000, 0.75);
below =  quantile(100*Tot_Infected_5'/1000, 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.93, 0.69, 0.13];

xlabel('Initial % of pupils immune (R_{init})');
ylabel('Total infected by end of half term(%)');

ylim([0 40]);
%xlim([0 70]);
set(gca, 'fontsize', 14);
legend('no control', 'serial contact testing', 'isolating year groups', 'twice weekly mass tests',  'twice weekly mass tests + serial contact testing', 'twice weekly mass tests + isolating year groups');



figure; %Figure S8b
set(gcf,'Position',[300,300,600,400] );
p = mean(Schooldays_missed_1');
q = mean(Schooldays_missed_2');
r = mean(Schooldays_missed_3');
s = mean(Schooldays_missed_4');
t = mean(Schooldays_missed_5');
u = mean(Schooldays_missed_6');

plot(B, r, 'linewidth', 1.5, 'color', [0.72, 0.27, 1]); hold on
above =  quantile(Schooldays_missed_3', 0.75);
below =  quantile(Schooldays_missed_3', 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.72, 0.27, 1];
plot(B, q, 'linewidth', 1.5, 'color', [0.30, 0.75, 0.93]); hold on
above =  quantile(Schooldays_missed_2', 0.75);
below =  quantile(Schooldays_missed_2', 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.30, 0.75, 0.93];

p = mean(Schooldays_missed_1');
plot(B, p, 'linewidth', 1.5, 'color', [1, 0.41, 0.16]); hold on
above =  quantile(Schooldays_missed_1', 0.75);
below =  quantile(Schooldays_missed_1', 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [1.00, 0.41, 0.16];

plot(B, s, 'linewidth', 1.5, 'color', [0.39, 0.83, 0.07]); hold on
above =  quantile(Schooldays_missed_4', 0.75);
below =  quantile(Schooldays_missed_4', 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.39, 0.83, 0.07];

plot(B, t, 'linewidth', 1.5, 'color', [0.93, 0.69, 0.13]); hold on
above =  quantile(Schooldays_missed_5', 0.75);
below =  quantile(Schooldays_missed_5', 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.93, 0.69, 0.13];


plot(B, u, 'linewidth', 1.5, 'color', [0.5, 0.5, 0.5]); hold on
above =  quantile(Schooldays_missed_6', 0.75);
below =  quantile(Schooldays_missed_6', 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.5, 0.5, 0.5];
xlabel('Initial % of pupils immune (R_{init})');
ylabel('Mean school days missed per pupil');
%xlim([0 70]);
set(gca, 'fontsize', 14);

%}

%% Figure S9 Impact of K %%
clearvars -except runs PCR_test_sym PCR_test_asym lat_test_sym lat_test_asym 
  for i = 1:runs

        isolationparams(i,1) = 1 +4*rand;       
        isolationparams(i,2) = 0.001026404983400;
        isolationparams(i,3) = 1;
        isolationparams(i,4) = 1;
        isolationparams(i,5) = 1;
        isolationparams(i,6) = 0.02;
        isolationparams(i,7) = 0.2;
        isolationparams(i,8) = 0;
        isolationparams(i,9) = 0.4*rand + 0.3;
        isolationparams(i,10) = 0.12 + 0.19*rand;
        isolationparams(i,11) = 0;
        isolationparams(i,12) = 200;
        isolationparams(i,13) = 5;
        isolationparams(i,14) = 0;
        isolationparams(i,15) = 0;

  end
    
for j = 1:101
   
    %parameters for isolation
    isolationparams(:,1) = (j-1)*0.1;

    %parameters for serial contact testing
    testingbaselineparams = isolationparams;
    testingbaselineparams(:,3) = 0;
    
    %parameters for weekly mass testing
    biweeklyparams = testingbaselineparams;
    biweeklyparams(:,4) = 6;
    
    %parameters for combined
    biweeklyscparams = testingbaselineparams;
    biweeklyscparams(:,4) = 8;
    
        
    testisolparams = testingbaselineparams;
    testisolparams(:,4) = 14;
        
    notestingparams = testingbaselineparams;
    notestingparams(:,4) = 4;
    
   tic
   parfor i = 1:runs
        %isolation
        history = Interactingyeargroups(isolationparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
       [~, ~, ~, ~, ~, Tot_Infected,Schooldays_missed , Infected_within_school, Infected_during_term, AsymsCaptured, AsymsTotal, Peak_infect, Tot_Isolating] =   Modeloutputs(history);
       Tot_Infected_1(j,i) = sum(Tot_Infected);
       AsymsCaptured_1(j,i) = sum(AsymsCaptured);
       AsymsTotal_1(j,i) = sum(AsymsTotal);
       Peak_infect_1(j,i) = Peak_infect;
       Tot_Isolating_1(j,i) = mean(Tot_Isolating(:));
       Infected_during_term_1(j,i) = sum(Infected_during_term);
       Infected_within_school_1(j,i) = sum(Infected_within_school);
       Schooldays_missed_1(j,i) = mean(Schooldays_missed(:));
       
        %baseline testing
        history = Interactingyeargroups(testingbaselineparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
       [~, ~, ~, ~, ~, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, Peak_infect, Tot_Isolating] =   Modeloutputs(history);
       Tot_Infected_2(j,i) = sum(Tot_Infected);
       AsymsCaptured_2(j,i) = sum(AsymsCaptured);
       AsymsTotal_2(j,i) = sum(AsymsTotal);
       Peak_infect_2(j,i) = Peak_infect;
       Tot_Isolating_2(j,i) = mean(Tot_Isolating(:));
       Schooldays_missed_2(j,i) = mean(Schooldays_missed(:));
       [Incidence, Num_lat, Num_PCR, Num_conf, School_Inc, Num_lat_t] = Moremodeloutputs(history);
       Num_lat_2(j,i) = Num_lat;


       
       
       %weekly mass testing
        history = Interactingyeargroups(biweeklyparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
       [~, ~, ~, ~, ~, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, Peak_infect, Tot_Isolating] =   Modeloutputs(history);
       Tot_Infected_3(j,i) = sum(Tot_Infected);
       AsymsCaptured_3(j,i) = sum(AsymsCaptured);
       AsymsTotal_3(j,i) = sum(AsymsTotal);
       Peak_infect_3(j,i) = Peak_infect;
       Tot_Isolating_3(j,i) = mean(Tot_Isolating(:));
       Schooldays_missed_3(j,i) = mean(Schooldays_missed(:));
       [Incidence, Num_lat, Num_PCR, Num_conf, School_Inc, Num_lat_t] = Moremodeloutputs(history);
       Num_lat_3(j,i) = Num_lat;

        %weekly mass testing and serial contact testing
        history = Interactingyeargroups(biweeklyscparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
       [~, ~, ~, ~, ~, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, Peak_infect, Tot_Isolating] =   Modeloutputs(history);
       Tot_Infected_4(j,i) = sum(Tot_Infected);
       AsymsCaptured_4(j,i) = sum(AsymsCaptured);
       AsymsTotal_4(j,i) = sum(AsymsTotal);
       Peak_infect_4(j,i) = Peak_infect;
       Tot_Isolating_4(j,i) = mean(Tot_Isolating(:));
       Schooldays_missed_4(j,i) = mean(Schooldays_missed(:));
       [Incidence, Num_lat, Num_PCR, Num_conf, School_Inc, Num_lat_t] = Moremodeloutputs(history);
       Num_lat_4(j,i) = Num_lat;

       
       
        % testing and isol
       history = Interactingyeargroups(testisolparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
       [~, ~, ~, ~, ~, Tot_Infected,Schooldays_missed , Infected_within_school, Infected_during_term, AsymsCaptured, AsymsTotal, Peak_infect, Tot_Isolating] =   Modeloutputs(history);
       Tot_Infected_5(j,i) = sum(Tot_Infected);
       AsymsCaptured_5(j,i) = sum(AsymsCaptured);
       AsymsTotal_5(j,i) = sum(AsymsTotal);
       Peak_infect_5(j,i) = Peak_infect;
       Tot_Isolating_5(j,i) = mean(Tot_Isolating(:));
       Infected_during_term_5(j,i) = sum(Infected_during_term);
       Infected_within_school_5(j,i) = sum(Infected_within_school);
       Schooldays_missed_5(j,i) = mean(Schooldays_missed(:));
       [Incidence, Num_lat, Num_PCR, Num_conf, School_Inc, Num_lat_t] = Moremodeloutputs(history);
       Num_lat_5(j,i) = Num_lat;

       
     
       %no testing
       history = Interactingyeargroups(notestingparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
       [~, ~, ~, ~, ~, Tot_Infected,Schooldays_missed , Infected_within_school, Infected_during_term, AsymsCaptured, AsymsTotal, Peak_infect, Tot_Isolating] =   Modeloutputs(history);
       Tot_Infected_6(j,i) = sum(Tot_Infected);
       AsymsCaptured_6(j,i) = sum(AsymsCaptured);
       AsymsTotal_6(j,i) = sum(AsymsTotal);
       Peak_infect_6(j,i) = Peak_infect;
       Tot_Isolating_6(j,i) = mean(Tot_Isolating(:));
       Infected_during_term_6(j,i) = sum(Infected_during_term);
       Infected_within_school_6(j,i) = sum(Infected_within_school);
       Schooldays_missed_6(j,i) = mean(Schooldays_missed(:));
       
   end
   toc
   
   j
end

B = 0:0.1:10;
 
figure; %Figure S9A
set(gcf,'Position',[300,300,600,400] );
p = mean(Tot_Infected_1');
q = mean(Tot_Infected_2');
r = mean(Tot_Infected_3');
s = mean(Tot_Infected_4');
t = mean(Tot_Infected_5');
u = mean(Tot_Infected_6');


plot(B, u/10, 'linewidth', 1.5, 'color', [0.5, 0.5, 0.5]); hold on
above =  quantile(100*Tot_Infected_6'/1000, 0.75);
below =  quantile(100*Tot_Infected_6'/1000, 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.5, 0.5, 0.5];

plot(B, q/10, 'linewidth', 1.5, 'color', [0.30, 0.75, 0.93]); hold on
above =  quantile(100*Tot_Infected_2'/1000, 0.75);
below =  quantile(100*Tot_Infected_2'/1000, 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.30, 0.75, 0.93];

p = mean(Tot_Infected_1');
plot(B, p/10, 'linewidth', 1.5, 'color', [1, 0.41, 0.16]); hold on
above =  quantile(100*Tot_Infected_1'/1000, 0.75);
below =  quantile(100*Tot_Infected_1'/1000, 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [1.00, 0.41, 0.16];

 
plot(B, r/10, 'linewidth', 1.5, 'color', [0.72, 0.27, 1.00]); hold on
above =  quantile(100*Tot_Infected_3'/1000, 0.75);
below =  quantile(100*Tot_Infected_3'/1000, 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.72, 0.27, 1.00]; 

plot(B, s/10, 'linewidth', 1.5, 'color', [0.39, 0.83, 0.07]); hold on
above =  quantile(100*Tot_Infected_4'/1000, 0.75);
below =  quantile(100*Tot_Infected_4'/1000, 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.39, 0.83, 0.07];

plot(B, t/10, 'linewidth', 1.5, 'color', [0.93, 0.69, 0.13]); hold on
above =  quantile(100*Tot_Infected_5'/1000, 0.75);
below =  quantile(100*Tot_Infected_5'/1000, 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.93, 0.69, 0.13];

xlabel('Within-school transmission, K');
ylabel('Total infected by end of half term(%)');

set(gca, 'fontsize', 14);
legend('no control', 'serial contact testing', 'isolating year groups', 'twice weekly mass tests',  'twice weekly mass tests + serial contact testing', 'twice weekly mass tests + isolating year groups');



figure; %Figure S9b
set(gcf,'Position',[300,300,600,400] );
p = mean(Schooldays_missed_1');
q = mean(Schooldays_missed_2');
r = mean(Schooldays_missed_3');
s = mean(Schooldays_missed_4');
t = mean(Schooldays_missed_5');
u = mean(Schooldays_missed_6');

plot(B, r, 'linewidth', 1.5, 'color', [0.72, 0.27, 1]); hold on
above =  quantile(Schooldays_missed_3', 0.75);
below =  quantile(Schooldays_missed_3', 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.72, 0.27, 1];
plot(B, q, 'linewidth', 1.5, 'color', [0.30, 0.75, 0.93]); hold on
above =  quantile(Schooldays_missed_2', 0.75);
below =  quantile(Schooldays_missed_2', 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.30, 0.75, 0.93];

p = mean(Schooldays_missed_1');
plot(B, p, 'linewidth', 1.5, 'color', [1, 0.41, 0.16]); hold on
above =  quantile(Schooldays_missed_1', 0.75);
below =  quantile(Schooldays_missed_1', 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [1.00, 0.41, 0.16];

plot(B, s, 'linewidth', 1.5, 'color', [0.39, 0.83, 0.07]); hold on
above =  quantile(Schooldays_missed_4', 0.75);
below =  quantile(Schooldays_missed_4', 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.39, 0.83, 0.07];

plot(B, t, 'linewidth', 1.5, 'color', [0.93, 0.69, 0.13]); hold on
above =  quantile(Schooldays_missed_5', 0.75);
below =  quantile(Schooldays_missed_5', 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.93, 0.69, 0.13];


plot(B, u, 'linewidth', 1.5, 'color', [0.5, 0.5, 0.5]); hold on
above =  quantile(Schooldays_missed_6', 0.75);
below =  quantile(Schooldays_missed_6', 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.5, 0.5, 0.5];
xlabel('Within-school transmission, K')
ylabel('Mean school days missed per pupil');
%xlim([0 70]);
set(gca, 'fontsize', 14);


clearvars -except runs PCR_test_sym PCR_test_asym lat_test_sym lat_test_asym 

%% Figures S11 and S12 -Impact of cohort size

%Currently set to density dependent transmission, for S12, change
%isolationparams(:,15)
  for i = 1:runs

        isolationparams(i,1) = 1 +4*rand;       
        isolationparams(i,2) = 0.001026404983400;
        isolationparams(i,3) = 1;
        isolationparams(i,4) = 1;
        isolationparams(i,5) = 1;
        isolationparams(i,6) = 0.02;
        isolationparams(i,7) = 0.2;
        isolationparams(i,8) = 0;
        isolationparams(i,9) = 0.4*rand + 0.3;
        isolationparams(i,10) = 0.12 + 0.19*rand;
        isolationparams(i,11) = 0;
        isolationparams(i,12) = 200;
        isolationparams(i,13) = 5;
        isolationparams(i,14) = 0;
        isolationparams(i,15) = 0;

    end
tic


for j = 1:4
    
    j
    
    if j == 1        
        isolationparams(:,12) = 25;
        isolationparams(:,13) = 40;
    elseif j == 2
        isolationparams(:,12) = 50;
        isolationparams(:,13) = 20;
    elseif j == 3
        isolationparams(:,12) = 100;
        isolationparams(:,13) = 10;
    elseif j == 4
        isolationparams(:,12) = 200;
        isolationparams(:,13) = 5;
    end
    
        isolationparams(:,3) = 1; 
        isolationparams(:,14) = 0;    
        isolationparams(:,15) = 1; %density
        %isolationparams(:,15) = 0; %frequency dependent transmission
    
    
    testingbaselineparams = isolationparams;
    testingbaselineparams(:,3) = 0;
    
    testinglowparams = testingbaselineparams;
    testinglowparams(:,4) = 2;
    testinghighparams = testingbaselineparams;
    testinghighparams(:,4) = 3;
    
    notestingparams = testingbaselineparams;
    notestingparams(:,4) = 4;
    
    biweeklyparams = testingbaselineparams;
    biweeklyparams(:,4) = 6;
      
    biweeklyscparams = testingbaselineparams;
    biweeklyscparams(:,4) = 8;
    
    testisolparams = testingbaselineparams;
    testisolparams(:,4) = 14;
    
    parfor i = 1:runs
           %baseline testing
                   
                   
                   
            history = Interactingyeargroups(isolationparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
            [Prev_true, Prev_asym, Known_asym, pos_tests, Absences, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, ~, ~, Days_tested, PresymsCaptured, day_caught_vec] = Modeloutputs(history); 
            Tot_Infected_1(i,j) = sum(Tot_Infected);
            Schooldays_missed_1(i,j) = sum(Schooldays_missed(:));
            [Incidence, Num_lat, Num_PCR, Num_conf, School_Inc, Num_lat_t] = Moremodeloutputs(history);
            Num_lat_1(i,j) = NaN;


            
            history = Interactingyeargroups(testingbaselineparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
            [Prev_true, Prev_asym, Known_asym, pos_tests, Absences, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, ~, ~, Days_tested, PresymsCaptured, day_caught_vec] = Modeloutputs(history); 
            Tot_Infected_2(i,j) = sum(Tot_Infected);
            Schooldays_missed_2(i,j) = sum(Schooldays_missed(:));
            [Incidence, Num_lat, Num_PCR, Num_conf, School_Inc, Num_lat_t] = Moremodeloutputs(history);
            Num_lat_2(i,j) = Num_lat;
                   
            
            %weekly mass testing 
            history = Interactingyeargroups(biweeklyparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
            [Prev_true, Prev_asym, Known_asym, pos_tests, Absences, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, ~, ~, Days_tested, PresymsCaptured, day_caught_vec] = Modeloutputs(history); 
            Tot_Infected_6(i,j) = sum(Tot_Infected);
            Schooldays_missed_6(i,j) = sum(Schooldays_missed(:));
            [Incidence, Num_lat, Num_PCR, Num_conf, School_Inc, Num_lat_t] = Moremodeloutputs(history);
            Num_lat_6(i,j) = Num_lat;
            
            
            %weekly mass testing + SCT
             history = Interactingyeargroups(biweeklyscparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
            [Prev_true, Prev_asym, Known_asym, pos_tests, Absences, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, ~, ~, Days_tested, PresymsCaptured, day_caught_vec] = Modeloutputs(history); 
            Tot_Infected_8(i,j) = sum(Tot_Infected);
            Schooldays_missed_8(i,j) = sum(Schooldays_missed(:));
            [Incidence, Num_lat, Num_PCR, Num_conf, School_Inc, Num_lat_t] = Moremodeloutputs(history);
            Num_lat_8(i,j) = Num_lat;
            
            %biweekly mass testing + SCT
            history = Interactingyeargroups(testisolparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
            [Prev_true, Prev_asym, Known_asym, pos_tests, Absences, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, ~, ~, Days_tested, PresymsCaptured, day_caught_vec] = Modeloutputs(history); 
            Tot_Infected_9(i,j) = sum(Tot_Infected);
            Schooldays_missed_9(i,j) = sum(Schooldays_missed(:));
            [Incidence, Num_lat, Num_PCR, Num_conf, School_Inc, Num_lat_t] = Moremodeloutputs(history);
            Num_lat_9(i,j) = Num_lat;
          
            history = Interactingyeargroups(notestingparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
            [Prev_true, Prev_asym, Known_asym, pos_tests, Absences, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, ~, ~, Days_tested, PresymsCaptured, day_caught_vec] = Modeloutputs(history); 
            Tot_Infected_5(i,j) = sum(Tot_Infected);
            Schooldays_missed_5(i,j) = sum(Schooldays_missed(:)); 
            [Incidence, Num_lat, Num_PCR, Num_conf, School_Inc, Num_lat_t] = Moremodeloutputs(history);
            Num_lat_5(i,j) = NaN;
            
            
    end

    
    
    
end
toc

figure;
A = [100*(Tot_Infected_1/1000),   NaN*ones(runs,1), 100*(Tot_Infected_9/1000), NaN*ones(runs,1), 100*(Tot_Infected_2/1000), NaN*ones(runs,1), 100*(Tot_Infected_8/1000), NaN*ones(runs,1), 100*(Tot_Infected_6/1000), NaN*ones(runs,1), 100*(Tot_Infected_5/1000)]';
    plot(100, 100, '.', 'MarkerSize', 20, 'color', [0.9, 0.9, 0.9]); hold on
    plot(100, 100, '.', 'MarkerSize', 20, 'color', [0.7, 0.7, 0.7]);
    plot(100, 100, '.', 'MarkerSize', 20, 'color', [0.5, 0.5, 0.5]); hold on
    plot(100, 100, '.', 'MarkerSize', 20, 'color', [0, 0, 0]);

v = violinplot(A'); hold on


xticks([2.5, 7.5, 12.5, 17.5, 22.5, 27.5]);
 xticklabels({'isolating cohorts', 'mass testing + isol', 'serial contact testing', 'mass testing + SCT', 'mass testing', 'no control'});

for i = 1:size(A,1)
   if mod(i,5) == 1
           v(i).ViolinColor = [0.9, 0.9, 0.9];
   elseif mod(i,5) == 2
       v(i).ViolinColor = [0.7, 0.7, 0.7];
   elseif mod(i,5) == 3
       v(i).ViolinColor = [0.5, 0.5, 0.5];
   elseif mod(i,5) == 4
       v(i).ViolinColor = [0, 0, 0];
   end
end

set(gcf, 'Position', [300, 300, 1000, 450]);
set(gca, 'fontsize', 14);
xlabel('Reopening strategy');
ylabel('Total infected by end of half term (%)');
legend('Group size 25', 'Group size 50', 'Group size 100', 'Group size 200');
     
ylim([0 60]);


figure;
A = [(Schooldays_missed_1/1000),   NaN*ones(runs,1), Schooldays_missed_9/1000, NaN*ones(runs,1), Schooldays_missed_2/1000, NaN*ones(runs,1), Schooldays_missed_8/1000, NaN*ones(runs,1), Schooldays_missed_6/1000, NaN*ones(runs,1), Schooldays_missed_5/1000]';
    plot(100, 100, '.', 'MarkerSize', 20, 'color', [0.9, 0.9, 0.9]); hold on
    plot(100, 100, '.', 'MarkerSize', 20, 'color', [0.7, 0.7, 0.7]);
    plot(100, 100, '.', 'MarkerSize', 20, 'color', [0.5, 0.5, 0.5]); hold on
    plot(100, 100, '.', 'MarkerSize', 20, 'color', [0, 0, 0]);

v = violinplot(A'); hold on


xticks([2.5, 7.5, 12.5, 17.5, 22.5, 27.5]);
 xticklabels({'isolating cohorts', 'mass testing + isol', 'serial contact testing', 'mass testing + SCT', 'mass testing', 'no control'});

for i = 1:size(A,1)
   if mod(i,5) == 1
           v(i).ViolinColor = [0.9, 0.9, 0.9];
   elseif mod(i,5) == 2
       v(i).ViolinColor = [0.7, 0.7, 0.7];
   elseif mod(i,5) == 3
       v(i).ViolinColor = [0.5, 0.5, 0.5];
   elseif mod(i,5) == 4
       v(i).ViolinColor = [0, 0, 0];
   end
end

set(gcf, 'Position', [300, 300, 1000, 450]);
set(gca, 'fontsize', 14);
xlabel('Reopening strategy');
ylabel('Mean school days missed per pupil');
legend('Group size 25', 'Group size 50', 'Group size 100', 'Group size 200', 'Location', 'northwest');
ylim([0 25]);

    axes('Position', [.58 .58 .3 .3])
     box on
    v = violinplot(A(11:end,:)');  
for i = 1:size(A(11:end,:),1)
   if mod(i,5) == 1
           v(i).ViolinColor = [0.9, 0.9, 0.9];
   elseif mod(i,5) == 2
       v(i).ViolinColor = [0.7, 0.7, 0.7];
   elseif mod(i,5) == 3
       v(i).ViolinColor = [0.5, 0.5, 0.5];
   elseif mod(i,5) == 4
       v(i).ViolinColor = [0, 0, 0];
   end
end 

    xticks([2.5, 7.5, 12.5, 17.5]);

   xticklabels({'serial contact testing', 'mass testing + SCT', 'mass testing', 'no control'});     

       set(gcf, 'Position', [300, 300, 1000, 450]);
            set(gca, 'fontsize', 14);

            %{
figure;
A = [(Num_lat_1/1000),   NaN*ones(runs,1), (Num_lat_9/1000), NaN*ones(runs,1), (Num_lat_2/1000), NaN*ones(runs,1), (Num_lat_8/1000), NaN*ones(runs,1), (Num_lat_6/1000), NaN*ones(runs,1), (Num_lat_5/1000)]';
    plot(100, 100, '.', 'MarkerSize', 20, 'color', [0.9, 0.9, 0.9]); hold on
    plot(100, 100, '.', 'MarkerSize', 20, 'color', [0.7, 0.7, 0.7]);
    plot(100, 100, '.', 'MarkerSize', 20, 'color', [0.5, 0.5, 0.5]); hold on
    plot(100, 100, '.', 'MarkerSize', 20, 'color', [0, 0, 0]);

v = violinplot(A'); hold on


xticks([2.5, 7.5, 12.5, 17.5, 22.5, 27.5]);
 xticklabels({'', 'mass testing + isol', 'serial contact testing', 'mass testing + SCT', 'mass testing', ''});

for i = 1:size(A,1)
   if mod(i,5) == 1
           v(i).ViolinColor = [0.9, 0.9, 0.9];
   elseif mod(i,5) == 2
       v(i).ViolinColor = [0.7, 0.7, 0.7];
   elseif mod(i,5) == 3
       v(i).ViolinColor = [0.5, 0.5, 0.5];
   elseif mod(i,5) == 4
       v(i).ViolinColor = [0, 0, 0];
   end
end

xlim([0 30]);

xlabel('Reopening strategy');
ylabel('Mean number of LFTs per pupil');
legend('Group size 25', 'Group size 50', 'Group size 100', 'Group size 200');
ylim([0 35]);
set(gcf, 'Position', [300, 300, 1000, 450]);
set(gca, 'fontsize', 14);


%}

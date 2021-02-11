%This .m file produces Supplementary Figure S5


%Positive test profiles
    PCR_test_sym = readtable('PCR_Curve_summary.csv');
    PCR_test_sym = table2array(PCR_test_sym(:, 2:4));
    PCR_test_asym = csvread('PCR_Curve_asym.csv');

    lat_test_sym = readtable('lat_Curve_summary.csv');
    lat_test_sym = table2array(lat_test_sym(:, 2:4));
    lat_test_asym = csvread('lat_Curve_asym.csv');

%{    
%Find external infection values to satisfy 10% total infected
for j = 1:10:41
    K = 0.1*(j-1);
    g = @(x) (interactingyeargroups_k(K,x, PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym) - 0.1);
    
    if j == 1
        x0 = 0.0025;
    else
        x0 = xmin_totinf(j-10);
    end
    
    xmin_totinf(j) = fzero(g, x0);
    j
end
%}
    
%generate parameters 
rng(1);
runs = 100;
for i = 1:runs
    isolationparams(i,1) = 0; % to change for each
    isolationparams(i,2) = 0; % to change
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


%Vary K from 0 to 4
for j = 1:41
   
    %parameters for isolation
    isolationparams(:,1) = (j-1)*0.1;
    isolationparams(:,2) = xmin_totinf_interp(j);
    
    %parameters for serial contact testing
    testingbaselineparams = isolationparams;
    testingbaselineparams(:,3) = 0;
    
    %parameters for weekly mass testing
    weeklyparams = testingbaselineparams;
    weeklyparams(:,4) = 5;
    
    %parameters for combined
    weeklyscparams = testingbaselineparams;
    weeklyscparams(:,4) = 7;  
    
   tic
   for i = 1:runs
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
       
       
       %weekly mass testing
        history = Interactingyeargroups(weeklyparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
       [~, ~, ~, ~, ~, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, Peak_infect, Tot_Isolating] =   Modeloutputs(history);
       Tot_Infected_3(j,i) = sum(Tot_Infected);
       AsymsCaptured_3(j,i) = sum(AsymsCaptured);
       AsymsTotal_3(j,i) = sum(AsymsTotal);
       Peak_infect_3(j,i) = Peak_infect;
       Tot_Isolating_3(j,i) = mean(Tot_Isolating(:));
       Schooldays_missed_3(j,i) = mean(Schooldays_missed(:));

        %weekly mass testing and serial contact testing
        history = Interactingyeargroups(weeklyscparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
       [~, ~, ~, ~, ~, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, Peak_infect, Tot_Isolating] =   Modeloutputs(history);
       Tot_Infected_4(j,i) = sum(Tot_Infected);
       AsymsCaptured_4(j,i) = sum(AsymsCaptured);
       AsymsTotal_4(j,i) = sum(AsymsTotal);
       Peak_infect_4(j,i) = Peak_infect;
       Tot_Isolating_4(j,i) = mean(Tot_Isolating(:));
       Schooldays_missed_4(j,i) = mean(Schooldays_missed(:));       
       
   end
   toc
   
   j
end
%}   
%proportion of infections occurring within school under isolation strategy   
B = 100*mean(Infected_within_school_1')./mean(Infected_during_term_1');
    
figure; %Figure S5C
set(gcf,'Position',[300,300,600,400] );
p = mean(Tot_Infected_1');
q = mean(Tot_Infected_2');
r = mean(Tot_Infected_3');
s = mean(Tot_Infected_4');
 
plot(B, r/10, 'linewidth', 1.5, 'color', [0.72, 0.27, 1.00]); hold on
above =  quantile(100*Tot_Infected_3'/1000, 0.75);
below =  quantile(100*Tot_Infected_3'/1000, 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.72, 0.27, 1.00]; 

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

plot(B, s/10, 'linewidth', 1.5, 'color', [0.39, 0.83, 0.07]); hold on
above =  quantile(100*Tot_Infected_4'/1000, 0.75);
below =  quantile(100*Tot_Infected_4'/1000, 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.39, 0.83, 0.07];

xlabel('% of infections within school under isolation strategy')
ylabel('Total Infected (%)')
xlim([0 70]);
set(gca, 'fontsize', 12);
legend( 'weekly mass testing', 'serial contact testing', 'isolation of year group bubbles', 'combined');


figure; %Figure S5D
set(gcf,'Position',[300,300,600,400] );
p = mean(Peak_infect_1');
q = mean(Peak_infect_2');
r = mean(Peak_infect_3');
s = mean(Peak_infect_4');

plot(B, r/10, 'linewidth', 1.5, 'color', [0.72, 0.27, 1.00]); hold on
above =  quantile(100*Peak_infect_3'/1000, 0.75);
below =  quantile(100*Peak_infect_3'/1000, 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.72, 0.27, 1.00]; 

plot(B, q/10, 'linewidth', 1.5, 'color', [0.30, 0.75, 0.93]); hold on
above =  quantile(100*Peak_infect_2'/1000, 0.75);
below =  quantile(100*Peak_infect_2'/1000, 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.30, 0.75, 0.93];

p = mean(Peak_infect_1');
plot(B, p/10, 'linewidth', 1.5, 'color', [1, 0.41, 0.16]); hold on
above =  quantile(100*Peak_infect_1'/1000, 0.75);
below =  quantile(100*Peak_infect_1'/1000, 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [1.00, 0.41, 0.16];

plot(B, s/10, 'linewidth', 1.5, 'color', [0.39, 0.83, 0.07]); hold on
above =  quantile(100*Peak_infect_4'/1000, 0.75);
below =  quantile(100*Peak_infect_4'/1000, 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.39, 0.83, 0.07];

xlabel('% of infections within school under isolation strategy')
ylabel('Peak prevalence (%)')
xlim([0 70]);
set(gca, 'fontsize', 12);
%legend( 'weekly mass testing', 'serial contact testing', 'isolation of year group bubbles', 'combined');


figure; %Figure S5E
set(gcf,'Position',[300,300,600,400] );
q = sum(AsymsCaptured_2')./sum(AsymsTotal_2');
r = sum(AsymsCaptured_3')./sum(AsymsTotal_3');
s = sum(AsymsCaptured_4')./sum(AsymsTotal_4');

plot(B, 100*r, 'linewidth', 1.5, 'color', [0.72, 0.27, 1.00]); hold on
above =  quantile(100*(AsymsCaptured_3'./AsymsTotal_3'), 0.75);
below =  quantile(100*(AsymsCaptured_3'./AsymsTotal_3'), 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.72, 0.27, 1.00];

plot(B, 100*q, 'linewidth', 1.5, 'color', [0.30, 0.75, 0.93]); hold on
above =  quantile(100*(AsymsCaptured_2'./AsymsTotal_2'), 0.75);
below =  quantile(100*(AsymsCaptured_2'./AsymsTotal_2'), 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.30, 0.75, 0.93];

plot(B, 100*s, 'linewidth', 1.5, 'color', [0.39, 0.83, 0.07]); hold on
above =  quantile(100*(AsymsCaptured_4'./AsymsTotal_4'), 0.75);
below =  quantile(100*(AsymsCaptured_4'./AsymsTotal_4'), 0.25);
p = patch([B fliplr(B)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
p.FaceColor = [0.39, 0.83, 0.07];

xlabel('% of infections within school under isolation strategy')
ylabel('% of asymptomatics identified');
xlim([0 70]);
set(gca, 'fontsize', 12);
%legend( 'weekly mass testing', 'serial contact testing', 'combined');


figure; %Figure S5F
set(gcf,'Position',[300,300,600,400] );
p = mean(Schooldays_missed_1');
q = mean(Schooldays_missed_2');
r = mean(Schooldays_missed_3');
s = mean(Schooldays_missed_4');

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
xlabel('% of infections within school under isolation strategy')
ylabel('Average days isolating per student');
xlim([0 70]);
set(gca, 'fontsize', 12);
%legend( 'weekly mass testing', 'serial contact testing', 'isolation of year group bubbles', 'combined');



%Results through time at K = 0 and at K = 4
for j = [1 41]
    isolationparams(:,1) = (j-1)*0.1;
    isolationparams(:,2) = xmin_totinf_interp(j);
    
    testingbaselineparams = isolationparams;
    testingbaselineparams(:,3) = 0;
    
    notestingparams = testingbaselineparams;
    notestingparams(:,4) = 4;
    
    weeklyparams = testingbaselineparams;
    weeklyparams(:,4) = 5;
    
    weeklyscparams = testingbaselineparams;
    weeklyscparams(:,4) = 7;
    
   for i = 1:runs
        %isolation
        history = Interactingyeargroups(isolationparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
       Prev_true =   Modeloutputs(history);
        Prevs_1(i,:) = sum(Prev_true);


        %baseline testing
        history = Interactingyeargroups(testingbaselineparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
        Prev_true =   Modeloutputs(history);
        Prevs_2(i,:) = sum(Prev_true);
       
       
       %weekly mass testing
        history = Interactingyeargroups(weeklyparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
       Prev_true =   Modeloutputs(history);
       Prevs_3(i,:) = sum(Prev_true);



        %weekly mass testing and serial contact testing
        history = Interactingyeargroups(weeklyscparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
       Prev_true =   Modeloutputs(history);
       Prevs_4(i,:) = sum(Prev_true);

       
       
   end
   
     figure; %Figures 2A and 2B
     set(gcf, 'Position', [300, 300, 600, 400]); 
     plot((1:56)/7 - 1, 100*mean(Prevs_3)/1000, 'linewidth', 1.5, 'color', [0.72, 0.27, 1.00]); hold on
     above =  quantile(100*Prevs_3/1000, 0.75);
     below =  quantile(100*Prevs_3/1000, 0.25);
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
     plot((1:56)/7 - 1, 100*mean(Prevs_4)/1000, 'linewidth', 1.5, 'color', [0.39, 0.83, 0.07]); hold on
     above =  quantile(100*Prevs_4/1000, 0.75);
     below =  quantile(100*Prevs_4/1000, 0.25);
     p = patch([(1:56)/7 - 1 fliplr((1:56)/7 - 1)], [below above(end:-1:1)], 'b', 'FaceAlpha',0.25, 'EdgeAlpha', 0, 'HandleVisibility', 'off'); hold on
     p.FaceColor = [0.39, 0.83, 0.07];
     xlim([(1/7)-1 7]);
     plot([0 0], [0.5, 4.5], 'k--');
     legend('no control', 'weekly mass testing', 'serial contact testing', 'Isolating year groups', 'combined', 'no control');
     xlabel('Week');
     ylabel('Prevalence (%)');
     set(gca, 'fontsize', 12);
   
    
end




function Z = interactingyeargroups_k(K,Ext, PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym)
%Interacting year groups model averaged over 1000 runs, with K and Ext as parameters       
    runs = 1000;    
    isolationparams = zeros(runs,11); 
    %generate parameters
    
    rng(1);
    for i = 1:runs
        isolationparams(i,1) = K;
        isolationparams(i,2) = Ext;
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
    
    for i = 1:runs
        history = Interactingyeargroups(isolationparams(i,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
        [~, ~, ~, ~, ~, Tot_Infected] = Modeloutputs(history);       
        Tot_Inf(i) = sum(Tot_Infected);        
        
    end

    Z = mean(Tot_Inf)/1000;


end
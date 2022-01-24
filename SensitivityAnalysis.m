%This code generates the figures for the sensitivity analysis (Supplementary Figures S3-S6)

clear

%Positive test profiles
    PCR_test_sym = readtable('PCR_Curve_summary.csv');
    PCR_test_sym = table2array(PCR_test_sym(:, 2:4));
    PCR_test_asym = csvread('PCR_Curve_asym.csv');

    lat_test_sym = readtable('lat_Curve_summary.csv');
    lat_test_sym = table2array(lat_test_sym(:, 2:4));
    lat_test_asym = csvread('lat_Curve_asym.csv');
    
    runs = 20;  %For quick run
   % runs = 2000; %Number of simulations for paper
    
    %set baseline parameters
    isolationparams = zeros(16, 11);
    
       isolationparams(:,1) = 3;
      isolationparams(:,2) = 0.001026404983400;      
        isolationparams(:,3) = 1;
        isolationparams(:,4) = 1;
        isolationparams(:,5) = 1;
        isolationparams(:,6) = 0.02;
        isolationparams(:,7) = 0.2;
        isolationparams(:,8) = 0;
        isolationparams(:,9) = 0.5;
        isolationparams(:,10) = 0.215;
        isolationparams(:,11) = 0;
        isolationparams(:,12) = 200;
        isolationparams(:,13) = 5;
        isolationparams(:,14) = 0;
        isolationparams(:,15) = 0;

    
    testingbaselineparams = isolationparams;
    testingbaselineparams(:,3) = 0;

    
    notestingparams = testingbaselineparams;
    notestingparams(:,4) = 4;
    
    biweeklyparams = testingbaselineparams;
    biweeklyparams(:,4) = 6;
    
    
    
    biweeklyscparams = testingbaselineparams;
    biweeklyscparams(:,4) = 8;
    
    
    testisolparams = testingbaselineparams;
    testisolparams(:,4) = 14;
    
    
    
    %sensitivity analyses
    
    %(2 and 3): lateral flow tests, 2 - low, 3 - high
    isolationparams(3,4) = 2;
    isolationparams(2,4) = 3; 

    %(4 and 5): PCR test positivity, 4 - low, 5 - high
    isolationparams(5,5) = 2;
    isolationparams(4,5) = 3;
    
    %(6 and 7) within school infection K, 6 - K = 1, 7 - K = 4
    isolationparams(7,1) = 1;
    isolationparams(6,1) = 5;
    %isolationparams(6,1) = 10;
    
    %(8 and 9) Population level immunity, 8 - 0.1%, 9 - 30%
    isolationparams(9, 7) = 0.1;
    isolationparams(8,7) = 0.3;
    
    %(10 and 11) %symptomatic, 10 - 12% symptomatics, 11 - 31% symptomatic
    isolationparams(11, 10) = 0.12;
    isolationparams(10, 10) = 0.31;
    
    %(12 and 13) infectiousness of asymptomatics 12 - 30% as infectious,
    %13 - 70% as infectious  
    isolationparams(13, 9) = 0.3;
    isolationparams(12, 9) = 0.7;
    
    %(14 and 15) Population level transmission 14 - half as much, 15 -
    %twice as much
    isolationparams(15,2) = 0.5*isolationparams(1,2);
    isolationparams(14,2) = 2*isolationparams(1,2);
    
    %16 - interaction between year groups
    isolationparams(16, 11) = 1;
            
        
    testingbaselineparams = isolationparams;
    testingbaselineparams(:,3) = 0;

    
    notestingparams = testingbaselineparams;
    notestingparams(:,4) = 4;
    
    biweeklyparams = testingbaselineparams;
    biweeklyparams(:,4) = 6;
    biweeklyparams(2,4) = 10; %for weekly mass testing
    biweeklyparams(3,4) = 9; %for weekly mass testing
    
   
    biweeklyscparams = testingbaselineparams;
    biweeklyscparams(:,4) = 8;
    biweeklyscparams(2,4) = 12; %for combined
    biweeklyscparams(3,4) = 11; %for combined
    
    testisolparams = testingbaselineparams;
    testisolparams(:,4) = 14;
    testisolparams(2,4) = 16;
    testisolparams(3,4) = 15;
    
    for j = 1:16
        j
        
        tic
        parfor i = 1:runs
            
            %For isolation of year groups
            
            history = Interactingyeargroups(isolationparams(j,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
            [~,~,~,~,~, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, ~]  = Modeloutputs(history);

            Tot_Infected_1(j, i) = sum(Tot_Infected);
            Schooldays_missed_1(j, i) = sum(Schooldays_missed(:));
            AsymsCaptured_1(j,i) = sum(AsymsCaptured);
            AsymsTotal_1(j,i) = sum(AsymsTotal);
            
                        [Incidence, Num_lat, Num_PCR, Num_conf, School_Inc, Num_lat_t] = Moremodeloutputs(history);
            Num_lat_1(j,i) = 0;
            
            %For rapid testing strategy
            history = Interactingyeargroups(testingbaselineparams(j,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
            [~,~,~,~,~, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, ~]  = Modeloutputs(history);

            Tot_Infected_2(j, i) = sum(Tot_Infected);
            Schooldays_missed_2(j, i) = sum(Schooldays_missed(:));
            AsymsCaptured_2(j,i) = sum(AsymsCaptured);
            AsymsTotal_2(j,i) = sum(AsymsTotal);  
            
            [Incidence, Num_lat, Num_PCR, Num_conf, School_Inc, Num_lat_t] = Moremodeloutputs(history);
            Num_lat_2(j,i) = Num_lat;
            
            
            %For twice weekly mass testing
            history = Interactingyeargroups(biweeklyparams(j,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
            [~,~,~,~,~, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, ~]  = Modeloutputs(history);

            Tot_Infected_6(j, i) = sum(Tot_Infected);
            Schooldays_missed_6(j, i) = sum(Schooldays_missed(:));
            AsymsCaptured_6(j,i) = sum(AsymsCaptured);
            AsymsTotal_6(j,i) = sum(AsymsTotal);
            
           [Incidence, Num_lat, Num_PCR, Num_conf, School_Inc, Num_lat_t] = Moremodeloutputs(history);
            Num_lat_6(j,i) = Num_lat;
            
           %For twice weekly mass testing + SCT
            history = Interactingyeargroups(biweeklyscparams(j,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
            [~,~,~,~,~, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, ~]  = Modeloutputs(history);

            Tot_Infected_8(j, i) = sum(Tot_Infected);
            Schooldays_missed_8(j, i) = sum(Schooldays_missed(:));
            AsymsCaptured_8(j,i) = sum(AsymsCaptured);
            AsymsTotal_8(j,i) = sum(AsymsTotal);  
            
            [Incidence, Num_lat, Num_PCR, Num_conf, School_Inc, Num_lat_t] = Moremodeloutputs(history);
            Num_lat_8(j,i) = Num_lat;
            
            %For twice weekly mass testing + isolation
            history = Interactingyeargroups(testisolparams(j,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
            [~,~,~,~,~, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, ~]  = Modeloutputs(history);

            Tot_Infected_9(j, i) = sum(Tot_Infected);
            Schooldays_missed_9(j, i) = sum(Schooldays_missed(:));
            AsymsCaptured_9(j,i) = sum(AsymsCaptured);
            AsymsTotal_9(j,i) = sum(AsymsTotal);  
            
            [Incidence, Num_lat, Num_PCR, Num_conf, School_Inc, Num_lat_t] = Moremodeloutputs(history);
            Num_lat_9(j,i) = Num_lat;
                       
        end
        
        toc
    end
    
    
    
    schools_more_cases = [];
    increase_in_cases = [];
    Asyms_known_prop = [];
    days_missed_reduction = [];
    
    for j = 1:16
        
        %mean Total infected                       
        Tot_inf_diff(j) = mean(Tot_Infected_1(j,:));        
        Tot_inf_diff9(j) = mean(Tot_Infected_9(j,:));
        Tot_inf_diff2(j) = mean(Tot_Infected_2(j,:));        
        Tot_inf_diff8(j) = mean(Tot_Infected_8(j,:));
        Tot_inf_diff6(j) = mean(Tot_Infected_6(j,:));        
        
        %mean days missed per pupil
        days_missed1(j) = mean(Schooldays_missed_1(j,:))/1000;
        days_missed9(j) = mean(Schooldays_missed_9(j,:))/1000;
        days_missed2(j) = mean(Schooldays_missed_2(j,:))/1000;
        days_missed8(j) = mean(Schooldays_missed_8(j,:))/1000;
        days_missed6(j) = mean(Schooldays_missed_6(j,:))/1000;
        
 
        %LFTs taken per pupil
        lfts1(j) = mean(Num_lat_1(j,:))/1000;
        lfts2(j) = mean(Num_lat_2(j,:))/1000;
        lfts6(j) = mean(Num_lat_6(j,:))/1000;
        lfts8(j) = mean(Num_lat_8(j,:))/1000;
        lfts9(j) = mean(Num_lat_9(j,:))/1000;
        
        %{
        asymprop2(j) = sum(AsymsCaptured_2(j,:))./sum(AsymsTotal_2(j,:));
        asymprop6(j) = sum(AsymsCaptured_6(j,:))./sum(AsymsTotal_6(j,:));
        asymprop8(j) = sum(AsymsCaptured_8(j,:))./sum(AsymsTotal_8(j,:));
        asymprop9(j) = sum(AsymsCaptured_9(j,:))./sum(AsymsTotal_9(j,:));
        %}

    end
   %} 

   
   cdata = [1.00, 0.41, 0.16; 0.30, 0.75, 0.93; 0.72, 0.27, 1.00; 0.39, 0.83, 0.07; 0.93, 0.69, 0.13]';
   
   
Tot_Inf_collective = [Tot_inf_diff; Tot_inf_diff2; Tot_inf_diff6; Tot_inf_diff8; Tot_inf_diff9]/10;

days_missed_collective = [days_missed1; days_missed2; days_missed6; days_missed8; days_missed9];

lfts_collective = [lfts1; lfts2; lfts6; lfts8;lfts9];

for i = 1:16
%temp = days_missed_collective(:,i); %To do days missed
%temp = Tot_Inf_collective(:,i); %To do total infected
temp = lfts_collective(:,i); %To do lfts
temp = [temp, (1:5)'];
temp = sortrows(temp);
infvals(:,i) = temp(:,1);
inforders(:,i) = temp(:,2);
end

inforderslong = inforders(:);
infvalslong = infvals(:);

counter = 1;
combinedlong = zeros(1,2);
for i = 1:length(inforderslong)
    if ~any(combinedlong(:,1) == infvalslong(i))
       combinedlong(counter,:) = [infvalslong(i), inforderslong(i)]; 
       counter = counter+1;
    end
end




%% Code to make colour data for the heat map (so that colours show the strategy)
%This is done by making a colour matrix with colour specifying strategy
%Check colours are unique in each row
combinedlong = sortrows(combinedlong);
%}

combinedverylong = zeros(round(20000*(combinedlong(end,1) - combinedlong(1,1))),2);

counter = 1;
combinedverylong(1,1) =  (combinedlong(1,1) + i/20000);
combinedverylong(1,2) = combinedlong(counter,2);  

for i = 2:length(combinedverylong)
    combinedverylong(i,1) =  (combinedlong(1,1) + i/20000);
    combinedverylong(i,2) = combinedlong(counter,2);    
    if  ~((combinedlong(1,1) + i/20000) < combinedlong(counter+1,1)  || i == length(combinedverylong))
        counter = counter+1;
        combinedverylong(i,2) = combinedlong(counter,2);        
    end
   
end
cdata3 = [];
for i = 1:length(combinedverylong)
    cdata3(:,i) = cdata(:,combinedverylong(i,2));
end

H = heatmap(infvals');
H.Colormap = cdata3';
H.CellLabelFormat = '%.2f';
H.ColorbarVisible = 'off';
ax = gca;
ax.XData = ["Lowest" "2nd lowest" "3rd lowest" "4th lowest" "Highest"];
ax.YData = ["Baseline", "high LFT sensitivity", "low LFT sensitivity", "high PCR sensitivity", "low PCR sensitivity", "high within-school transmission ", "low within-school transmission", "high immunity ", "low immunity", "high % symptomatic", "low % symptomatic" , "high infectiousness of asymptomatics" , "low infectiousness of asymptomatics", "high probability of community infection", "low probability of community infection", "interaction between year groups"];
set(gca, 'fontsize', 14);
xlabel('Total infected by end of half term (%)');
%xlabel('Mean school days missed per pupil');
xlabel('Mean number of LFTs per pupil');
%ylabel('Sensitivity scenario');
set(gcf, 'Position', [300, 300, 1000, 800]);

%For Isol
Tot_inf_diff_now =  Tot_inf_diff/10;
days_missednow = days_missed1;
lftsnow = lfts1;
%}

%For SCT
%{
Tot_inf_diff_now =  Tot_inf_diff2/10;
days_missednow = days_missed2;
lftsnow = lfts2;
%}

%For twice weekly mass testing
%{
Tot_inf_diff_now =  Tot_inf_diff6/10;
days_missednow = days_missed6;
lftsnow = lfts6;
%}

%For mass testing + SCT
%{
Tot_inf_diff_now =  Tot_inf_diff8/10;
days_missednow = days_missed8;
lftsnow = lfts8;
%}

%For mass testing + isol
%{
Tot_inf_diff_now =  Tot_inf_diff9/10;
days_missednow = days_missed9;
lftsnow = lfts9;
%}

   figure;  
  % set(gcf,'Position',[300,300,700,400]);
  set(gcf,'Position',[300,300,1000,250]);
  pos = [0.2 0.2 0.25 0.7];
  subplot('Position',pos); hold on
  box on;

    barh(fliplr(Tot_inf_diff_now(2:2:end)), 'BaseValue', Tot_inf_diff_now(1), 'FaceColor', '#edf8b1', 'FaceAlpha', 0.75); hold on
    barh([ Tot_inf_diff_now(1), fliplr(Tot_inf_diff_now(3:2:end))], 'BaseValue', Tot_inf_diff_now(1), 'FaceColor', '#2c7fb8' , 'FaceAlpha', 0.75);
   yticks([0:9]);
    yticklabels({ '',  'interaction between year groups', 'community infection', 'infectiousness of asymptomatics',  '% symptomatic',  'population immunity', 'within school transmission', 'PCR sensitivity', 'LFT sensitivity'})

   %yticklabels({});
    ylim([0 9]);
    xlim([4 16]);
    plot([50 50], [0 10], 'k--');
    set(gca, 'fontsize', 12);

   xlabeldat = {'Total infected by end of half term (%)'};

   xlabel(xlabeldat);

    

  pos = [0.47 0.2 0.25 0.7];
  subplot('Position',pos); hold on
  box on;

   barh(fliplr(days_missednow(2:2:end)), 'BaseValue', days_missednow(1), 'FaceColor', '#edf8b1', 'FaceAlpha', 0.75); hold on
   barh([days_missednow(1) fliplr(days_missednow(3:2:end))], 'BaseValue', days_missednow(1), 'FaceColor', '#2c7fb8', 'FaceAlpha', 0.75);
  yticks([0:9]);
   
   yticklabels({});
   ylim([0 9]);
   
   xlim([0 20]);
   set(gca, 'fontsize', 12);
   legend('higher than baseline', 'lower than baseline');
    xlabeldat = {'Mean school days missed per pupil'};

   xlabel(xlabeldat);
   
   
   pos = [0.74 0.2 0.25 0.7];
  subplot('Position',pos); hold on
      box on;   
   barh(fliplr(lftsnow(2:2:end)), 'BaseValue', lftsnow(1), 'FaceColor', '#edf8b1' , 'FaceAlpha', 0.75); hold on
   barh([lftsnow(1) fliplr(lftsnow(3:2:end))], 'BaseValue', lftsnow(1), 'FaceColor', '#2c7fb8' , 'FaceAlpha', 0.75);
    yticks([0:9]);
   
   yticklabels({});
   ylim([0 9]);
    xlim([5 35]);
   set(gca, 'fontsize', 12);   
  xlabeldat = {'Mean number of LFTs per pupil'};
   xlabel(xlabeldat);

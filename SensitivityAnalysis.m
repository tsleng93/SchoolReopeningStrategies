%This code generates the figures for the sensitivity analysis (Supplementary Figures S2-S4)
%N.B. to generate the sensitivity analysis for weekly mass testing or
%combined testing, uncomment the relevant parts of the script


%Positive test profiles
    PCR_test_sym = readtable('PCR_Curve_summary.csv');
    PCR_test_sym = table2array(PCR_test_sym(:, 2:4));
    PCR_test_asym = csvread('PCR_Curve_asym.csv');

    lat_test_sym = readtable('lat_Curve_summary.csv');
    lat_test_sym = table2array(lat_test_sym(:, 2:4));
    lat_test_asym = csvread('lat_Curve_asym.csv');
    
    runs = 100;  
    
    %set baseline parameters
    isolationparams = zeros(16, 11);
    
    isolationparams(:,1) = 2.5;
    isolationparams(:,2) = 0.0013133;
    isolationparams(:,3) = 1;
    isolationparams(:,4) = 1; %for serial contact testing
    %isolationparams(:,4) = 5; %for weekly mass testing
    %isolationparams(:,4) = 7; %for combined
    
    isolationparams(:,5) = 1;
    isolationparams(:,6) = 0.02;
    isolationparams(:,7) = 0.2;
    isolationparams(:,8) = 0;
    isolationparams(:,9) = 0.5;
    isolationparams(:,10) = 0.215;
    isolationparams(:,11) = 0;
    
    %sensitivity analyses
    
    %(2 and 3): lateral flow tests, 2 - low, 3 - high
    isolationparams(2,4) = 2; %for serial contact testing
    isolationparams(3,4) = 3; %for serial contact testing
    %isolationparams(2,4) = 9; %for weekly mass testing
    %isolationparams(3,4) = 10; %for weekly mass testing
    %isolationparams(2,4) = 11; %for combined
    %isolationparams(3,4) = 12; %for combined
    
    
    %(4 and 5): PCR test positivity, 4 - low, 5 - high
    isolationparams(4,5) = 2;
    isolationparams(5,5) = 3;
    
    %(6 and 7) within school infection K, 6 - K = 1, 7 - K = 4
    isolationparams(6,1) = 1;
    isolationparams(7,1) = 4;
    
    %(8 and 9) Population level immunity, 8 - 10%, 9 - 30%
    isolationparams(8, 7) = 0.1;
    isolationparams(9,7) = 0.3;
    
    %(10 and 11) %symptomatic, 10 - 12% symptomatics, 11 - 31% symptomatic
    isolationparams(10, 10) = 0.12;
    isolationparams(11, 10) = 0.31;
    
    %(12 and 13) infectiousness of asymptomatics 12 - 30% as infectious,
    %13 - 70% as infectious  
    isolationparams(12, 9) = 0.3;
    isolationparams(13, 9) = 0.7;
    
    %(14 and 15) Population level transmission 14 - half as much, 15 -
    %twice as much
    isolationparams(14,2) = 0.5*isolationparams(14,2);
    isolationparams(15,2) = 2*isolationparams(15,2);
    
    %16 - interaction between year groups
    isolationparams(16, 11) = 1;
    
    testingparams = isolationparams;
    testingparams(:,3) = 0;

    
    
    
    for j = 1:16
        j
        
        tic
        for i = 1:runs
            
            %For isolation of year groups
            history = Interactingyeargroups(isolationparams(j,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
            [~,~,~,~,~, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, ~]  = Modeloutputs(history);

            Tot_Infected_1(j, i) = sum(Tot_Infected);
            Schooldays_missed_1(j, i) = sum(Schooldays_missed(:));
            AsymsCaptured_1(j,i) = sum(AsymsCaptured);
            AsymsTotal_1(j,i) = sum(AsymsTotal);
            
            %For rapid testing strategy
            history = Interactingyeargroups(testingparams(j,:), PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, i);
            [~,~,~,~,~, Tot_Infected, Schooldays_missed, ~, ~, AsymsCaptured, AsymsTotal, ~]  = Modeloutputs(history);

            Tot_Infected_2(j, i) = sum(Tot_Infected);
            Schooldays_missed_2(j, i) = sum(Schooldays_missed(:));
            AsymsCaptured_2(j,i) = sum(AsymsCaptured);
            AsymsTotal_2(j,i) = sum(AsymsTotal);          
            
            
        end
        
        toc
    end
    
    %}
    schools_more_cases = [];
    increase_in_cases = [];
    Asyms_known_prop = [];
    days_missed_reduction = [];
    
    for j = 1:16
        %schools with more cases
        schools_more_cases(j) = sum(Tot_Infected_2(j,:) > Tot_Infected_1(j,:));
        
        %increase in infections
        increase_in_cases(j) = (sum(Tot_Infected_2(j,:)) - sum(Tot_Infected_1(j,:)))/sum(Tot_Infected_1(j,:));
        
        %Asymptomatics
        Asyms_known_prop(j) = sum(AsymsCaptured_2(j,:))/sum(AsymsTotal_2(j,:));
        
        %Reduction in school days missed 
        days_missed_reduction(j) = 1 - sum(Schooldays_missed_2(j,:))/sum(Schooldays_missed_1(j,:));
    end
   %} 
   
   %{
   schools_more_cases(end+1) = schools_more_cases(end);
   increase_in_cases(end+1) = increase_in_cases(end);
   Asyms_known_prop(end+1) = Asyms_known_prop(end);
   days_missed_reduction(end+1) = days_missed_reduction(end);
   
   schools_more_cases(end-1) = schools_more_cases(end);
   increase_in_cases(end-1) = increase_in_cases(end);
   Asyms_known_prop(end-1) = Asyms_known_prop(end);
   days_missed_reduction(end-1) = days_missed_reduction(end);
%}

   
   figure;  
  % set(gcf,'Position',[300,300,700,400]);
  set(gcf,'Position',[300,300,1000,800]);
  pos = [0.25 0.55 0.3 0.4];
  subplot('Position',pos); hold on
  box on;
   %ax1 = axes('Position',[0.33 0.2 0.6 0.8]);
   %ax1.ActivePositionProperty = 'position';
    barh(100*schools_more_cases(2:2:end)/runs, 'BaseValue', 100*schools_more_cases(1)/runs, 'FaceColor', [0.49, 0.18, 0.56]); hold on
    barh(100*schools_more_cases(3:2:end)/runs, 'BaseValue', 100*schools_more_cases(1)/runs, 'FaceColor', [0.93, 0.69, 0.13]);
   yticklabels({ '', 'LFT sensitivity','PCR sensitivity','within school transmission', 'population immunity', '% symptomatic', 'infectiousness of asymptomatics', 'community infection', 'interaction between year groups'})
   %yticklabels({});
    ylim([0 9]);
    xlim([0 100]);
    plot([50 50], [0 10], 'k--');
    set(gca, 'fontsize', 10);
    legend('lower than baseline', 'higher than baseline');
    xlabel('more cases with serial contact testing testing (%)');
    %xlabel('more cases with weekly mass testing (%)'); % for weekly mass testing
    %xlabel('more cases with combined testing (%)'); % for combined testing

    
   %figure;
  pos = [0.6 0.55 0.3 0.4];
  subplot('Position',pos); hold on
  box on;
   %set(gcf,'Position',[300,300,700,400]);
   barh(100*increase_in_cases(2:2:end), 'BaseValue', 100*increase_in_cases(1), 'FaceColor', [0.49, 0.18, 0.56]); hold on
   barh(100*increase_in_cases(3:2:end), 'BaseValue', 100*increase_in_cases(1), 'FaceColor', [0.93, 0.69, 0.13]);
   yticklabels({});
   ylim([0 9]);
   xlim([-30 40]);
   plot([0 0], [0 10], 'k--');
   %yticklabels({'LFT sensitivity','PCR sensitivity','within school transmission', 'population immunity', '% symptomatic', 'infectiousness of asymptomatics', 'community infection', 'interaction between year groups'})
   set(gca, 'fontsize', 10);
   %legend('lower than baseline', 'higher than baseline');
   xlabel('relative increase in total infected (%)');
   
     
     pos = [0.25 0.05 0.3 0.4];
  subplot('Position',pos); hold on
     box on;
      %set(gcf,'Position',[300,300,700,400]);
   %ax1 = axes('Position',[0.33 0.2 0.6 0.8]);
   %ax1.ActivePositionProperty = 'position';
    barh(100*Asyms_known_prop(2:2:end), 'BaseValue', 100*Asyms_known_prop(1), 'FaceColor', [0.49, 0.18, 0.56]); hold on
    barh(100*Asyms_known_prop(3:2:end), 'BaseValue', 100*Asyms_known_prop(1), 'FaceColor', [0.93, 0.69, 0.13]);
   yticklabels({'', 'LFT sensitivity','PCR sensitivity','within school transmission', 'population immunity', '% symptomatic', 'infectiousness of asymptomatics', 'community infection', 'interaction between year groups'})
   %yticklabels({});
    ylim([0 9]);
    xlim([10 80]);
    set(gca, 'fontsize', 10);
    %legend('lower than baseline', 'higher than baseline');
    xlabel('% asymptomatics identified');
    
   
   pos = [0.6 0.05 0.3 0.4];
  subplot('Position',pos); hold on
      box on;   
   barh(100*days_missed_reduction(2:2:end), 'BaseValue', 100*days_missed_reduction(1), 'FaceColor', [0.49, 0.18, 0.56]); hold on
   barh(100*days_missed_reduction(3:2:end), 'BaseValue', 100*days_missed_reduction(1), 'FaceColor', [0.93, 0.69, 0.13]);
   yticklabels({});
   ylim([0 9]);
   xlim([90 100]);
  % yticklabels({'LFT sensitivity','PCR sensitivity','within school transmission', 'population immunity', '% symptomatic', 'infectiousness of asymptomatics', 'community infection', 'interaction between year groups'})
   set(gca, 'fontsize', 10);
   %legend('lower than baseline', 'higher than baseline');
   xlabel('relative reduction in absences (%)');
   
    


  
    
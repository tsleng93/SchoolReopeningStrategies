function Fittingmodelfunction

%ABC method used for fitting the vection

%Initialise
Pillar_data_weeks = [];
Pillar_data_hist_2020 = [];
Pillar_data_5_weeks = [];
Pillar_data_hist_2021 = [];
Adherence = [];
Strategy = [];
Testing =  [];

%load('Dataforfitting.mat'); %%Commented out because of data availability



%load relevant testing and infection profiles
PCR_test_sym = readtable('PCR_Curve_summary.csv');
PCR_test_sym = table2array(PCR_test_sym(:, 2:4));
PCR_test_asym = csvread('PCR_Curve_asym.csv');
lat_test_sym = readtable('lat_Curve_summary.csv');
lat_test_sym = table2array(lat_test_sym(:, 2:4));
lat_test_asym = csvread('lat_Curve_asym.csv');    

Infectivity_since_infection = [0.0063 0.0563 0.1320 0.1798 0.1817 0.1521 0.1117 0.0746 0.0464 0.0272 0.0152 0.0082 0.0043 0.0022 0.0011 0.0005 0.0002 0.0001 0.0001 0.0000];
Infectivity_since_infection = Infectivity_since_infection'/sum(Infectivity_since_infection);
Infectivity_since_infection(15) = sum(Infectivity_since_infection(15:end)); Infectivity_since_infection = Infectivity_since_infection(1:15);

Symptom_onset =  [0.0055 0.0534 0.1307 0.1814 0.1847 0.1545 0.1129 0.0747 0.0458 0.0265 0.0146 0.0077 0.0039 0.0020 0.0010 0.0005 0.0002 0.0001 0 0];
Symptom_onset = Symptom_onset/sum(Symptom_onset);
Symptom_onset(15) = sum(Symptom_onset(15:end)); Symptom_onset = Symptom_onset(1:15);

Prob_profiles.PCRsym = PCR_test_sym;
Prob_profiles.PCRasym = PCR_test_asym;
Prob_profiles.latsym = lat_test_sym;
Prob_profiles.latasym = lat_test_asym;
Prob_profiles.Infectivity = Infectivity_since_infection;
Prob_profiles.Symptom_onset = Symptom_onset;





params = zeros(100,12);
LLs = zeros(100,1);




%bounds for parameter values
paramsmin = [1 1 1 1 0 0 0 0 0 1 0.999 1 ];
paramsmax = [2 4 2 3 1 1 1 1 1.5 2 1 4];

%Set function
numschools = 100;
LLfun2 = @(infparams,data_t,data_hist, data_t2, data_hist2)LLfunLTLAwithEaster(infparams, data_t, data_hist, Prob_profiles, Testing, Adherence, Strategy, combinedschools(:,4), combinedschools(:,3), combinedschools(:,2), LTLAEasterHoliday, communityprevforschool, numtestingposforschool, numtestingnegforschool, numpopLTLA, data_t2, data_hist2, schoolLTLA, Sigmoidfit, initsLTLA, UrbanRuralType2, combinedschools(:,8), schoolregion, numschools);


%Initial runs
LLprevgens = [];
paramsandLLkeep = [];
tic
parfor i = 1:100
    i
    
    %Initial search - force model to find finite values
    LLs(i) = Inf;
    
    while isinf(LLs(i)) 
        
    p = zeros(1,12);    
    
    p(1) = 2*rand + 1; %K
    p(2) = 3*rand + 1; %Ext
    p(3) = rand + 1; %After first half-term
    p(4) = 2*rand+1; %Impact of Holidays
    p(5) = rand; %rural or urban
    p(6) = rand; %Infectiousness of asymptomatics
    p(7) = rand; %Proportion symptomatic
    p(8) = rand; %sigma parameter lognormal Ext
    p(9) = rand; %sigma parameter lognormal K
    p(10) = 1 + 1*rand; %new variant impact
    p(11) = 0.999 + 0.001*rand; %LFT specificity
    p(12) = 3*rand +1; %LFT underreporting
    
    params(i,:) = p;

    
   % p = params(i,:);
   
   
   
    if any(p < paramsmin) || any(p > paramsmax)
        LLs(i) = Inf; %%Infinite log-likelihood if out of bounds
    else        

    LLs(i) = LLfun2(p, Pillar_data_weeks, Pillar_data_hist_2020, Pillar_data_5_weeks, Pillar_data_hist_2021); %record log likelihood   
    end
    
    
    end
end
toc
median(LLs(~isinf(LLs))) %median finite LLs from initial runs



%Middle runs - with 100 sampled schools
paramsandLL = [];
 for j = 1:9
     

     
     LLfun2 = @(infparams,data_t,data_hist, data_t2, data_hist2)LLfunLTLAwithEaster(infparams, data_t, data_hist, Prob_profiles, Testing, Adherence, Strategy, combinedschools(:,4), combinedschools(:,3), combinedschools(:,2), LTLAEasterHoliday, communityprevforschool, numtestingposforschool, numtestingnegforschool, numpopLTLA, data_t2, data_hist2, schoolLTLA, Sigmoidfit, initsLTLA, UrbanRuralType2, combinedschools(:,8), schoolregion, numschools);

    LLprevgens = [LLprevgens, LLs];    
    paramsandLLkeep = [paramsandLLkeep; paramsandLL];    
    paramsandLL = [params, LLs];
    paramsandLL = sortrows(paramsandLL, size(params, 2)+1);
    
    
    %Keep top 20%    
    params = zeros(100,size(params,2));
    params(1:20,:) = paramsandLL(1:20,1:size(params,2));    
    LLs(1:20) = paramsandLL(1:20, size(params,2)+1);
    
   
    %Sample rest from covariance    
    Covs = cov(paramsandLL(1:20,1:size(params,2)));
   params(21:end,:) =  mvnrnd(mean(paramsandLL(1:20,1:size(params,2))),2*Covs,80);   
   meanp = mean(paramsandLL(1:20,1:size(params,2)));
   
       %parpool(4);
       tic
      parfor i = 1:100 

          p = params(i,:);

            while any(p < paramsmin) || any(p > paramsmax)
               %choose new parameters if out of bounds
               params(i,:) = mvnrnd(meanp,2*Covs,1);
               p = params(i,:);                         
            end       

          LLs(i) = LLfun2(p, Pillar_data_weeks, Pillar_data_hist_2020, Pillar_data_5_weeks, Pillar_data_hist_2021);

      end
       
      toc
      median(LLs(~isinf(LLs)))

 end
 %}

 save('JulXfit');


%Runs on final size
minvalue = 0;
cutoff = 0.05;
threshold  = 1;
medianLL = median(LLs);

LLtemp = sort(LLs);

LLthresh = (LLtemp(20));

%longer runs - with 1000 sampled schools
for j = 1:20
    
    
    numschools = 1000;
    
   LLprevgens = [LLprevgens, LLs(1:100)];
    
   paramsandLLkeep = [paramsandLLkeep; paramsandLL];
    
    paramsandLL = [params, LLs];
    paramsandLL = sortrows(paramsandLL, size(params, 2)+1);
    
    %Keep top 20%    
    params = zeros(100,size(params,2));
    params(1:20,:) = paramsandLL(1:20,1:size(params,2));        
    LLs(1:20) = paramsandLL(1:20, size(params,2)+1);
    
    %Sample rest from covariance
    Covs = cov(paramsandLL(1:20,1:size(params,2)));   
    params(21:end,:) =  mvnrnd(mean(paramsandLL(1:20,1:size(params,2))),2*Covs,80);
    meanp = mean(paramsandLL(1:20,1:size(params,2)));    
    
       LLfun2 = @(infparams,data_t,data_hist, data_t2, data_hist2)LLfunLTLAwithEaster(infparams, data_t, data_hist, Prob_profiles, Testing, Adherence, Strategy, combinedschools(:,4), combinedschools(:,3), combinedschools(:,2), LTLAEasterHoliday, communityprevforschool, numtestingposforschool, numtestingnegforschool, numpopLTLA, data_t2, data_hist2, schoolLTLA, Sigmoidfit, initsLTLA, UrbanRuralType2, combinedschools(:,8), schoolregion, numschools);


        
        tic
        parfor i = 1:100 

          p = params(i,:);

            while any(p < paramsmin) || any(p > paramsmax)
               params(i,:) = mvnrnd(meanp,2*Covs,1);
               p = params(i,:);                         
            end       

          LLs(i) = LLfun2(p, Pillar_data_weeks, Pillar_data_hist_2020, Pillar_data_5_weeks, Pillar_data_hist_2021);

        end
        toc    
end


save('JulXfit');

 


end

%Likelihood function for model fitting
 function [LL, tests_model_t, tests_model_t2, tests_model_hist, tests_absences_t, model_poslft, model_posboth, model_numlft, model_prev, model_inc, tests_model_hist1, tests_model_hist2]  = LLfunLTLAwithEaster(infparams, data_t, data_hist, Prob_profiles, Testing, Adherence, Strategy, closecontactsizes, yeargroupsizes, numyeargroups, LTLAEasterHoliday, communityprevforschool, numtestingposforschool, numtestingnegforschool, numpopLTLA, data_t2, data_hist2, schoolLTLA, Sigmoidfit, initsLTLA, urbanorrural, closecontactsizesnewterm, schoolregion, numschools)
   
    %Infection parameters (that don't depend on randomness or LTLA)
    Infection.Weeks = length(data_t)+1; %Number of weeks (including week before term    
    %Infection.alpha_betweenyears = 0.1;%Interaction between year groups
    Infection.alpha_betweenyears = 0.01;
    Infection.alpha_withinyear = 0.1;

    Infection.leak_infect = 0; %0 or 1, whether infection occurs on test days    
    Infection.Rec_0 = 0.0625;
    Infection.Sym_proportion = infparams(7);         
    Infection.K_asym = infparams(6);

    
    %Testing_parameters 
    Testing.sens_PCR = 1; % parameter deciding which column of the PCR sensitivity profile to read
    %(1) - baseline sensitivity, (2) - low sensitivity, (3) - high sensitivity 
    Testing.sens_lat = 1; %parameter deciding which column of the LFT sensitivity profile to read.
    %(1) - baseline sensitivity, (2) - low sensitivity, (3) - high sensitivity 
    Testing.spec_PCR = 1; %specificity of PCR tests
    Testing.spec_lat = infparams(11); %specificity of LFT tests
    Testing.PCR_delay = 2; %delay on PCR tests
    
    %Adherence parameters
    Adherence.C_isolate = 0; % constant scaling the external force of infection on isolating individuals
    Adherence.serialnonadherence = 0; %proportion of individuals who agree to be tested but never actually take home tests
    Adherence.takenproperly = 1;
    Adherence.isolatingproperly = zeros(1,500);

    
    
    %Strategy params (set as isolation for baselines
    Strategy.isolation = 1; %(0) - not isolating, (1) - isolating close contacts, (2) - isolating year groups
    Strategy.SCT = 0; %(0) - no serial contact testing, (1) - serial contact testing close contacts, (2) - serial contact testing year groups %%overrides isolation
    Strategy.SCTuptake = 0; %Proportion of individuals agreed to participate in serial contact testing
    Strategy.masstesting = 2; %(0) - no mass testing, (1) - weekly mass testing, (2) - twice weekly mass testing
    Strategy.masstestinguptake = 1; %Proportion of individuals agreed to be mass tested
    Strategy.initialtestdays = []; %initial testing days
    Strategy.initialuptake = 1; %uptake to initial testing
        

    p_rewire = 0; %rewiring of close contacts (always has been 0)
    
    
    
    %relative frequency of S-gene negatives for school
    for i = 1:length(schoolLTLA)
        communitySgeneforschool(i,:) = 1*(1-Sigmoidfit( 237:(238+ 7*Infection.Weeks), schoolLTLA(i))) + infparams(10)*Sigmoidfit(237:(238+7*Infection.Weeks), schoolLTLA(i));  
    end
    
    
  
  %choose schools
  eees = randsample(length(schoolLTLA), numschools);




              
        for i = 1:numschools
                      
                 eee = eees(i);
              
            
               %Set infection parameters 
               Infection.K = lognrnd(log(infparams(1)) - (infparams(9)^2/2), infparams(9));
               Infection.K2 = infparams(3)*Infection.K;    
               Infection.K3 = Infection.K2;

           
          %assigning weeks of easter holiday dependent on LTLA
          if LTLAEasterHoliday(schoolLTLA(eee)) == 1
             Infection.HolidayWeek = [10, 18:28, 32, 33];
          else
             Infection.HolidayWeek = [10, 18:28, 33, 34]; 
          end



          Infection.Inf_0 = initsLTLA(schoolLTLA(eee)); %assigning initial infection dependent on LTLA

           %external infection impacted by being urban or rural
          if urbanorrural(eee) == 2
            Infection.Ext = lognrnd((log(infparams(2)) - (infparams(8)^2/2)), infparams(8));
          else
            Infection.Ext = infparams(5)*lognrnd((log(infparams(2)) - (infparams(8)^2/2)), infparams(8));      
          end

          Infection.HolidayExt = infparams(4); %scaling of external infection during holidays/schoolclosures 

          Underreporting = ones(1,Infection.Weeks*7); %Underreporting
          Underreporting(30*7:end) = infparams(12);
          Adherence.probtakelatflow = ((Underreporting.*numtestingnegforschool(i,:)) + numtestingposforschool(i,:))/numpopLTLA(i); %LFT participation           


           
          Infection.commext = communityprevforschool(eee,:); %positive testing rates in the community          
          Prob_profiles.newvar = communitySgeneforschool(eee,:); %proportion of new variant
            
         [school_pop, ~] = SchoolPopulation(yeargroupsizes(eee), numyeargroups(eee), closecontactsizes(eee), p_rewire);
         [school_pop2, ~] = SchoolPopulation(yeargroupsizes(eee), numyeargroups(eee), closecontactsizesnewterm(eee), p_rewire);
                    
         history = interactingyeargroupsquicker(school_pop, Infection, Testing, Strategy, Adherence, Prob_profiles, i, school_pop2);                   
         [pos_PCRsday, peakcovabs, Absences, Prev_true, Num_LFTs, Pos_LFTs, Pos_LFTsandPCRs, Inc] = Modeloutputscondensed(history, Infection.Weeks);  
         
         
        
         
           Prevs_t(i,:) = sum(Prev_true);
           Incs_t(i,:) = sum(Inc);
           tests_t(i,:) = pos_PCRsday;
           size_school(i) = yeargroupsizes(eee)*numyeargroups(eee);            
           Peak_covid_absences(i) = max(peakcovabs);
                        
            Rinfs1(i,:) = history.Rday(1,:);
            Rinfs2(i,:) = history.Rday(2,:);
            
            ExternalIncs(i,:) = squeeze(sum(sum((history.ext_or_int == 1))));
            InternalIncs(i,:) = squeeze(sum(sum((history.ext_or_int == 2))));
            
            numlft_t(i,:) = Num_LFTs;
            poslft_t(i,:) = Pos_LFTs;
           posboth_t(i,:) = Pos_LFTsandPCRs;
            
           Absences_t(i,:) = sum(Absences);
           
            Peak_covid_absences1(i) = max(peakcovabs(8:119));
            Peak_covid_absences2(i) = max(peakcovabs(end-76:end));


        end
       
        
        %model outputs
        
        tests_model_t = sum(tests_t)/sum(size_school);
        tests_model_hist = histcounts(Peak_covid_absences, [0:15 Inf]);
        tests_model_hist = tests_model_hist/sum(tests_model_hist);
        
        
        
       tests_model_hist1 = histcounts(Peak_covid_absences1, [0:15 Inf]);
       tests_model_hist2 = histcounts(Peak_covid_absences2, [0:8 Inf]);
       tests_model_hist1 = tests_model_hist1/sum(tests_model_hist1);
       tests_model_hist2 = tests_model_hist2/sum(tests_model_hist2);
                
        tests_absences_t = sum(Absences_t)/sum(size_school);        
        model_prev = sum(Prevs_t)/sum(size_school);
        model_inc = sum(Incs_t)/sum(size_school);        
        model_numlft = sum(numlft_t)/sum(size_school);
        model_poslft = sum(poslft_t)/sum(size_school);        
        model_posboth = sum(posboth_t)/sum(size_school);
        %}
        

 
 num_3 = 3896599; %for secondary schools
 
 tests_model_t2 = model_poslft;
 for j = 1:(Infection.Weeks-1)
     tests_model_t_weeks(j) = sum(tests_model_t(1 +(j-1)*7: 7 + (j-1)*7));
     tests_model_t_weeks2(j) = sum(tests_model_t2(1 +(j-1)*7: 7 + (j-1)*7));

 end
 %}
 
 tests_model_t_weeks2 = tests_model_t_weeks2((end-10):end);
 
 
 
 
 optllt = sum(data_t.*log(data_t/num_3) + (num_3 - data_t).*log(1 - (data_t/num_3)));
 optllt2 = sum(data_t2.*log(data_t2/num_3) + (num_3 - data_t2).*log(1 - (data_t2/num_3))); 
 optllhist = sum(data_hist.*log(data_hist/sum(data_hist)));
 optllhist2 = sum(data_hist2.*log(data_hist2/sum(data_hist2)));
 
 
 %Log likelihood functions
 
 LL_t = sum(data_t.*log(tests_model_t_weeks) + (num_3 - data_t).*log(1 - tests_model_t_weeks)) - optllt;
 LL_t2 = (sum(data_t2.*log(tests_model_t_weeks2) + (num_3 - data_t2).*log(1 - tests_model_t_weeks2)) - optllt2);
 LL_hist1 = (sum(data_hist.*log(tests_model_hist1)) - (optllhist));% 
 LL_hist2 = (sum(data_hist2.*log(tests_model_hist2)) - (optllhist2));%
 LL = -1*(LL_t + LL_hist1 + LL_t2 + LL_hist2);
 

end
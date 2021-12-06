function Fittingmodelfunction

%ABC-SMC method used for fitting the vection


%initialise
Pillar_data_weeks = [];
Pillar_data_hist_2020 = [];
Pillar_data_5_weeks = [];
Pillar_data_hist_2021 = [];
Adherence = [];
Strategy = [];
Testing =  [];

%load('Dataforfitting.mat'); %%Commented out because of data availability

%% ---- load relevant testing and infection profiles
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

%% ---- set initial parameters and data

%min and max parameter values
paramsmin = [0 0 0.5 1 0 0 0 0 1 1 ];
paramsmax = [4 4 3 3 1 1 1 1  2 4];

N = 200; %number of particles
LLprevgens = []; %log likelihoods from previous generations
wprevgens = []; %weightings from previous generations
paramsandLLkeep = []; %params, LLs, and weightings from previous generations
params = zeros(N,10); %parameters (will be filled in in for loop)
LLs = Inf*ones(N,1); %Log likelihoods (start all as Inf)
Numgens = 20;
timespergen = zeros(1,Numgens+1); %number of times loop is repeated each generation
iterspergen = zeros(1,Numgens+1); %number of parameter sets tried each generation

esnow = eees; %subset of schools to consider (100 to begin with)
closecontactsizesnow = combinedschools(esnow,4); %close contact sizes
yeargroupsizesnow = combinedschools(esnow,3); %year group sizes
numyeargroupsnow = combinedschools(esnow,2); %number of year groups
communityprevnow = communityprevforschool(esnow,:); %community testing rates for each school
numtestingposnow = numtestingposforschool(esnow,:); %number of 10-19 year olds testing positive in school's LTLA
numtestingnegnow = numtestingnegforschool(esnow,:); %number of 10-19 year olds testing positive in school's LTLA
numpopLTLAnow = numpopLTLA(esnow); %number of 10-19 year olds in each school's LTLA
schoolLTLAnow = schoolLTLA(esnow); %the LTLA of each school
Urbannow = UrbanRuralType2(esnow); %whether school is rural or urban
closecontactsnewtermnow = combinedschools(esnow,8); %close contact sizes after new term


Lmax = 0; %irrelevant when 100 schools - set to 0 for now
%model Log Likelihood function
LLfun2 = @(infparams,data_t,data_hist, data_t2, data_hist2)LLfunLTLAwithEaster(infparams, data_t, data_hist, Prob_profiles, Testing, Adherence, Strategy, closecontactsizesnow, yeargroupsizesnow, numyeargroupsnow, LTLAEasterHoliday, communityprevnow, numtestingposnow, numtestingnegnow, numpopLTLAnow, data_t2, data_hist2, schoolLTLAnow, Sigmoidfit, initsLTLA, Urbannow, closecontactsnewtermnow, Lmax);


%% ----- Initial run -----
tic
while any(isinf(LLs))   
  
   timespergen(1) = timespergen(1)+1; 
   iterspergen(1) = iterspergen(1) + sum(isinf(LLs));  
   timespergentemp = timespergen(1)
   
   rng(N^3 + timespergentemp*(N));   

    newparams = (paramsmax - paramsmin).*rand(N,10) + paramsmin; %set new parameter
    parfor i = 1:N   

        p = newparams(i,:);

        if isinf(LLs(i)) % only enter into loop if LLs(i) is infinite

            rng(N^3 + timespergentemp*(N) + (i-1));    

            if (any(p < paramsmin) || any(p > paramsmax)) 
                LL = Inf; %out of prior range
            else        
                LL = LLfun2(p, Pillar_data_weeks, Pillar_data_hist_2020, Pillar_data_5_weeks, Pillar_data_hist_2021);
            end

            LLs(i) = LL; %update Log likelihood
            weighting(i) = 1; %set weighting as 1 for initial runs
            params(i,:) = p; %update parameters
        end


    end

end
toc

 %% --- Main model fitting loop -----

 B = 1; %number of runs for each parameter set
 KernelConst = 0.68; %perturbation kernel constant
 paramsandLL = []; %combined parameters, LLs, and weighting

 
  for j = 1:Numgens
    j 
    
    tic
         
     sum(LLs > eps) %number of parameter sets above threshold - will be 0        
     if sum(LLs > eps) == 0
         
        eps = quantile(LLs, 0.5); %threshold for each generation          
        epskeep(j) = eps; %store the threshold of each generation

        LLprevgens = [LLprevgens, LLs]; %store log likelihoods from preveious generations      
        paramsandLLkeep = [paramsandLLkeep; paramsandLL];  %store params, LLs, and weightings  
        paramsandLL = [params, LLs, weighting'/sum(weighting)]; %update params and LLs
    
        w = weighting'/sum(weighting); %normalise weightings
        wprevgens = [wprevgens, w']; %store weightings
        ESS = sum(1/(w.^2)); %calculate effective sample size
        ESSkeep(j) = ESS; %store effective sample sizes
        X = paramsandLL(:, 1:size(params,2)); %X - previous generations accepted params, to be used in parfor loop
        mu = sum(w.*X) / sum(w); % row containing weighted mean of each column
        Kernel = KernelConst*((w/sum(w)).*(X-mu))'*(X - mu)/(1 - sum((w/sum(w)).^2)); %Perturbation Kernel
        LLs = Inf*ones(N,1); %Set LLs to be -Infinity for now 
        

     end
         
    if j > 10
       %sample of 600 schools
       esnow = eees600;
       Lmax = epskeep(11); %epsilon value after 10 generations at smaller sample
    else
       %sample of 100 schools
       esnow = eees;
       Lmax = 0; %not used 
    end
        
    closecontactsizesnow = combinedschools(esnow,4); %close contact sizes
    yeargroupsizesnow = combinedschools(esnow,3); %year group sizes
    numyeargroupsnow = combinedschools(esnow,2); %number of year groups
    communityprevnow = communityprevforschool(esnow,:); %community testing rates for each school
    numtestingposnow = numtestingposforschool(esnow,:); %number of 10-19 year olds testing positive in school's LTLA
    numtestingnegnow = numtestingnegforschool(esnow,:); %number of 10-19 year olds testing positive in school's LTLA
    numpopLTLAnow = numpopLTLA(esnow); %number of 10-19 year olds in each school's LTLA
    schoolLTLAnow = schoolLTLA(esnow); %the LTLA of each school
    Urbannow = UrbanRuralType2(esnow); %whether school is rural or urban
    closecontactsnewtermnow = combinedschools(esnow,8); %close contact sizes after new term

    %model Log Likelihood function
    LLfun2 = @(infparams,data_t,data_hist, data_t2, data_hist2)LLfunLTLAwithEaster(infparams, data_t, data_hist, Prob_profiles, Testing, Adherence, Strategy, closecontactsizesnow, yeargroupsizesnow, numyeargroupsnow, LTLAEasterHoliday, communityprevnow, numtestingposnow, numtestingnegnow, numpopLTLAnow, data_t2, data_hist2, schoolLTLAnow, Sigmoidfit, initsLTLA, Urbannow, closecontactsnewtermnow, Lmax);
        
     eps 
          
    while(any(LLs > eps)) %while any Log-likelihoods remain above threshold do run

        rng((j+1)*N^3 + timespergentemp*(N))

        sum(LLs > eps) %number of parameter sets above threshold


        paramchoices = randsample(1:N, N, true, w); %choose parameters to perturb

        for k = 1:N            
            newparams(k,:) =  mvnrnd(X(paramchoices(k), :), Kernel, 1); %geenrate new parameters
        end

        iterspergen(j+1) = iterspergen(j+1) + sum(isinf(LLs));  
        timespergen(j+1) = timespergen(j+1)+1;
        timespergentemp = timespergen(j+1); 

        if mod(timespergentemp, 20) == 0
            delete(gcp('nocreate')); %close parpool to reduce memory use
        end

        parfor k = 1:N

            if LLs(k) > eps

                Successes  = 0;  %number of runs below threshold - set to 0 initially             
                params(k,:) = newparams(k,:);               

                if (any(params(k,:) < paramsmin) || any(params(k,:) > paramsmax))
                    %out of prior range
                    Successes = 0;
                    LLs(k) = Inf;
                else

                    rng((j+1)*N^3 + timespergentemp*(N) + (k-1)); 

                    weighting(k) =  1/sum((w/sum(w)).*mvnpdf(params(k,:), X, Kernel)); %calculate weighting
                    tempLL = zeros(1,B); %log likelihoods for each run

                    for jj = 1:B
                        
                        tempLL(jj) = LLfun2(params(k,:), Pillar_data_weeks, Pillar_data_hist_2020, Pillar_data_5_weeks, Pillar_data_hist_2021);

                        if tempLL(jj) < eps
                            Successes = Successes+1;
                        end

                    end

                    LLs(k) = min(tempLL); %set Log likelihood as LL from best run
                    weighting(k) = Successes*weighting(k)/B; %set weighting according to number of successes

                end
            end

        end

    end

      
    save('Modelfittingworkspace');   
    delete(gcp('nocreate'));
  
 toc
  end
 
end

%Likelihood function for model fitting
 function [LL, LL_t, LL_t2, LL_hist1, LL_hist2]  = LLfunLTLAwithEaster(inputparams, data_t, data_hist, Prob_profiles, Testing, Adherence, Strategy, closecontactsizes, yeargroupsizes, numyeargroups, LTLAEasterHoliday, communityprevforschool, numtestingposforschool, numtestingnegforschool, numpopLTLA, data_t2, data_hist2, schoolLTLA, Sigmoidfit, initsLTLA, urbanorrural, closecontactsizesnewterm, Lmax)


 %% adjust parameters to be in previous order 
 infparams(1:7) = inputparams(1:7);
 infparams(8) = 0;
 infparams(9) = inputparams(8);
 infparams(10) = inputparams(9);
 infparams(11) = 0.9997;
 infparams(12) = inputparams(10);
    
    %% Infection params (that don't depend on randomness or LTLA)
    Infection.Weeks = length(data_t)+1; %Number of weeks (including week before term    
    Infection.alpha_betweenyears = 0.01;
    Infection.leak_infect = 0; %0 or 1, whether infection occurs on test days    
    Infection.Rec_0 = 0.0625;
    Infection.Sym_proportion = infparams(7);         
    Infection.K_asym = infparams(6);
    Infection.alpha_withinyear = 1;
   
    %% Testing_parameters 
    Testing.sens_PCR = 1; % parameter deciding which column of the PCR sensitivity profile to read
    %(1) - baseline sensitivity, (2) - low sensitivity, (3) - high sensitivity 
    Testing.sens_lat = 1; %parameter deciding which column of the LFT sensitivity profile to read.
    %(1) - baseline sensitivity, (2) - low sensitivity, (3) - high sensitivity 
    Testing.spec_PCR = 1; %specificity of PCR tests
    Testing.spec_lat = infparams(11); %specificity of LFT tests
    Testing.PCR_delay = 2; %delay on PCR tests
    
    %% Adherence parameters
    Adherence.C_isolate = 0; % constant scaling the external force of infection on isolating individuals
    Adherence.serialnonadherence = 0; %proportion of individuals who agree to be tested but never actually take home tests
    Adherence.takenproperly = 1;
    Adherence.isolatingproperly = zeros(1,500);

      
    %% Strategy params (set as isolation for baselines
    Strategy.isolation = 1; %(0) - not isolating, (1) - isolating close contacts, (2) - isolating year groups
    Strategy.SCT = 0; %(0) - no serial contact testing, (1) - serial contact testing close contacts, (2) - serial contact testing year groups %%overrides isolation
    Strategy.masstesting = 2; %(0) - no mass testing, (1) - weekly mass testing, (2) - twice weekly mass testing
    Strategy.initialtestdays = []; %initial testing days

    %new variant in each school's LTLA
    communitySgeneforschool = zeros(length(schoolLTLA), length(237:(238+ 7*Infection.Weeks)));
    for i = 1:length(schoolLTLA)
        communitySgeneforschool(i,:) = 1*(1-Sigmoidfit( 237:(238+ 7*Infection.Weeks), schoolLTLA(i))) + infparams(10)*Sigmoidfit(237:(238+7*Infection.Weeks), schoolLTLA(i));  
    end
    

    %% initialise output storage
    tests_t = zeros(length(schoolLTLA), 266);
    poslft_t = zeros(length(schoolLTLA), 266);
    size_school = zeros(length(schoolLTLA),1);
    Peak_covid_absences = zeros(length(schoolLTLA),1);
    Peak_covid_absences1 = Peak_covid_absences;
    Peak_covid_absences2 = Peak_covid_absences;


              
    %% do first 100
        for i = 1:100
            
           Infection.K = lognrnd((log(infparams(1)) - (infparams(9)^2/2)), infparams(9));
           Infection.K2 = infparams(3)*Infection.K;    
           Infection.K3 = Infection.K2;
            eee = i;

           %Set Easter Holiday according to LTLA
          if LTLAEasterHoliday(schoolLTLA(eee)) == 1
             Infection.HolidayWeek = [10, 18:28, 32, 33];
          else
             Infection.HolidayWeek = [10, 18:28, 33, 34]; 
          end

            Infection.Inf_0 = initsLTLA(schoolLTLA(eee));

          %adjust Ext according to school's urban or rural status
          if urbanorrural(eee) == 2
            Infection.Ext = lognrnd((log(infparams(2)) - (infparams(8)^2/2)), infparams(8));
          else
               Infection.Ext = infparams(5)*lognrnd((log(infparams(2)) - (infparams(8)^2/2)), infparams(8));      
          end

          Infection.HolidayExt = infparams(4);  

          Underreporting = ones(1,Infection.Weeks*7);
          Underreporting(30*7:end) = infparams(12);
          Adherence.probtakelatflow = ((Underreporting.*numtestingnegforschool(eee,:)) + numtestingposforschool(eee,:))/numpopLTLA(eee);                      
          Infection.commext = communityprevforschool(eee,:);            
          Prob_profiles.newvar = communitySgeneforschool(eee,:);           
         school_pop = SchoolPopulationquicker(yeargroupsizes(eee), numyeargroups(eee), closecontactsizes(eee));
         school_pop2 = SchoolPopulationquicker(yeargroupsizes(eee), numyeargroups(eee), closecontactsizesnewterm(eee));
         
         history = Interactingyeargroupsmuchquicker(school_pop, Infection, Testing, Strategy, Adherence, Prob_profiles, school_pop2); %do model run                 
         [pos_PCRsday, peakcovabs, Pos_LFTs] = Modeloutputssupercondensed(history, Infection.Weeks);  %obtain model outputs         
            
            %store outputs for each school
            tests_t(i,:) = pos_PCRsday;
            size_school(i) = yeargroupsizes(eee)*numyeargroups(eee);            
            Peak_covid_absences(i) = max(peakcovabs);
            poslft_t(i,:) = Pos_LFTs;
            Peak_covid_absences1(i) = max(peakcovabs(8:119));
            Peak_covid_absences2(i) = max(peakcovabs(end-76:end));


        end
       
        
        %model outputs for fitting
        
        tests_model_t = sum(tests_t)/sum(size_school);
        tests_model_hist = histcounts(Peak_covid_absences, [0:15 Inf]);
        tests_model_hist = tests_model_hist/sum(tests_model_hist);        
        tests_model_hist1 = histcounts(Peak_covid_absences1, [0:15 Inf]);
        tests_model_hist2 = histcounts(Peak_covid_absences2, [0:8 Inf]);
       
        
       %adjust tests_model_hist1 and 2 if there are any 0 values
       if any(~tests_model_hist1)
            tests_model_hist1 = tests_model_hist1+1;
       end
       
       if any(~tests_model_hist2)
            tests_model_hist2 = tests_model_hist2+1;
       end
       
       tests_model_hist1 = tests_model_hist1/sum(tests_model_hist1);
       tests_model_hist2 = tests_model_hist2/sum(tests_model_hist2);
       model_poslft = sum(poslft_t)/sum(size_school);        
 
         num_3 = 3896599; %for secondary schools

         tests_model_t2 = model_poslft;
         
         tests_model_t_weeks = zeros(1,length(Infection.Weeks)-1);
         tests_model_t_weeks2 = tests_model_t_weeks;

         %testing data in weeks rather than days
         for j = 1:(Infection.Weeks-1)
             tests_model_t_weeks(j) = sum(tests_model_t(1 +(j-1)*7: 7 + (j-1)*7));
             tests_model_t_weeks2(j) = sum(tests_model_t2(1 +(j-1)*7: 7 + (j-1)*7));

         end

 tests_model_t_weeks2 = tests_model_t_weeks2((end-10):end);
 
 %Log likelihood at optimum
 optllt = sum(data_t.*log(data_t/num_3) + (num_3 - data_t).*log(1 - (data_t/num_3)));
 optllt2 = sum(data_t2.*log(data_t2/num_3) + (num_3 - data_t2).*log(1 - (data_t2/num_3))); 
 optllhist = sum(data_hist.*log(data_hist/sum(data_hist)));
 optllhist2 = sum(data_hist2.*log(data_hist2/sum(data_hist2)));
 
 %Calculate Log likelihoods
 LL_t = sum(data_t.*log(tests_model_t_weeks) + (num_3 - data_t).*log(1 - tests_model_t_weeks)) - optllt;
 LL_t2 = (sum(data_t2.*log(tests_model_t_weeks2) + (num_3 - data_t2).*log(1 - tests_model_t_weeks2)) - optllt2);
 LL_hist1 = (sum(data_hist.*log(tests_model_hist1)) - (optllhist));% 
 LL_hist2 = (sum(data_hist2.*log(tests_model_hist2)) - (optllhist2));%

 
 %Calculate - total Log likelihood
 LL = -1*(LL_t + LL_hist1 + LL_t2 + LL_hist2);
 
 
 if length(schoolLTLA) > 100  && LL < Lmax
 %if under threshold do for the next 500 schools
     
        for i = 101:length(schoolLTLA)
            
               Infection.K = lognrnd((log(infparams(1)) - (infparams(9)^2/2)), infparams(9));
               Infection.K2 = infparams(3)*Infection.K;    
               Infection.K3 = Infection.K2;            
               eee = i;
           
          %Set Easter Holiday according to LTLA
          if LTLAEasterHoliday(schoolLTLA(eee)) == 1
             Infection.HolidayWeek = [10, 18:28, 32, 33];
          else
             Infection.HolidayWeek = [10, 18:28, 33, 34]; 
          end

          Infection.Inf_0 = initsLTLA(schoolLTLA(eee));

          %adjust Ext according to school's urban or rural status
          if urbanorrural(eee) == 2
            Infection.Ext = lognrnd((log(infparams(2)) - (infparams(8)^2/2)), infparams(8));
          else
               Infection.Ext = infparams(5)*lognrnd((log(infparams(2)) - (infparams(8)^2/2)), infparams(8));      
          end

          Infection.HolidayExt = infparams(4);  
          Underreporting = ones(1,Infection.Weeks*7);
          Underreporting(30*7:end) = infparams(12);
          Adherence.probtakelatflow = ((Underreporting.*numtestingnegforschool(eee,:)) + numtestingposforschool(eee,:))/numpopLTLA(eee);                
          Infection.commext = communityprevforschool(eee,:);           
          Prob_profiles.newvar = communitySgeneforschool(eee,:); 
         school_pop = SchoolPopulationquicker(yeargroupsizes(eee), numyeargroups(eee), closecontactsizes(eee));
         school_pop2 = SchoolPopulationquicker(yeargroupsizes(eee), numyeargroups(eee), closecontactsizesnewterm(eee));
         
         
         history = Interactingyeargroupsmuchquicker(school_pop, Infection, Testing, Strategy, Adherence, Prob_profiles, school_pop2); %do model run                  
         [pos_PCRsday, peakcovabs, Pos_LFTs] = Modeloutputssupercondensed(history, Infection.Weeks); %model outputs          

          %store outputs for each school
            tests_t(i,:) = pos_PCRsday;
            size_school(i) = yeargroupsizes(eee)*numyeargroups(eee);            
            Peak_covid_absences(i) = max(peakcovabs);
            poslft_t(i,:) = Pos_LFTs;         
            Peak_covid_absences1(i) = max(peakcovabs(8:119));
            Peak_covid_absences2(i) = max(peakcovabs(end-76:end));

        end
       
        
      %model outputs for fitting
        
       tests_model_t = sum(tests_t)/sum(size_school);      
       tests_model_hist1 = histcounts(Peak_covid_absences1, [0:15 Inf]);
       tests_model_hist2 = histcounts(Peak_covid_absences2, [0:8 Inf]);
       
       %adjust tests_model_hist1 and 2 if there are any 0 values
       if any(~tests_model_hist1)
            tests_model_hist1 = tests_model_hist1+1;
       end
       
       if any(~tests_model_hist2)
            tests_model_hist2 = tests_model_hist2+1;
       end
       
       tests_model_hist1 = tests_model_hist1/sum(tests_model_hist1);
       tests_model_hist2 = tests_model_hist2/sum(tests_model_hist2);
       model_poslft = sum(poslft_t)/sum(size_school);        

     num_3 = 3896599; %for secondary schools

     tests_model_t2 = model_poslft;
 
  %testing data in weeks rather than days
 for j = 1:(Infection.Weeks-1)
     tests_model_t_weeks(j) = sum(tests_model_t(1 +(j-1)*7: 7 + (j-1)*7));
     tests_model_t_weeks2(j) = sum(tests_model_t2(1 +(j-1)*7: 7 + (j-1)*7));

 end
 %}
 
 tests_model_t_weeks2 = tests_model_t_weeks2((end-10):end);
 
 
 

  %Log likelihood at optimum
 optllt = sum(data_t.*log(data_t/num_3) + (num_3 - data_t).*log(1 - (data_t/num_3)));
 optllt2 = sum(data_t2.*log(data_t2/num_3) + (num_3 - data_t2).*log(1 - (data_t2/num_3))); 
 optllhist = sum(data_hist.*log(data_hist/sum(data_hist)));
 optllhist2 = sum(data_hist2.*log(data_hist2/sum(data_hist2)));
 
  %Calculate Log likelihoods
 LL_t = sum(data_t.*log(tests_model_t_weeks) + (num_3 - data_t).*log(1 - tests_model_t_weeks)) - optllt;
 LL_t2 = (sum(data_t2.*log(tests_model_t_weeks2) + (num_3 - data_t2).*log(1 - tests_model_t_weeks2)) - optllt2);
 LL_hist1 = (sum(data_hist.*log(tests_model_hist1)) - (optllhist));% 
 LL_hist2 = (sum(data_hist2.*log(tests_model_hist2)) - (optllhist2));%

 
 %Calculate - total Log likelihood
 LL = -1*(LL_t + LL_hist1 + LL_t2 + LL_hist2);
     
     
 end

end

 
 


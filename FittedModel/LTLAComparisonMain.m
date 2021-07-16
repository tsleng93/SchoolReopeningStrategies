function LTLAComparisonMain


numtestingnegforschool = [];
numtestingposforschool = [];
communityprevforschool = [];
combinedschools = [];
schoolLTLA = [];
schoolregion = [];
load('exampleworkspace2.mat');




tests_vec1region = zeros(8, 100, 266);
poslft_vec1region = zeros(8, 100, 266);


        

        p_rewire = 0;
        


       %Probability Profiles
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


 poolobj = gcp;    
 
 
 paramsin = params;
 schoolLTLAin = schoolLTLA;
 Sigmoidfitin = Sigmoidfit;
 LTLAEasterHolidayin = LTLAEasterHoliday;
 initsLTLAin = initsLTLA;
 urbanorruralin = urbanorrural;
 numpopLTLAin = numpopLTLA;
 
    tic
  parfor j = 1:100

            j

      infparams = paramsin(j,1:12);
      
      communitySgeneforschool = zeros(2979, 274);
    
   for i = 1:length(schoolLTLAin)
    communitySgeneforschool(i,:) = 1*(1-Sigmoidfitin( 237:(489+21), schoolLTLAin(i))) + infparams(10)*Sigmoidfitin(237:(489+21), schoolLTLAin(i));  

    

    end
      
      
      
      
      
     LTLAtests = zeros(1,317);
     LTLAsize = zeros(1,317);
     LTLALFTs = zeros(1,317);
    
     
     Infection = [];
     Strategy = [];
     Adherence = [];
     Prob_profiles = [];
     Testing = [];
     

    
    Prob_profiles.PCRsym = PCR_test_sym;
    Prob_profiles.PCRasym = PCR_test_asym;
    Prob_profiles.latsym = lat_test_sym;
    Prob_profiles.latasym = lat_test_asym;
    Prob_profiles.Infectivity = Infectivity_since_infection;
    Prob_profiles.Symptom_onset = Symptom_onset;

     
     
     size_school = zeros(1,2979);
     tests_t = zeros(2979, 266);
     LFTs_t =  zeros(2979, 266);
     
     Absences_t = zeros(2979, 266);
     Peak_covid_absences1 = zeros(1,2979);
     Peak_covid_absences2 = zeros(1,2979);
     
     
     Prevs_t = zeros(2979, 266);
     Isolated_Infs = zeros(2979, 266);
     
     Rinfs1 = zeros(2979,266);
     Rinfs2 = zeros(2979,266);
     

    ExternalIncs = zeros(2979,266);
    InternalIncs = zeros(2979,266);
     
    
    
    %Set strategy
    Strategy.masstesting = 2; %(0) - no mass testing, (1) - weekly mass testing, (2) - twice weekly mass testing
    Strategy.isolation = 1;
    Strategy.SCT = 0; %(0) - no serial contact testing, (1) - serial contact testing close contacts, (2) - serial contact testing year groups 
    Strategy.SCTuptake = 1; %Proportion of individuals agreed to participate in serial contact testing
    Strategy.masstestinguptake = 1; %Proportion of individuals agreed to be mass tested
    Strategy.initialtestdays = []; %initial testing days
    Strategy.initialuptake = 1; %uptake to initial testing
     
    
    %Set infection
    Infection.Weeks = 39; %Number of weeks (including week before term)       
   Infection.alpha_betweenyears = 0.01; %Interaction with other year groups
   Infection.alpha_withinyear = 0.1; %Interaction with non-close contact members of year group
   Infection.leak_infect = 0;
   Infection.Rec_0 = 0.0625;
   Infection.HolidayExt = infparams(4);
   Infection.K_asym = infparams(6);
   Infection.Sym_proportion = infparams(7);


  %Set adherence
    Adherence.isolatingproperly = zeros(1,500);    
    Adherence.takenproperly = 1;
    Adherence.C_isolate = 0; % constant scaling the external force of infection on isolating individuals
    Adherence.serialnonadherence = 0; %proportion of individuals who agree to be tested but never actually take home tests
    
    
    %Set testing
    Testing.spec_lat = infparams(11); %specificity of LFT tests
    %Testing_parameters 
    Testing.sens_PCR = 1; % parameter deciding which column of the PCR sensitivity profile to read
    %(1) - baseline sensitivity, (2) - low sensitivity, (3) - high sensitivity 
    Testing.sens_lat = 1; %parameter deciding which column of the LFT sensitivity profile to read.
    %(1) - baseline sensitivity, (2) - low sensitivity, (3) - high sensitivity 
    Testing.spec_PCR = 1; %specificity of PCR tests
   %Testing.spec_lat = infparams(8);
    Testing.PCR_delay = 2; %delay on PCR tests
    

       %loop through schools       
        for i = 1:2979
     
        
       Infection.K = lognrnd(log(infparams(1)) - (infparams(9)^2/2), infparams(9));
      

    
    eee = i;
   

      %Set Easter Holidays
      if LTLAEasterHolidayin(schoolLTLAin(eee)) == 1
         Infection.HolidayWeek = [10, 18:28, 32, 33];
      else
         Infection.HolidayWeek = [10, 18:28, 33, 34]; 
      end
  
  
    
   Infection.Inf_0 = initsLTLAin(schoolLTLAin(eee));    
   Infection.K2 = infparams(3)*Infection.K;   
   Infection.K3 = Infection.K2;
  
  if urbanorruralin(eee) == 2
    Infection.Ext = lognrnd((log(infparams(2)) - (infparams(8)^2/2)), infparams(8));

  else
       Infection.Ext = infparams(5)*lognrnd((log(infparams(2)) - (infparams(8)^2/2)), infparams(8));      
  end
      
  
  %Using positive LFT data
  %{
  Underreporting = ones(1,Infection.Weeks*7);
  Underreporting(30*7:end) = infparams(12);
  Adherence.probtakelatflow = ((Underreporting.*numtestingnegforschool(i,:)) + numtestingposforschool(i,:))/numpopLTLAin(i);
  %}
  
  
  %set uptake to be 36%
  Adherence.probtakelatflow = zeros(1, Infection.Weeks*7);
  Adherence.probtakelatflow(28*7:end) = 0.36*(2/7);
  
           
           Infection.commext = communityprevforschool(eee,:);
           
          Prob_profiles.newvar = communitySgeneforschool(eee,:);
         % Prob_profiles.newvar = ones(1,(Infection.Weeks*7) +1); %If no
         % impact of new variant
           
         

            
            [school_pop, ~] = SchoolPopulation(combinedschools(eee,3), combinedschools(eee,2), combinedschools(eee,4), p_rewire);
            [school_pop2, ~] = SchoolPopulation(combinedschools(eee,3), combinedschools(eee,2), combinedschools(eee,8), p_rewire);  


         
            
            history = Interactingyeargroupsextended(school_pop, Infection, Testing, Strategy, Adherence, Prob_profiles, randi(1e6), school_pop2);
            %history = Interactingyeargroupsquicker(school_pop, Infection, Testing, Strategy, Adherence, Prob_profiles, randi(1e6), school_pop2);

           [pos_PCRsday, peakcovabs, Absences, Prev_true, Num_LFTs, Pos_LFTs, Pos_LFTsandPCRs, Inc] = Modeloutputscondensed(history, Infection.Weeks);           

            
           
      
           
            %record information at LTLA level
            LTLAtests(schoolLTLA(eee)) = LTLAtests(schoolLTLA(eee)) + sum(pos_PCRsday);
            LTLALFTs(schoolLTLA(eee)) = LTLALFTs(schoolLTLA(eee)) + sum(Pos_LFTs);
            LTLAsize(schoolLTLA(eee)) = LTLAsize(schoolLTLA(eee)) + combinedschools(eee,3)*combinedschools(eee,2);

           %record information at school level
            tests_t(i,:) = pos_PCRsday;
            size_school(i) = combinedschools(eee,3)*combinedschools(eee,2);  
            Peak_covid_absences1(i) = max(peakcovabs(8:119));
            Peak_covid_absences2(i) = max(peakcovabs(end-76:end));
           LFTs_t(i,:) = Pos_LFTs;           
            Absences_t(i,:) = sum(Absences(:, 8:end));
           Prevs_t(i,:) = sum(Prev_true(:, 8:end));
           Isolated_Infs(i,:) = sum(history.Isolated_Infecteds(:, 8:end));           
            Rinfs1(i,:) = history.Rday(1,8:end-10);
            Rinfs2(i,:) = history.Rday(2,8:end-10);           
            ExternalIncs(i,:) = squeeze(sum(sum((history.ext_or_int(:,:,8:end-10) == 1))));
            InternalIncs(i,:) = squeeze(sum(sum((history.ext_or_int(:,:,8:end-10) == 2))));
         

        end
        
        %record aggregated information for each parameter set
        
        tests_vec1(j,:) = sum(tests_t)/sum(size_school(1:2979));
        hist_vec1(j,:) = histcounts(Peak_covid_absences1, [0:50 Inf]);
        hist_vec2(j,:) = histcounts(Peak_covid_absences2, [0:50 Inf]);
        poslft_vec1(j,:) = sum(LFTs_t)/sum(size_school(1:2979));        
        abs_vec1(j,:) = sum(Absences_t)/sum(size_school(1:2979));        
        isolinfprop_vec1(j,:) = sum(Isolated_Infs)./sum(Prevs_t);
                
        for k = 2:8
            tests_vec1region(k,j,:) = sum(tests_t(schoolregion == k,:))/sum(size_school(schoolregion == k));
            poslft_vec1region(k,j,:) = sum(LFTs_t(schoolregion ==k, :))/sum(size_school(schoolregion == k));         
        end
        
        Rinfs1_vec1(j,:) = sum(Rinfs1);
        Rinfs2_vec1(j,:) = sum(Rinfs2);
        
        ExternalInc_vec1(j,:) = sum(ExternalIncs)/sum(size_school(1:2979));
        InternalInc_vec1(j,:) = sum(InternalIncs)/sum(size_school(1:2979));
        %}
      
        Prevs_vec1(j,:) = sum(Prevs_t)/sum(size_school(1:2979));        
        LTLAtestsmatrix(j,:) = LTLAtests;
        LTLALFTsmatrix(j,:) = LTLALFTs;
        LTLAsizematrix(j,:) = LTLAsize;
        
        
        

              
        
  end

  
  toc
  
  save('Exampleoutput');
  
end

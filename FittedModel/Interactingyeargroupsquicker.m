function [history] = Interactingyeargroupsquicker(school_pop, Infection, Testing, Strategy, Adherence, Prob_profiles, randnum, school_pop2)
%This function simulates the spread of infection within a secondary school
%formed of year groups.

%Input:
%   - params is a parameter vector (meaning of each entry is described
%   below)
%   - PCR_test_sym is the PCR test positivity profile for symptomatic individuals
%   - PCR_test_asym is the PCR test positivity profile for asymptomatic individuals
%   - lat_test_sym is the LFT positivity profile for symptomatic individuals
%   - lat_test_asym is the LFT positivity profile for asymptomatic individuals
%   - randnum is a positive integer to seed the random numbers chosen for reproducibility

%Output:
%   - history is a structure, containing the infection, testing, and
%   isolation details throughout the course of the simulation. Useful
%   outputs from this can be obtained using the files Modeloutputs.m and
%   Moremodeloutputs.m

%Authors: Trystan Leng and Edward M. Hill 
%Last update 15/07/2021.

if nargin == 0
    %random number   
    randnum = 5;
    rng(randnum);
    
    %School population
    school_pop = SchoolPopulation(200, 5, 200, 0);
    
    
    %Infection_parameters
    Infection.K = 1; %secondary infections expected from symptomatic individual (without depletion of susceptibles)
    Infection.K2 = 1;%after October break
    Infection.K3 = 1;%after Easter break
       
    Infection.Ext = 1; %community infection scaling factor
    Infection.Inf_0 = 0.02; %Initial population prevalence
    Infection.Rec_0 = 0.2;  %Initial population immunity
    Infection.K_asym = 0.4*rand + 0.3; % Relative infectiousness of asymptomatic individuals
    Infection.Sym_proportion = 0.12 + 0.19*rand; %Proportion of school population who are symptomatic
    Infection.alpha_withinyear = 0; %Interaction with other members of the year group (not close contacts)
    Infection.alpha_betweenyears = 0; %Interaction between year groups
    Infection.leak_infect = 0; %0 or 1, whether infection occurs on test days
    
    
    Infection.Weeks = 39; %Number of weeks (including week before term)
    
    Infection.HolidayWeek = [10, 18:28, 33, 34];% Weeks of school holidays
    
    %Testing_parameters
    Testing.sens_PCR = 1; % parameter deciding which column of the PCR sensitivity profile to read
    %(1) - baseline sensitivity, (2) - low sensitivity, (3) - high sensitivity 
    Testing.sens_lat = 1; %parameter deciding which column of the LFT sensitivity profile to read.
    %(1) - baseline sensitivity, (2) - low sensitivity, (3) - high sensitivity 
    Testing.spec_PCR = 1; %specificity of PCR tests
    Testing.spec_lat = 0.9999; %specificity of LFT tests
    Testing.PCR_delay = 2; %delay on PCR tests
    
    %Strategy parameters
    Strategy.isolation = 1; %(0) - not isolating, (1) - isolating close contacts, (2) - isolating year groups
    Strategy.SCT = 0; %(0) - no serial contact testing, (1) - serial contact testing close contacts, (2) - serial contact testing year groups %%overrides isolation
    Strategy.SCTuptake = 1; %Proportion of individuals agreed to participate in serial contact testing
    Strategy.masstesting = 2; %(0) - no mass testing, (1) - weekly mass testing, (2) - twice weekly mass testing
    Strategy.masstestinguptake = 1; %Proportion of individuals agreed to be mass tested
    Strategy.initialtestdays = [1 4 8]; %initial testing days
    Strategy.initialuptake = 0; %uptake to initial testing
    
    %Adherence parameters
    Adherence.C_isolate = 0; % constant scaling the external force of infection on isolating individuals
    Adherence.probtakelatflow = 0.05*ones(1, Infection.Weeks*7); %probability of taking lat flow (0 - 1), larger integers could specify the inclusion of correlations into this   
    Adherence.serialnonadherence = 0; %proportion of individuals who agree to be tested but never actually take home tests    
    Adherence.takenproperly = 1; %proportion of LFTs that are not taken properly
        
    Adherence.isolatingproperly = zeros(1, Infection.Weeks*7); %probability of schools *NOT* isolating students properly upon confirmation of a positive case, always been set to 0
    
    %Probability Profiles
    PCR_test_sym = readtable('PCR_Curve_summary.csv');
    PCR_test_sym = table2array(PCR_test_sym(:, 2:4));
    PCR_test_asym = csvread('PCR_Curve_asym.csv');
    lat_test_sym = readtable('lat_Curve_summary.csv');
    lat_test_sym = table2array(lat_test_sym(:, 2:4));
    lat_test_asym = csvread('lat_Curve_asym.csv');
    
    Prob_profiles.PCRsym = PCR_test_sym;
    Prob_profiles.PCRasym = PCR_test_asym;
    Prob_profiles.latsym = lat_test_sym;
    Prob_profiles.latasym = lat_test_asym;
        
    Infectivity_since_infection = [0.0063 0.0563 0.1320 0.1798 0.1817 0.1521 0.1117 0.0746 0.0464 0.0272 0.0152 0.0082 0.0043 0.0022 0.0011 0.0005 0.0002 0.0001 0.0001 0.0000];
    Infectivity_since_infection = Infectivity_since_infection'/sum(Infectivity_since_infection);
    Infectivity_since_infection(15) = sum(Infectivity_since_infection(15:end)); Infectivity_since_infection = Infectivity_since_infection(1:15);
    Prob_profiles.Infectivity = Infectivity_since_infection;
    
    Symptom_onset =  [0.0055 0.0534 0.1307 0.1814 0.1847 0.1545 0.1129 0.0747 0.0458 0.0265 0.0146 0.0077 0.0039 0.0020 0.0010 0.0005 0.0002 0.0001 0 0];
    Symptom_onset = Symptom_onset/sum(Symptom_onset);
    Symptom_onset(15) = sum(Symptom_onset(15:end)); Symptom_onset = Symptom_onset(1:15);
    Prob_profiles.Symptom_onset = Symptom_onset;
    
    
    %time-varying community infection
    Infection.commext = ones(1,Infection.Weeks*7);
    
    %effect of new variant
    Prob_profiles.newvar = ones(1,Infection.Weeks*7 +1);
    
end


%%Parameters%%

%Directly from input

%get infection parameters used often
Ext = Infection.Ext; 
alpha_withinyear = Infection.alpha_withinyear; 
alpha_betweenyears = Infection.alpha_betweenyears;
leak_infect = Infection.leak_infect;
Weeks = Infection.Weeks;

%get testing parameters used often
PCR_delay = Testing.PCR_delay;

%%Inferred from other parameters%%

%1 = baseline, 2 = lower bound, 3 = upper bound
%Probability of PCR tests
Prob_profiles.PCR_symday = [0;  Prob_profiles.PCRsym(10*(1:30) + 1, Testing.sens_PCR); zeros(max(51,Weeks*7-31) + 50,1)];
Prob_profiles.PCR_asymday = [0; Prob_profiles.PCRasym(10*(1:30) + 1, Testing.sens_PCR); zeros(max(51,Weeks*7-31) + 50,1)];

%Probability of lateral flow tests
Prob_profiles.lat_symday = [0; Prob_profiles.latsym(10*(1:30) + 1, Testing.sens_lat); zeros(max(51,Weeks*7-31) + 50,1)];
Prob_profiles.lat_asymday = [0; Prob_profiles.latasym(10*(1:30) + 1, Testing.sens_lat); zeros(max(51,Weeks*7-31) + 50,1)];

%Adjust for specificity
Prob_profiles.lat_symday(Prob_profiles.lat_symday < (1- Testing.spec_lat)) = (1-Testing.spec_lat);
Prob_profiles.lat_asymday(Prob_profiles.lat_asymday < (1-Testing.spec_lat)) = (1-Testing.spec_lat);


Prob_profiles.lat_symday = Adherence.takenproperly*Prob_profiles.lat_symday;
Prob_profiles.lat_asymday = Adherence.takenproperly*Prob_profiles.lat_asymday;

%also set PCR specificity
Prob_profiles.PCR_symday(Prob_profiles.PCR_symday < (1- Testing.spec_PCR)) = (1-Testing.spec_PCR);
Prob_profiles.PCR_asymday(Prob_profiles.PCR_asymday < (1-Testing.spec_PCR)) = (1-Testing.spec_PCR);


%set year sizes and year group sizes
YearSize =  school_pop.yearsize;
YearGroup = school_pop.yeargroups;


%%Infectivity profile%%
Prob_profiles.Inf_symday = [0;Prob_profiles.Infectivity; zeros(max(51,Weeks*7-15 + 50),1)];
Prob_profiles.Inf_asymday = Infection.K_asym*Prob_profiles.Inf_symday;

%Keep entire history
history.Infection = zeros(YearGroup, YearSize, Weeks*7+10); %Infection status of individuals through time
history.Isolation = zeros(YearGroup, YearSize, Weeks*7+10); %Whether individuals are isolated or not through time


history.Rday = zeros(2, Weeks*7+10); %for calculating R
history.DayInf = zeros(YearGroup,YearSize); %day someone gets infected


history.TotInfstudent = zeros(YearGroup, YearSize, Weeks*7+10); %total infection to a student

history.ext_or_int = zeros(YearGroup, YearSize,Weeks*7+10); %Infected externally or from schoo

history.Isolated_Infecteds = zeros(YearGroup, Weeks*7 + 10); %Number of isolated infecteds each day

history.pos_test_day = zeros(YearGroup, Weeks*7+10); %Number of positive tests that day

history.taken_test = zeros(YearGroup, YearSize, Weeks*7+10);
history.posLFT = zeros(YearGroup, YearSize, Weeks*7+10);
history.posLFTandPCR = zeros(YearGroup, YearSize, Weeks*7+10);


%{
history.ever_test_pos = zeros(YearGroup, YearSize); %2 if positive to PCR, 1 if positive to Lateral flow, 0 if never tested positive
history.presyms = zeros(YearGroup, YearSize); %Number Presymptomatic individuals captured
history.PCRteststaken = zeros(YearGroup, Weeks*7+10); %count PCRs taken
history.LFTteststaken = zeros(YearGroup, Weeks*7+10); %count LFTs taken


%External infection
history.ExtInfstudent = zeros(YearGroup, YearSize, Weeks*7 + 10);
%}


%{
history.ext_or_int = [];
history.Isolated_Infecteds = [];
history.Known_Infecteds = [];
history.Known_asymptomatics = [];
history.ever_test_pos = [];
history.presyms = [];
history.PCRteststaken = [];
history.LFTteststaken = [];
history.ExtInfstudent = [];
%}

history.IsolatingthroughCovid = zeros(YearGroup, YearSize,Weeks*7 +10);
%history.IsolatingthroughCovid = [];


%Testing days
history.gettingatest = zeros(YearGroup, YearSize, Weeks*7+10); %2 for in-school test, 1 for home test


%% Things that depend on randomness
rng(randnum);
%generate random numbers beforehand
%for PCR tests, for LF tests, and for infection       
r_lat = rand(YearGroup, YearSize, Weeks*7);
r_inf = rand(YearGroup, YearSize, Weeks*7);
r_PCR = rand(YearGroup, YearSize, Weeks*7);

%for confirmatory PCR tests
r_conf = rand(YearGroup, YearSize, Weeks*7);


%random for adherence

r_ad1 = rand(YearGroup,YearSize, Weeks*7);
r_ad2 = rand(YearGroup, YearSize, Weeks*7);


Ks = Infection.K*ones(YearGroup,YearSize);
Ks2 = Infection.K2*ones(YearGroup, YearSize);
Ks3 = Infection.K3*ones(YearGroup, YearSize);

%for testing day
r_testingday = randi(5, YearGroup, YearSize) - 1;
r_taketest = rand(YearGroup, YearSize, Weeks*7);

%Set who is symptomatic
Symptomatic = zeros(YearGroup, YearSize);
%who is symptomatic?
for yr = 1:YearGroup
    Symptomatic(yr, randsample(YearSize, binornd(YearSize, Infection.Sym_proportion))) = 1;
end

%day an individual will seek PCR test if symptomatic i.e. during of
%presymptomatic period    
day_of_PCR = zeros(YearGroup, YearSize);
PCR_days = randsample(1:15, sum(sum(Symptomatic)), true, Prob_profiles.Symptom_onset); 
day_of_PCR(Symptomatic == 1) = PCR_days;

will_take_PCR = zeros(YearGroup, YearSize);

%who is initially infected and recovered?
for yr = 1:YearGroup
    init = mnrnd(YearSize, [Infection.Inf_0; Infection.Rec_0; 1 - Infection.Inf_0 - Infection.Rec_0]);
    init_not_sus = randsample(YearSize, init(1)+init(2));
    init_Inf = init_not_sus(1:init(1));   
    init_Rec = init_not_sus((init(1)+1):(init(1)+init(2)));
               
    %randomly assign which day of infection infecteds are at
    history.Infection(yr, init_Inf,1) = randi(15, 1, init(1));
    %assume recovereds do not test positive to any tests
    history.Infection(yr, init_Rec,1) = 32;
   
   if ~isempty(init_Inf)
       rinit = rand(length(init_Inf),1);



        % testing for initially infecteds
        for student = init_Inf'
            
            
            history.DayInf(yr,student) = 1;

            if Symptomatic(yr, student)
                 will_take_PCR(yr,student) = day_of_PCR(yr, student) - history.Infection(yr, student,1 ) +1; 





                 if will_take_PCR(yr, student) < 1 
                     if rinit(end) <  Prob_profiles.PCR_asymday(day_of_PCR(yr, student))
                         %tests positive before simulation begins
                         history.Isolation(yr, student, 1: (1 + 8 + will_take_PCR(student))) = 1; %if tests is on day 0, then student should isolate  from day 1 to day 9                                                                                            %if test is on day -1, then student should isolate from day 1 to day 8 
                     end
                 end

            end  
                rinit(end) = [];
        end
   end
end

%set populations who agree to being tested
testing_pop = zeros(YearGroup, YearSize, 4); %1 for initial tests, 2 for mass tests, 3 for SCTs, 4 for adhering to mass tests
testing_props = [ Strategy.initialuptake,  1; Strategy.masstestinguptake, 2; Strategy.SCTuptake, 3];
testing_orders = sortrows(testing_props);





%set days of week pupils seek home LFTs
rdays = zeros(5, Weeks*7+10);
for rday = 0:4
for day = 1:Weeks*7
    if mod(day,7) == rday || mod(day,7) == mod(rday+3, 7) 
               rdays(rday+1,day) = 1;
    end
end
end





%Assume that the min of the 3 agree to all types of testing (serial contact
%testing, mass testing, initial testing))
for yr = 1:YearGroup   
    all_students = 1:YearSize;
    core_group = randsample(all_students, round(testing_orders(1,1)*YearSize));    
    %set core group who get tested to all
    testing_pop(yr, core_group, testing_orders(1,2)) = 1;    
    testing_pop(yr,:, testing_orders(2,2)) = testing_pop(yr,:, testing_orders(1,2));
    
    remaining_students = all_students; 
    remaining_students(core_group) = [];

    next_group = randsample(remaining_students, round((testing_orders(2,1) - testing_orders(1,1))*YearSize));
    
    %set next group who get tested to two of the tests
    testing_pop(yr, next_group, testing_orders(2,2)) = 1; 
    testing_pop(yr,:, testing_orders(3,2)) = testing_pop(yr,:, testing_orders(2,2));

    remaining_students = all_students; 
    remaining_students([core_group, next_group]) = [];
    
    final_group = randsample(remaining_students, round((testing_orders(3,1) - testing_orders(2,1))*YearSize));
    
    %set final group to be those who agree to all three tests
    testing_pop(yr, final_group, testing_orders(3,2)) = 1;
            
    %non-adhering students
    testing_pop(yr, :, 4) = testing_pop(yr,:,2);
    
    remaining_students = all_students;
    remaining_students(testing_pop(yr,:,2) == 0) = [];
    
    %remove non-ahdering students
    nonadhering_group = randsample(remaining_students, round(Adherence.serialnonadherence*Strategy.masstestinguptake*YearSize));
    testing_pop(yr, nonadhering_group, 4) = 0;
    
end


%set test days for mass tests for adhering students
if Strategy.masstesting == 1
  %Test Monday to Friday
  for day = 1:Weeks*7
    
     
    for yr = 1:YearGroup
        dow = mod(yr, 5);
        if dow == 0
            dow = 5;
        end

        if mod(day, 7) == dow
            history.gettingatest(yr, testing_pop(yr,:,4) == 1, day) = 1;
        end
    end

   end
    
elseif Strategy.masstesting == 2
  
    %tests occur on one of the following:
    %Mon-Thur, Tue-Fri, Wed-Sat, Thurs-Sun, Fri-Mon
    
    for yr = 1:YearGroup
        for student = find(testing_pop(yr,:,4) == 1)
            r = r_testingday(yr, student);
         
            history.gettingatest(yr,student,:) = rdays(r+1,:);
            
            
        end
    end
   %}
end

%find proportion of students being tested each day
prop_tested = squeeze(sum(sum(history.gettingatest > 0)))/(YearSize*YearGroup);

%set test days for initial opt-ins (overrides mass test days so has to be
%set after)
for yr = 1:YearGroup
   history.gettingatest(yr, testing_pop(yr,:,1) == 1, Strategy.initialtestdays) = 2;   
end





%Individual's infection probability and testing probabilities
Now.Infect = zeros(YearGroup, YearSize); 
Now.PCR = zeros(YearGroup, YearSize);
Now.lat = zeros(YearGroup, YearSize);
Now.had_PCR = zeros(YearGroup, YearSize); %have they had a +ve PCR?

%Initial test and infection probabilities 
Now = Update2(Now, history, Prob_profiles, Symptomatic, 1, Ks);


%define all_students vec
all_students = 1:YearSize;





YearGroupMatrix1 = full(subMatrix(school_pop.school_matrix, yr, YearSize));
YearGroupMatrix2 = full(subMatrix(school_pop2.school_matrix, yr, YearSize));

%For quicker
%YearGroupMatrix1 = (subMatrix(school_pop.school_matrix, yr, YearSize));
%YearGroupMatrix2 = (subMatrix(school_pop2.school_matrix, yr, YearSize));



Stratisolinit = Strategy.isolation;
Stratmasstestinit = Strategy.masstesting;

%Simulation
for Week = 1:Weeks
    
    
   %Strategy of isolating close contacts before then
    if Week < 28
       YearGroupMatrix = YearGroupMatrix1;
       NumCloseContacts = sum(YearGroupMatrix(1,:));
       
       
        Strategy.isolation = 1;
        Strategy.masstesting = 0;
       
       
    else
    %Inputted strategy afterwards,    
        YearGroupMatrix = YearGroupMatrix2;
        NumCloseContacts = sum(YearGroupMatrix(1,:));
        
        Strategy.isolation = Stratisolinit;
        Strategy.masstesting = Stratmasstestinit;

    end
%}
    
    for day_of_week = 1:7
        
        
        

 

        
        
        day = (Week-1)*7 + day_of_week;
        
        
       if any(Week == Infection.HolidayWeek)
           
          %don't count this during holidays
          history.IsolatingthroughCovid(:,:,day) = 0;

       end
        
        


       %External force of infection depending on community prevalence    
       Ext = Infection.Ext*Infection.commext(day);
       Ext_isolate = Adherence.C_isolate*Ext;
       
       
       
             
        %Test students via PCR (self seeking) with probability       
        for yr = 1:YearGroup
            

           
            
            for student = all_students(will_take_PCR(yr,:) == day)
                    if Now.had_PCR(yr, student) == 0 
                        %student seeks test
                        if r_PCR(yr, student, day) < Now.PCR(yr, student)
                            %student tests positive
                            Isolation_period = 10;
                            
                            notisolmarker = 0;
                            Now.had_PCR(yr, student) = 1;

                            
                            
                            if history.Isolation(yr,student,day) == 0
                                notisolmarker = 1;
                            end
                            
                            history.Isolation(yr, student, day:(day+Isolation_period)) = 1; %isolates for 10 days
                            history.IsolatingthroughCovid(yr,student,(day+PCR_delay):(day + Isolation_period)) = 1;

                            
                           history.pos_test_day(yr, day) = history.pos_test_day(yr, day) +1;

                            if (Week > 1 && ~any(Week == Infection.HolidayWeek)) && notisolmarker
                                    %school only isolates and tests  after children are back in
                                    %school

                                if Strategy.isolation == 1
                                    %only isolating close contacts
                                    isolgroup = all_students(YearGroupMatrix(student,:) == 1);  
                                    %isolation strategy not compatible with SCT
                                    sctgroup = [];
                                elseif Strategy.isolation == 2
                                    %isolating whole year group
                                    isolgroup = all_students;
                                    %isolation strategy not compatible with SCT
                                    sctgroup = [];
                                elseif Strategy.SCT == 1
                                    %close contacts not agreed to SCT
                                    isolgroup = all_students(YearGroupMatrix(student,:) == 1 & testing_pop(yr,:,3) == 0);
                                    %close contacts for sct
                                    sctgroup = all_students(YearGroupMatrix(student,:) == 1 & testing_pop(yr,:,3) == 1);
                                elseif Strategy.SCT == 2
                                    %whole year group not agreed to SCT
                                    isolgroup = all_students( testing_pop(yr,:,3) == 0);
                                    %whole year group for sct
                                    sctgroup = all_students(testing_pop(yr,:,3) == 1);
                                else
                                    %no isolation or sct
                                    isolgroup = [];
                                    sctgroup = [];
                                end
                                
                                
                                
                              if sum(sum(history.IsolatingthroughCovid(:,:,day-1))) > 0 || r_ad1(yr,student,day) > Adherence.isolatingproperly(day)
                                %isolate other individuals in group (who were
                                    %actually in school the day before)
                                    for otherpupil = isolgroup
                                           history.Isolation(yr, otherpupil, (day+PCR_delay+1):(day+Isolation_period)) = 1; %isolates indivs not getting tested
                                                                              % end
                                    end
                                
                               end

                                %change test days for those in SCT group
                                for otherpupil = sctgroup
                                    history.gettingatest(yr,otherpupil, (day+PCR_delay + 1): day+7) = 2;
                                end                         
                            end                            
                                                                                 
                        else
                            %student tests negative
                            history.Isolation(yr, student, day:(day+PCR_delay)) = 1; %isolates until PCR returns
                        end
                    end
            end
        end
         
       
       
       
        

         
        if Week > 27
        
         %mass testing   
            
              for yr = 1:YearGroup
                 
                 %find mass testing group
                 masstestinggroup =  all_students(history.gettingatest(yr,:,day) > 0 & will_take_PCR(yr,:) ~= day);
                 
                 
                 
                 for student = masstestinggroup
                     
                     
                   %To match Testing uptake
                   if history.gettingatest(yr,student,day) == 1
                         P_taketest = min(Adherence.probtakelatflow(day)/prop_tested(day), 1); %tested at home
                      else
                         P_taketest = 1; %tested in school
                   end
               %}
                   
                   
                   %For % uptake plots
                   %P_taketest = Adherence.probtakelatflow(day);

                      
                      if r_taketest(yr,student,day) < P_taketest
                         history.taken_test(yr,student,day) = 1; 
                    
                            if r_lat(yr, student, day) < (1 - history.Isolation(yr,student,day))*Now.lat(yr,student)
                             
                             %tests LFT positive   
                                
                             history.posLFT(yr,student,day) = 1;   


                               %confirmatory PCR
                                  %lateral flow tests positive
                                  history.ever_test_pos(yr,student) = 1;

                                 if history.Infection(yr, student, day) == 0 || history.Infection(yr, student, day) > 15
                                      %false positive, but recovered, so may
                                      %test positive by chance                                 
                                     if  r_conf(yr, student, day) < Now.PCR(yr, student)
                                          %tests positive
                                          Isolation_period = 10;
                                          Testing_period = 7;                                 
                                          %Assume they will not seek another PCR test
                                          Now.had_PCR(yr,student) = 1; 
                                           
                                     history.posLFTandPCR(yr,student,day) = 1;
                                     history.IsolatingthroughCovid(yr,student,(day+PCR_delay):(day + Isolation_period)) = 1;

                                     else
                                         %tests negative 
                                         Isolation_period = PCR_delay;
                                         Testing_period = PCR_delay;
                                     end

                                  else
                                      %true positive so confirmatory PCR comes
                                      %back positive
                                      Isolation_period = 10;
                                      Testing_period = 7;                                 
                                      %Assume they will not seek another PCR test
                                      Now.had_PCR(yr,student) = 1; 

                                       history.posLFTandPCR(yr,student,day) = 1;
                                       history.IsolatingthroughCovid(yr,student,(day+PCR_delay):(day + Isolation_period)) = 1;

                                  end


                                  history.Isolation(yr, student, (day):(day+Isolation_period)) = 1; %isolates +ve individual

                                  %set isolation groups and SCT groups
                                  %according to strategy
                                  if day > 8  && ~any(Week == Infection.HolidayWeek)
                                     if Strategy.isolation == 1
                                        %only isolating close contacts
                                        isolgroup = all_students(YearGroupMatrix(student,:) == 1);  
                                        %isolation strategy not compatible with SCT
                                        sctgroup = [];
                                    elseif Strategy.isolation == 2
                                        %isolating whole year group
                                        isolgroup = all_students;
                                        %isolation strategy not compatible with SCT
                                        sctgroup = [];
                                    elseif Strategy.SCT == 1
                                        %close contacts not agreed to SCT
                                        isolgroup = all_students(YearGroupMatrix(student,:) == 1 & testing_pop(yr,:,3) == 0);
                                        %close contacts for sct
                                        sctgroup = all_students(YearGroupMatrix(student,:) == 1 & testing_pop(yr,:,3) == 1);
                                    elseif Strategy.SCT == 2
                                        %whole year group not agreed to SCT
                                        isolgroup = all_students( testing_pop(yr,:,3) == 0);
                                        %whole year group for sct
                                        sctgroup = all_students(testing_pop(yr,:,3) == 1);
                                    else
                                        %no isolation or sct
                                        isolgroup = [];
                                        sctgroup = [];
                                     end
                                   
                                     
                                
                                 %Have always set Adherence.isolatingproperly == 0, so always satisfied 
                                 if sum(sum(history.IsolatingthroughCovid(:,:,day-1))) > 0 || r_ad2(yr,student,day) > Adherence.isolatingproperly(day)
                                     
                                   for otherpupil = isolgroup
                                           history.Isolation(yr, otherpupil, (day):(day+Isolation_period)) = 1; %isolates indivs not getting tested
                                   end
                                                
                                 end

                                    %change test days for those in SCT group
                                    for otherpupil = sctgroup
                                        history.gettingatest(yr,otherpupil, (day+1): day+Testing_period) = 2;
                                    end


                                    %do not isolate individual on that day if
                                    %leak_infect
                                   if leak_infect
                                       history.Isolation(yr, student, day) = 0;  %does not isolate on that day
                                    end

                                  end

                           end
                           
                      end
                     
                 end
                 
      
              end
             
       end
         
      
      %Susceptibles Not isolating during holiday weeks
      if any(Week == Infection.HolidayWeek)
            history.Isolation(:,:,day) = 0;
      end
      
      
      %Infections  - Tracking R
%{      
       Fromotheryears = zeros(1,YearGroup);

      
      %Within-school transmission
      if ~(Week == 1 ||day_of_week == 6 || day_of_week == 7 || any(Week == Infection.HolidayWeek))
          
          for yr = 1:YearGroup
              
              
              
              
              %find infectious pupils attending school
              j = find(Now.Infect(yr,:) > 0 & history.Isolation(yr,:,day) == 0);
             


              for student = j
                 
                  
                  
                  
                  %close contact infections            
                  closecontactprob = Now.Infect(yr,student)/(NumCloseContacts + alpha_withinyear*(YearSize-NumCloseContacts - 1) + alpha_betweenyears*(YearSize*(YearGroup-1)));
                  contactprob = closecontactprob;
                  
                  for contact = all_students((YearGroupMatrix(student,:) == 1 & history.Infection(yr,:,day) == 0 & history.Isolation(yr,:,day) == 0 ) & history.Infection(yr,:,day+1) == 0)
                            if rand < contactprob
                            %if rand < contactprob

                                history.Infection(yr,contact,day+1) = 1;
                                history.DayInf(yr,contact) = day;
                                history.Rday(1,day) = history.Rday(1,day)+1;
                                
                                history.ext_or_int(yr, contact, day+1) = 2;


                                if Symptomatic(yr, contact)
                                    %set day will seek and take PCR
                                    will_take_PCR(yr, contact) = day + day_of_PCR(yr, contact);                                                  
                                end

                                history.Rday(2, history.DayInf(yr,student)) = history.Rday(2, history.DayInf(yr,student)) + 1;


                            end
                      
                  end
                  
                  

                  
                  
                  %within year infections
                  withinyearprob = alpha_withinyear*closecontactprob;
                  contactprob = withinyearprob;
                  
                  for contact  = all_students(((YearGroupMatrix(student,:) == 0 & history.Infection(yr,:,day) == 0 & history.Isolation(yr,:,day) == 0 )) & history.Infection(yr,:,day+1) == 0)
                            if rand < contactprob
                            %if rand < contactprob

                                history.Infection(yr,contact,day+1) = 1;
                                history.DayInf(yr,contact) = day;
                                history.Rday(1,day) = history.Rday(1,day)+1;
                                
                                history.ext_or_int(yr, contact, day+1) = 2;


                                if Symptomatic(yr, contact)
                                    %set day will seek and take PCR
                                    will_take_PCR(yr, contact) = day + day_of_PCR(yr, contact);                                                  
                                end

                                history.Rday(2, history.DayInf(yr,student)) = history.Rday(2, history.DayInf(yr,student)) + 1;
                            end

                      
                  end


                                  
          %between year infections
             
                 if alpha_betweenyears > 0   
            
                  closecontactprob = Now.Infect(yr,student)/(NumCloseContacts + alpha_withinyear*(YearSize-NumCloseContacts - 1) + alpha_betweenyears*(YearSize*(YearGroup-1)));  
                  %rest of school infections
                  restofschoolprob = alpha_betweenyears*closecontactprob;
                  other_years = 1:YearGroup;
                    other_years(yr) = []; 
                 
                 for whichyear = other_years
                     
                     suspop = length(all_students((history.Infection(whichyear,:,day) == 0 & history.Infection(whichyear,:, day+1) == 0) & history.Isolation(whichyear,:,day) == 0));
                     
                     binominfect = mybinornd(suspop, restofschoolprob);
                                     
                     Fromotheryears(whichyear) = Fromotheryears(whichyear) + binominfect;                       
                     history.Rday(2, history.DayInf(yr,student)) = history.Rday(2, history.DayInf(yr,student)) + binominfect;
                     
                 end
                 
                 
                 end  
                  
                 end
               
          end
          
          %choose WHO gets infected from between year infections
          
          for whichyear = 1:YearGroup
                    
                     susgroup = all_students((history.Infection(whichyear,:,day) == 0 & history.Infection(whichyear,:, day+1) == 0) & history.Isolation(whichyear,:,day) == 0);
                           
                     actualinfecteds = randsample(susgroup, min(Fromotheryears(whichyear), length(susgroup)));
                     
                    for contact = actualinfecteds 

                             history.Infection(whichyear,contact,day+1) = 1;
                             history.DayInf(whichyear,contact) = day;
                             history.Rday(1,day) = history.Rday(1,day) + 1;
                             history.ext_or_int(whichyear,contact,day+1) = 2;

                            if Symptomatic(whichyear,contact)
                                %set day will seek and take PCR
                                will_take_PCR(whichyear,contact) = day + day_of_PCR(whichyear,contact);                                                  
                            end
                                            
                    end         
          end
                       
        end
          
      

            
      %now do external infections
      
      
      if any(Week == Infection.HolidayWeek)
          
          TodayExt  = Infection.HolidayExt*Ext;
      else
          TodayExt = Ext;
      end
      

      
      %mexed external infecctions      
      for yr = 1:YearGroup
          
           susgroup = all_students(history.Infection(yr,:,day) == 0 & history.Infection(yr,:,day+1) == 0);
           [histInfnextday, history.DayInf(yr,:), history.Rday, histExtorIntnextday, will_take_PCR(yr,:)] = Interactingyeargroups_externalinfection_mex(susgroup,  history.DayInf(yr,:),  Symptomatic(yr,:), will_take_PCR(yr,:), day_of_PCR(yr,:), history.Isolation(yr,:,day), history.Rday, rand(1,length(susgroup)), TodayExt, Ext_isolate, day);         
          
           %update infections
           history.Infection(yr,:,day+1) = 1*(histInfnextday|history.Infection(yr,:,day+1));           
           history.ext_or_int(yr,:,day+1) = max(histExtorIntnextday, history.ext_or_int(yr,:,day+1));
           
      end
      
      
     
      for yr = 1:YearGroup
          %update infection status of infected/recovered pupils
          infrecstudents = history.Infection(yr,:,day) > 0;
          history.Infection(yr, infrecstudents, day+1) = history.Infection(yr, infrecstudents, day) +1;
      end
          
%}      
    

      %%in use%% - not tracking R
     
      
            
      %Infection
      if any(Week == Infection.HolidayWeek)
     
          history.TotInfstudent(:,:,day) = Infection.HolidayExt*Ext;
      
      elseif Week == 1 ||day_of_week == 6 || day_of_week == 7 
            %weekend so no school infection or before kids are back to
            %school
            
            history.TotInfstudent(:,:,day) = Ext;
      
          if leak_infect %if infection happens on within-school test days
              
              if day_of_week == any(Strategy.initialtestdays)
                  %infection from each year group
                  for yr = 1:YearGroup
                      yg_inf(yr) = sum((Now.Infect(yr,:)).*(1 - history.Isolation(yr,:,day)));
                  end
                  
                  %infection to each student
                  for yr = 1:YearGroup
                      other_years = 1:YearGroup;
                      other_years(yr) = []; 
                      
                      %YearGroupMatrix = school_pop.school_matrix(((yr-1)*YearSize + 1): ((yr-1)*YearSize + YearSize), ((yr-1)*YearSize + 1): ((yr-1)*YearSize + YearSize));
                      %YearGroupMatrix = subMatrix(school_pop.school_matrix, yr, YearSize);
                  
                      
                 %{     
                      
                  if Week > Infection.HolidayWeek(end)

                        YearGroupMatrix = subMatrix(school_pop2.school_matrix, yr, YearSize);
                    else
                        YearGroupMatrix = subMatrix(school_pop.school_matrix, yr, YearSize);
                  end
                 %}   
                      closecontact_infvec = (Now.Infect(yr,:).*(1-history.Isolation(yr,:,day)))*YearGroupMatrix;
                      history.TotInfstudent(yr,:,day) =  Ext + (closecontact_infvec+alpha_withinyear*(yg_inf(yr) - closecontact_infvec)+alpha_betweenyears*sum(yg_inf(other_years)))/(NumCloseContacts + alpha_withinyear*(YearSize - NumCloseContacts - 1) + alpha_betweenyears*length(other_years)*YearSize); 
                  end
              
              end
              
          end
          
      else
          
      %infection from each year group
              for yr = 1:YearGroup
                  yg_inf(yr) = sum((Now.Infect(yr,:)).*(1 - history.Isolation(yr,:,day)));
              end

              %infection to each student
              for yr = 1:YearGroup
                  other_years = 1:YearGroup;
                  other_years(yr) = []; 

                  %YearGroupMatrix = school_pop.school_matrix(((yr-1)*YearSize + 1): ((yr-1)*YearSize + YearSize), ((yr-1)*YearSize + 1): ((yr-1)*YearSize + YearSize));
                  %YearGroupMatrix = subMatrix(school_pop.school_matrix, yr, YearSize);
                  
                  closecontact_infvec = (Now.Infect(yr,:).*(1-history.Isolation(yr,:,day)))*YearGroupMatrix;
                  history.TotInfstudent(yr,:,day) =  Ext + (closecontact_infvec+alpha_withinyear*(yg_inf(yr) - closecontact_infvec)+alpha_betweenyears*sum(yg_inf(other_years)))/(NumCloseContacts + alpha_withinyear*(YearSize - NumCloseContacts - 1) + alpha_betweenyears*length(other_years)*YearSize); 
              end
      end
      
      %}
      

   

      
%%In use%% - not tracking R


      %[history.Infection, will_take_PCR] = Interactingyeargroups_Infection_mex(history.Infection, history.TotInfstudent, history.Isolation,  will_take_PCR, Symptomatic, day,Ext_isolate, r_inf, day_of_PCR);
      [InfNextday, will_take_PCR] = Interactingyeargroups_Infection_mex(squeeze(history.Infection(:,:,day)), squeeze(history.TotInfstudent(:,:,day)), squeeze(history.Isolation(:,:,day)),  will_take_PCR, Symptomatic, day,Ext_isolate, squeeze(r_inf(:,:,day)), day_of_PCR);
      
      history.Infection(:,:,day+1) = InfNextday;
%}      
      

        
        
        %Updates  
        if ~isempty(Infection.HolidayWeek) && Week < Infection.HolidayWeek(1)
          %  Now = Update2(Now, history, Prob_profiles, Symptomatic, day+1, Ks);
            
            Now = Update3(Now, history.Infection(:,:,(day+1)), Prob_profiles, Symptomatic, day+1, Ks);            
        elseif isempty(Infection.HolidayWeek)
            %Now = Update2(Now, history, Prob_profiles, Symptomatic, day+1, Ks);
             Now = Update3(Now, history.Infection(:,:,(day+1)), Prob_profiles, Symptomatic, day+1, Ks);
        elseif ~isempty(Infection.HolidayWeek) && Week < Infection.HolidayWeek(end-1)
            %Now = Update2(Now, history, Prob_profiles, Symptomatic, day+1, Ks2);
             Now = Update3(Now, history.Infection(:,:,(day+1)), Prob_profiles, Symptomatic, day+1, Ks2);
        else
            % Now = Update2(Now, history, Prob_profiles, Symptomatic, day+1, Ks3);
             Now = Update3(Now, history.Infection(:,:,(day+1)), Prob_profiles, Symptomatic, day+1, Ks3);

        end
        
        
        %Quantities through time
        for yr = 1:YearGroup
            %Isolated infecteds       
            history.Isolated_Infecteds(yr,day) = sum(history.Isolation(yr, history.Infection(yr,:,day) > 0 & history.Infection(yr,:,day) < 16, day));
            history.Tot_Isolated(yr,day) = sum(history.Isolation(yr,:,day));
        end   
        %}    
    end
    
end


history.Infection = history.Infection(:,:,1:(Weeks*7));
history.Isolation = history.Isolation(:,:,1:(Weeks*7));



history.Isolated_Infecteds = history.Isolated_Infecteds(:, 1:(Weeks*7));


history.Symptomatic = Symptomatic;
history.gettingatest = history.gettingatest(:,:, 1:(Weeks*7));
%}
history.IsolatingthroughCovid = history.IsolatingthroughCovid(:,:,(1:Weeks*7));
history.pos_test_day = history.pos_test_day(:,1:(Weeks*7));

%to compare with previous version
history.day_of_PCR = day_of_PCR;
history.r_lat = r_lat;
history.r_PCR = r_PCR;
%}

history.YearGroups = YearGroup;
history.HolidayWeek = Infection.HolidayWeek;

end
    
function Now = Update2(Now, history, Prob_profiles, Symptomatic, Day,Ks)
    %Symptomatics need to be transposed when 1 year group 
    %{
    Now.Infect = (Symptomatic'.*Probs.Inf_sym(history.Infection(:,:,Day) +1) + (1-Symptomatic)'.*Probs.Inf_asym(history.Infection(:,:,Day) +1))';
    Now.PCR  =   (Symptomatic'.*Probs.PCR_sym(history.Infection(:,:,Day) +1) + (1-Symptomatic)'.*Probs.PCR_asym(history.Infection(:,:,Day) + 1))';
    Now.lat  =   (Symptomatic'.*Probs.lat_sym(history.Infection(:,:,Day) +1) + (1-Symptomatic)'.*Probs.lat_asym(history.Infection(:,:,Day) + 1))';
    %}
    %But not with larger year groups? figure out
    Now.Infect = Ks.*(Prob_profiles.newvar(Day)*(Symptomatic.*Prob_profiles.Inf_symday(history.Infection(:,:,Day) +1) + (1-Symptomatic).*Prob_profiles.Inf_asymday(history.Infection(:,:,Day) +1)));
   % Now.Infect = (1-history.Isolation(:,:,Day)).*Now.Infect;
    Now.PCR  =   (Symptomatic.*Prob_profiles.PCR_symday(history.Infection(:,:,Day) +1) + (1-Symptomatic).*Prob_profiles.PCR_asymday(history.Infection(:,:,Day) + 1));
    Now.lat  =   (Symptomatic.*Prob_profiles.lat_symday(history.Infection(:,:,Day) +1) + (1-Symptomatic).*Prob_profiles.lat_asymday(history.Infection(:,:,Day) + 1));
    
    
end

function Now = Update3(Now, histInf, Prob_profiles, Symptomatic, Day, Ks)
    Now.Infect = Ks.*(Prob_profiles.newvar(Day)*(Symptomatic.*Prob_profiles.Inf_symday(histInf +1) + (1-Symptomatic).*Prob_profiles.Inf_asymday(histInf +1)));
   % Now.Infect = (1-history.Isolation(:,:,Day)).*Now.Infect;
    Now.PCR  =   (Symptomatic.*Prob_profiles.PCR_symday(histInf +1) + (1-Symptomatic).*Prob_profiles.PCR_asymday(histInf + 1));
    Now.lat  =   (Symptomatic.*Prob_profiles.lat_symday(histInf +1) + (1-Symptomatic).*Prob_profiles.lat_asymday(histInf + 1));
    
    
end

function ret = mybinornd( n, p )
  ret = sum(rand(1,n)<p);
end
    
function subM = subMatrix(M, whichyear, YearSize)
            subM = M(((whichyear-1)*YearSize + 1): ((whichyear-1)*YearSize + YearSize), ((whichyear-1)*YearSize + 1): ((whichyear-1)*YearSize + YearSize));
end
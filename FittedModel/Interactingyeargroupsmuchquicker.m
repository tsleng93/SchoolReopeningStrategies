function [history] = Interactingyeargroupsquicker(school_pop, Infection, Testing, Strategy, Adherence, Prob_profiles, school_pop2)
%This function simulates the spread of infection between pupils within a
%secondary school formed of close contact groups in within year groups.

%This version of the model code is significantly faster than the version
%that tracks R_school. Speed gains are made by setting external infections
%at the beginning of the simulation, only simulating infection within
%schools if a pupil is infected, and by calculating the force of infection
%to each pupil from all other pupils (rather than working out who infects
%who)


%% --- Inputs --- %%
if nargin == 0
% - school_pop and school_pop2 are structs containing the close contact
    [school_pop] = SchoolPopulationquicker(200, 5, 25);
    school_pop2 = school_pop;

%   structure of year groups (.school_matrix) , the number of year groups (.YearGroups) and the size of each
%   year group (.YearSize). It is assumed each year group has the same close
%   contact structure, i.e. the same school_matrix is used for each year group

%- Infection is a struct containing parameters relevant to transmission:
        
        %used in paper
                
        %K (> 0) %secondary infections expected from symptomatic individual 
                 %over course of their infection (without depletion of 
                 %susceptibles)
        Infection.K = 1;                         
        %K2 (> 0) % value of K after October break
        Infection.K2 = 1;
        %Inf_0 (0-1) %Initial population prevalence
        Infection.Inf_0 = 0;
        %Rec_0 (0-1) %Initial population immunity
        Infection.Rec_0 = 0;
        %K_asym (0-1) Relative infectiousness of asymptomatic individuals
        Infection.K_asym = 0.5;
        %Sym_proportion (0-1) %Proportion of school population who are 
                              %symptomatic
        Infection.Sym_proportion = 0.3;
        %alpha_withinyear (0-1) %interaction with other members of the year 
                                %group (not close contacts)
        Infection.alpha_withinyear = 0.1;
        %alpha_betweenyears (0-1) %interaction between year groups
        Infection.alpha_betweenyears = 0.01;
        %Weeks (Positive integer) %Number of weeks (including week before 
                                  %term)
        Infection.Weeks = 39;
        %HolidayWeek (Subset of Weeks) %Weeks school is on holiday/not open
        Infection.HolidayWeek =  [10, 18:28, 32, 33];
        %Ext (> 0) %External transmission scaling constant (i.e. eps from 
                   % paper)
        Infection.Ext = 1;       
        %HolidayExt (> 0) % sclaing factor on Ext during school holidays
        Infection.HolidayExt = 1.5;        
        %commext (non-negative vector of length Weeks) %time-varying
                 %community infection (informed by community testing data)
        Infection.commext = 0.002*ones(1,Infection.Weeks*7);
        
        %not used in paper
        
        %K3 (> 0) value of K after schools return in March 
        Infection.K3 = 1;
        %leak_infect (0 or 1)  whether infection occurs on test days
        Infection.leak_infect = 0;

%- Testing is a struct containing parameters relevant to testing:
    
    %used in paper
    %sens_PCR ({1,2,3}) %the column of PCR sensitivity profile to read 
    Testing.sens_PCR = 1;
    %sens_lat ({1,2,3}) %the column of LFT sensitivity profile to read
    Testing.sens_lat = 1;
                        %(1) - baseline sensitivity
                        %(2) - low sensitivity
                        %(3) - high sensitivity
   %spec_PCR (0-1) %specificity of PCR tests
   Testing.spec_PCR = 1;
   %spec_lat (0-1) %specificity of LFT tests
   Testing.spec_lat = 0.9997;
   %PCR_delay (Positive integer) % delay on PCR tests
   Testing.PCR_delay = 2;   
   
%- Strategy is a struct containing parameters relevant to the strategy 
    % implemented in schools from Week 27
    % It is assumed until Week 27 schools implement an isolation of close
    % contacts policy without mass testing

    %N.B. if an isolation policy is implemented, this is assumed to override
    %the SCT policy

    %used in paper
    %isolation - (0) - not isolating, (1) - isolating close contacts, (2) - isolating year groups
    Strategy.isolation = 1;
    %SCT - (0) - no serial contact testing, (1) - serial contact testing close contacts, (2) - serial contact testing year groups 
    Strategy.SCT = 0;
    %masstesting - %(0) - no mass testing, (1) - weekly mass testing, (2) - twice weekly mass testing
    Strategy.masstesting = 2;
    
    %not used in paper
    %initialtestdays (Subset of 1:Weeks*7) %initial days pupils are tested
    Strategy.initialtestdays = [];
    
%- Adherence is a struct containing parameters relevant to adherence

    %used in paper
    %probtakelatflow % vector of daily values of the expected proportion of
                     % the school population who take an LFT each dsay
    Adherence.probtakelatflow = 0.5*ones(1,Infection.Weeks*7);
    
    %not used in paper
    % takenproperly (0-1) % of LFTs taken properly (scaling sensitivity)
    Adherence.takenproperly = 0;
    
%- Prob_profiles is a struct containing probability profiles used in the model
   %also contains the relative frequency of Alpha variant through time
   
   
       %Probability Profiles
    PCR_test_sym = readtable('PCR_Curve_summary.csv');
    PCR_test_sym = table2array(PCR_test_sym(:, 2:4));
    PCR_test_asym = csvread('PCR_Curve_asym.csv');
    lat_test_sym = readtable('lat_Curve_summary.csv');
    lat_test_sym = table2array(lat_test_sym(:, 2:4));
    lat_test_asym = csvread('lat_Curve_asym.csv');

  %used in paper
  %PCRsym - prob symptomatic individual tests positive to a PCR on day d of
            %their infectious period
  Prob_profiles.PCRsym = PCR_test_sym;
  %PCRasym - prob asymptomatic individual tests positive to a PCR on day d 
            %of their infectious period
  Prob_profiles.PCRasym = PCR_test_asym;
  %latsym - prob symptomatic individual tests positive to a PCR on day d of
            %their infectious period
  Prob_profiles.latsym = lat_test_sym;
  %latasym - prob asymptomatic individual tests positive to a PCR on day d 
            %of their infectious period
    Prob_profiles.latasym = lat_test_asym;
            
  %Infectivity - infectiousness profile through time
    Infectivity_since_infection = [0.0063 0.0563 0.1320 0.1798 0.1817 0.1521 0.1117 0.0746 0.0464 0.0272 0.0152 0.0082 0.0043 0.0022 0.0011 0.0005 0.0002 0.0001 0.0001 0.0000];
    Infectivity_since_infection = Infectivity_since_infection'/sum(Infectivity_since_infection);
    Infectivity_since_infection(15) = sum(Infectivity_since_infection(15:end)); Infectivity_since_infection = Infectivity_since_infection(1:15);
    Prob_profiles.Infectivity = Infectivity_since_infection;
  %Symptom_onset - incubation period (time to symptom onset) profile
                   %through time
    Symptom_onset =  [0.0055 0.0534 0.1307 0.1814 0.1847 0.1545 0.1129 0.0747 0.0458 0.0265 0.0146 0.0077 0.0039 0.0020 0.0010 0.0005 0.0002 0.0001 0 0];
    Symptom_onset = Symptom_onset/sum(Symptom_onset);
    Symptom_onset(15) = sum(Symptom_onset(15:end)); Symptom_onset = Symptom_onset(1:15);
    Prob_profiles.Symptom_onset = Symptom_onset;                 
                   
  %new_var - Relative frequency of new variant on given day
  Prob_profiles.newvar = ones(1,Infection.Weeks*7 +1);          
end
           

%% ----Output ---- 
%   - history is a structure, containing the infection, testing, and
%   isolation details throughout the course of the simulation. Useful
%   outputs from this can be obtained using Modeloutputscondensed.m 


%Authors: Trystan Leng and Edward M. Hill 
%Last update 11/11/2021.
%% ---- Parameters ----- 

%Directly from input

%get infection parameters used often
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
Prob_profiles.PCR_symday(Prob_profiles.PCR_symday < (1- Testing.spec_PCR)) = (1-Testing.spec_PCR);
Prob_profiles.PCR_asymday(Prob_profiles.PCR_asymday < (1-Testing.spec_PCR)) = (1-Testing.spec_PCR);

%Adjust for taken properly
Prob_profiles.lat_symday = Adherence.takenproperly*Prob_profiles.lat_symday;
Prob_profiles.lat_asymday = Adherence.takenproperly*Prob_profiles.lat_asymday;




YearSize =  school_pop.yearsize;
YearGroup = school_pop.yeargroups;
Ks = Infection.K*ones(YearGroup,YearSize);


%%Infectivity profile%%
Prob_profiles.Inf_symday = [0;Prob_profiles.Infectivity; zeros(max(51,Weeks*7-15 + 50),1)];
Prob_profiles.Inf_asymday = Infection.K_asym*Prob_profiles.Inf_symday;

%% ---- storing history ----


%Keep entire history
history.Infection = zeros(YearGroup, YearSize, Weeks*7+10); %Infection status of individuals through time
history.Isolation = zeros(YearGroup, YearSize, Weeks*7+10); %Whether individuals are isolated or not through time
history.pos_test_day = zeros(YearGroup, Weeks*7+10); %Number of positive tests that day
history.posLFT = zeros(YearGroup, YearSize, Weeks*7+10); %Whether an individual tests positive to an LFT that day
%history.posLFTandPCR = zeros(YearGroup, YearSize, Weeks*7+10);
history.IsolatingthroughCovid = zeros(YearGroup, YearSize,Weeks*7 +10); %whether an individual is isolating because they are confirmed covid positive
history.gettingatest = zeros(YearGroup, YearSize, Weeks*7+10); %Whether an individual is tested that day %2 for SCT, 1 for LFT 


%not actually a history
history.TotInfstudent = zeros(YearGroup,YearSize); %force of infection to each student (not stored as a history)


%% ---- Things that depend on randomness ----

%testing day
r_testingday = randi(5,YearGroup,YearSize) -1;

%Set who is symptomatic
Symptomatic = 1*(rand(YearGroup,YearSize) < Infection.Sym_proportion);

%day an individual will seek PCR test if symptomatic i.e. duration of
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
                     if rinit(end) <  Prob_profiles.PCR_symday(day_of_PCR(yr, student))
                         %tests positive before simulation begins
                         history.Isolation(yr, student, 1: (1 + 8 + will_take_PCR(yr,student))) = 1; %if tests is on day 0, then student should isolate  from day 1 to day 9                                                                                            %if test is on day -1, then student should isolate from day 1 to day 8 
                        % history.ever_test_pos(yr, student) = 1;
                     end
                 end

            end  
                rinit(end) = [];
        end
   end
end


%day individual will be infected from population at large
for i = 1:Infection.Weeks
    for j = 1:7
    day = (i-1)*7 + j;
        if any(Infection.HolidayWeek == i)
            Holiday(day) = Infection.HolidayExt;
        else
            Holiday(day) = 1;
        end
    end
end

%logical of whether someone is infected externally on a given day
ExtInfectLogical = 1*(rand(YearSize,Weeks*7, YearGroup) < Infection.Ext*(Holiday.*Infection.commext));
ExtInfectLogical = permute(ExtInfectLogical, [3 1 2]); %put in format (YearGroup, YearSize, Weeks*7)


%set populations who agree to being tested
%testing_pop = zeros(YearGroup, YearSize, 4); %1 for initial tests, 2 for mass tests, 3 for SCTs, 4 for adhering to mass tests
%testing_props = [ Strategy.initialuptake,  1; Strategy.masstestinguptake, 2; Strategy.SCTuptake, 3];
%testing_orders = sortrows(testing_props);
testing_pop = ones(YearGroup,YearSize,4); %Set everyone to take tests

%set days of week pupils seek home LFTs
rdays = zeros(5, Weeks*7+10);
for rday = 0:4
    for day = 1:Weeks*7
        if mod(day,7) == rday || mod(day,7) == mod(rday+3, 7) 
               rdays(rday+1,day) = 1;
        end
    end
end



%Used in longer version - but in shorter version we assume all agree to be
%tested
%{

%Assume that the min of the 3 agree to all
for yr = 1:YearGroup   
    all_students = 1:YearSize;
    core_group = randsample(all_students, round(testing_orders(1,1)*YearSize));    
    %set core group who get tested to all
    testing_pop(yr, core_group, testing_orders(1,2)) = 1;    
    testing_pop(yr,:, testing_orders(2,2)) = testing_pop(yr,:, testing_orders(1,2));
    
    remaining_students = all_students; 
    remaining_students(core_group) = [];

    next_group = randsample(remaining_students, min(length(remaining_students), round((testing_orders(2,1) - testing_orders(1,1))*YearSize)));
    
    %set next group who get tested to two of the tests
    testing_pop(yr, next_group, testing_orders(2,2)) = 1; 
    testing_pop(yr,:, testing_orders(3,2)) = testing_pop(yr,:, testing_orders(2,2));

    remaining_students = all_students; 
    remaining_students([core_group, next_group]) = [];
    
    final_group = randsample(remaining_students, min(length(remaining_students), round((testing_orders(3,1) - testing_orders(2,1))*YearSize)));
    
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
%}

%set test days for mass tests for adhering students
if Strategy.masstesting == 1
  % for day = 8:(Weeks*7)
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
    
    %{
    %3 years tested monday and thursday, 2 years tested tuesday and friday
    for day = 8:Weeks*7
      for yr = 1:YearGroup
          if yr <= round(0.5*YearGroup)
              if mod(day, 7) == 1 || mod(day,7) == 4
                  history.gettingatest(yr,testing_pop(yr,:,4) == 1,day) = 1;
                  history.gettingatest(yr,testing_pop(yr,:,4) == 1,day) = 1;
              end
          else
              if mod(day, 7) == 2 || mod(day,7) == 5
                  history.gettingatest(yr,testing_pop(yr,:,4) == 1,day) = 1;
                  history.gettingatest(yr,testing_pop(yr,:,4) == 1,day) = 1;
              end
          end
      end
    end
    %}
    
    %change so tests occur on random days
    %Mon-Thur, Tue-Fri, Wed-Sat, Thurs-Sun, Fri-Mon
    
    %{
    for yr = 1:YearGroup
        for student = find(testing_pop(yr,:,4) == 1)
            r = r_testingday(yr, student);
         
        %{    
          %  for day = 8:Weeks*7
            for day = 1:Weeks*7
               if mod(day,7) == r || mod(day,7) == mod(r+3, 7) 
                   history.gettingatest(yr,student,day) = 1;
               end
            
            end
         %}
            history.gettingatest(yr,student,:) = rdays(r+1,:);
            
            
            
        end
    end
   %}
   
    
   %%In use - comment out for bugfinding 
   
   gettingatest1 = rdays(r_testingday+1, :);
   history.gettingatest = reshape(gettingatest1, YearGroup, YearSize, 283);
    %}
end




%find proportion of students being tested each day
prop_tested = squeeze(sum(sum(history.gettingatest > 0)))/(YearSize*YearGroup);


%prop_tested = squeeze(sum(sum(history.gettingatest > 0)));


%set test days for initial opt-ins (overrides mass test days so has to be
%set after)
for yr = 1:YearGroup
   history.gettingatest(yr, testing_pop(yr,:,1) == 1, Strategy.initialtestdays) = 2;   
end

%Sample which tests are actually taken
 %%In use - comment out for bugfinding
actuallytakentest1 = zeros(YearGroup,YearSize,Weeks*7);

%sample whether an individual takes a test on a given day, set to match the
%proportion of the school population taking an LFT that day
actuallytakentest1(:,:,(27*7+1:end)) = 1*(rand(YearGroup, YearSize, (Weeks-27)*7) < reshape((Adherence.probtakelatflow((27*7+1):end)./prop_tested((27*7+1):273)'), 1,1, []));
history.actuallytakentest = actuallytakentest1.*history.gettingatest(:,:,1:273);



%Individual's transmission probability and testing probabilities
Now.Infect = zeros(YearGroup, YearSize); 
Now.had_PCR = zeros(YearGroup, YearSize); %have they had a +ve PCR?
%Not stored like this in quicker version
%{
Now.PCR = zeros(YearGroup, YearSize);
Now.lat = zeros(YearGroup, YearSize);
%}


%Initial test and infection probabilities 
Now = Update2(Now, history, Prob_profiles, Symptomatic, 1, Ks);


%define all_students vec
all_students = 1:YearSize;



history.closecontacttestpos = zeros(YearGroup, YearSize, Weeks*7);
%history.closecontacts = zeros(YearGroup, YearSize);

%Chosen isolation, mass testing, and SCT strategies
Stratisolinit = Strategy.isolation;
Stratmasstestinit = Strategy.masstesting;
StratSCTinit = Strategy.SCT;

%% ----- Simulation ------ %%
for Week = 1:Weeks
       
   %Strategy of isolating close contacts before then
    if Week < 28

        YearGroupMatrix = school_pop.school_matrix;
        NumCloseContacts = sum(YearGroupMatrix(1,:));
       
        Strategy.isolation = 1;
        Strategy.masstesting = 0;
        Strategy.SCT = 0;
       
    else
    %Inputted strategy afterwards,    
        YearGroupMatrix = school_pop2.school_matrix;
        NumCloseContacts = sum(YearGroupMatrix(1,:));
        
        Strategy.isolation = Stratisolinit;
        Strategy.masstesting = Stratmasstestinit;
        Strategy.SCT = StratSCTinit;

    end

    for day_of_week = 1:7
                      
        day = (Week-1)*7 + day_of_week;
        
        
       if any(Week == Infection.HolidayWeek)
           
          %don't count this during holidays
          history.IsolatingthroughCovid(:,:,day) = 0;

       end
        
        


            
     %  Ext = Infection.Ext*Infection.commext(day);
     %  Ext_isolate = Adherence.C_isolate*Ext;
       
     %% --- Testing students via PCR (self-seeking) --- %%
       
   if sum(sum(will_take_PCR == day)) > 0  % see if any students seek PCR that day before doing loop 
        for yr = 1:YearGroup
           if sum(will_take_PCR(yr,:) == day) > 0 % see if any students seek PCR that day before doing loop 
            
            for student = all_students(will_take_PCR(yr,:) == day)
                  
                    if Now.had_PCR(yr, student) == 0  %student does not seek PCR if already tested positive to a PCR
                        %student seeks test
                        if Symptomatic(yr,student)
                           nPCR = Prob_profiles.PCR_symday(history.Infection(yr,student,day) +1); %prob of testing +ve if symptomatic
                        else
                           nPCR = Prob_profiles.PCR_asymday(history.Infection(yr,student,day) +1); %prob of testing +ve is asymptomatic
                        end
                        
                        if rand < nPCR
                            %student tests positive
                            Isolation_period = 10;  %Isolates for 10 days from symptom onset                          
                            notisolmarker = 0; %not already isolating 
                            Now.had_PCR(yr, student) = 1; %has had a +ve PCR so set this to 1
                           
                            if history.Isolation(yr,student,day) == 0
                                notisolmarker = 1; %update this to 1 if already isolating
                            elseif history.Isolation(yr,student,day) == 1 || history.gettingatest(yr,student,day) == 2
                                history.closecontacttestpos(yr,student,day) = 1; %a 'close contact' tests +ve
                            end
                            
                            history.Isolation(yr, student, day:(day+Isolation_period)) = 1; %isolates for Isolation_period days
                            history.IsolatingthroughCovid(yr,student,(day+PCR_delay):(day + Isolation_period)) = 1; %isolating because confirmed covid positive
                            history.pos_test_day(yr, day) = history.pos_test_day(yr, day) +1; %positive PCR test that day so increase counter

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
                                
                                %isolate other individuals in group 
                                for otherpupil = isolgroup
                                   history.Isolation(yr, otherpupil, (day+PCR_delay+1):(day+Isolation_period)) = 1; %isolates indivs not getting tested
                                end                               

                                %change test days for those in SCT group
                                for otherpupil = sctgroup
                                    history.gettingatest(yr,otherpupil, (day+PCR_delay + 1): day+7) = 2;
                                    history.actuallytakentest(yr,otherpupil, (day+PCR_delay + 1): day+7) = 2;
                                end                         
                            end                            
                                                                                 
                        else
                            %student tests negative
                            history.Isolation(yr, student, day:(day+PCR_delay)) = 1; %isolates until PCR returns
                        end
                    end
            end
            
           end
        end
        
    end
         
       
       
       
        

         
        if Week > 27
        
         %% --- mass testing using LFTs -----
            
              for yr = 1:YearGroup
                 
              % find group who are mass testing -> those who actually take
              % a test that day but have not already taken a PCR
              masstestinggroup = all_students(history.actuallytakentest(yr,:,day) > 0 & will_take_PCR(yr,:) ~= day);
                 
                 for student = masstestinggroup

                    if history.Isolation(yr,student,day)
                        nLat = 0; %assume isolating individuals do not take twice weekly LFTs
                    else
                        if Symptomatic(yr,student)
                           nLat = Prob_profiles.lat_symday(history.Infection(yr,student,day) +1); %prob of testing +ve if symptomatic
                        else
                           nLat = Prob_profiles.lat_asymday(history.Infection(yr,student,day) +1); %prob of testing +ve if asymptomatic
                        end 
                    end
                     
                            if rand < nLat
                                %student tests positive to an LFT
                                
                             history.posLFT(yr,student,day) = 1;  %record student test positive to an LFT that day 

                               %confirmatory PCR
                                 if history.Infection(yr, student, day) == 0 || history.Infection(yr, student, day) > 15
                                      %false LFT positive, but recovered, so may
                                      %test positive by chance  
                                      
                                    if Symptomatic(yr,student)
                                       nPCR = Prob_profiles.PCR_symday(history.Infection(yr,student,day) +1); %prob of testing +ve if symptomatic
                                    else
                                       nPCR = Prob_profiles.PCR_asymday(history.Infection(yr,student,day) +1); %prob of testing +ve if asymptomatic
                                    end                                     
                                      
                                     if  rand < nPCR
                                          %tests positive - so close
                                          %contacts will isolate and those
                                          %agreeing to be tested will test
                                          Isolation_period = 10; 
                                          Testing_period = 7;                                 
                                          %Assume they will not seek another PCR test
                                          Now.had_PCR(yr,student) = 1; 

                                          history.IsolatingthroughCovid(yr,student,(day+PCR_delay):(day + Isolation_period)) = 1; %isolating because confirmed covid positive
                                                                     
                                            if history.Isolation(yr,student,day) == 1 || history.gettingatest(yr,student,day) == 2
                                                history.closecontacttestpos(yr,student,day) = 1; %a 'close contact' tests +ve
                                            end
                                     
                                     else
                                         %tests negative -> testing and
                                         %isolation only goes on until
                                         %negative result is returned
                                         Isolation_period = PCR_delay; 
                                         Testing_period = PCR_delay;
                                     end

                                  else
                                      %true positive so we assume confirmatory PCR comes
                                      %back positive
                                      Isolation_period = 10;
                                      Testing_period = 7;                                 
                                      %Assume they will not seek another PCR test
                                      Now.had_PCR(yr,student) = 1; 
                                                                            
                                        if history.Isolation(yr,student,day) == 1 || history.gettingatest(yr,student,day) == 2
                                            history.closecontacttestpos(yr,student,day) = 1; %a 'close contact' tests +ve
                                        end
 
                                       history.IsolatingthroughCovid(yr,student,(day+PCR_delay):(day + Isolation_period)) = 1;  %isolating because confirmed covid positive

                                  end


                                  history.Isolation(yr, student, (day):(day+Isolation_period)) = 1; %isolates +ve individual

                                  %set isolation groups and SCT groups
                                  %according to strategy
                                  if day > 8  && ~any(Week == Infection.HolidayWeek)
                                  %if day > 7
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
                                   
                                   %isolate other individuals in group 
                                   for otherpupil = isolgroup
                                           history.Isolation(yr, otherpupil, (day):(day+Isolation_period)) = 1; %isolates indivs not getting tested
                                   end

                                    %change test days for those in SCT group
                                    for otherpupil = sctgroup
                                        history.gettingatest(yr,otherpupil, (day+1): day+Testing_period) = 2;
                                        history.actuallytakentest(yr,otherpupil, (day+1): day+Testing_period) = 2;
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
  
      %Assume pupils do not isolate during holiday weeks 
      %Not isolating during holiday weeks
      if any(Week == Infection.HolidayWeek)
            history.Isolation(:,:,day) = 0;
      end
      
%% ---- Transmission ----


%much quicker version
%% ---- infections within school --- 
if ~(any(Week ==  Infection.HolidayWeek) || (Week == 1 || day_of_week == 6 || day_of_week == 7))
    %if pupils are in school
    
    if sum(sum(history.Infection(:,:,day) > 0 & history.Infection(:,:,day) < 16)) > 0
    %only do next bits if any students are actually infected     
        
       %Updates  
        if ~isempty(Infection.HolidayWeek) && Week < Infection.HolidayWeek(1)
          %  Now = Update2(Now, history, Prob_profiles, Symptomatic, day+1, Ks);
            
            Now = Update4(Now, history.Infection(:,:,(day)), Prob_profiles, Symptomatic, day, Infection.K, history.Isolation(:,:,day));            
        elseif isempty(Infection.HolidayWeek)
            %Now = Update2(Now, history, Prob_profiles, Symptomatic, day+1, Ks);
             Now = Update4(Now, history.Infection(:,:,(day)), Prob_profiles, Symptomatic, day, Infection.K, history.Isolation(:,:,day));
        elseif ~isempty(Infection.HolidayWeek) && Week < Infection.HolidayWeek(end-1)
            %Now = Update2(Now, history, Prob_profiles, Symptomatic, day+1, Ks2);
             Now = Update4(Now, history.Infection(:,:,(day)), Prob_profiles, Symptomatic, day, Infection.K2, history.Isolation(:,:,day));
        else
            % Now = Update2(Now, history, Prob_profiles, Symptomatic, day+1, Ks3);
             Now = Update4(Now, history.Infection(:,:,(day)), Prob_profiles, Symptomatic, day, Infection.K3, history.Isolation(:,:,day));

        end
        
        %if anyone is actually infected


        %{
        
              yg_inf = sum(Now.Infect, 2);  %force of infection from year group           
              rs_inf = sum(yg_inf) - yg_inf; %force of infection from rest of school to each year group

              
              closecontact_infmat = Now.Infect*YearGroupMatrix; %infection to each student
              
              
              %force of infection to each student (before accounting for if
              %they are isolating)
              history.TotInfstudentorig = (closecontact_infmat+ alpha_withinyear*(yg_inf - closecontact_infmat) + alpha_betweenyears*rs_inf)/(NumCloseContacts + alpha_withinyear*(YearSize - NumCloseContacts - 1) + alpha_betweenyears*(YearGroup-1)*YearSize);
           %}
        
               const = (NumCloseContacts + alpha_withinyear*(YearSize - NumCloseContacts - 1) + alpha_betweenyears*(YearGroup-1)*YearSize);
               
               
               for yr = 1:YearGroup
                   closecontact_infmat(yr,:) = prod(1 - (YearGroupMatrix.*Now.Infect(yr,:)')/const);
                   yg_infmat(yr,:) = prod(1 - (alpha_withinyear*(1-YearGroupMatrix).*Now.Infect(yr,:)')/const);
               end
               
               %closecontact_infmat = Now.Infect*YearGroupMatrix;
               
               
               %yg_infmat = Now.Infect*(1-YearGroupMatrix);
               rs_inf = (1 - (prod(prod(1 - (alpha_betweenyears/const)*Now.Infect)))./prod(1 - (alpha_betweenyears/const)*Now.Infect'))';
               rs_inf(isnan(rs_inf)) = 0;
               
               history.TotInfstudent = 1 - (closecontact_infmat).*(yg_infmat).*(1-rs_inf);
        
        
             %Set force of infection to isolating students to be 0 
             isolating = find(history.Isolation(:,:,day));                 
             history.TotInfstudent(isolating) = 0;    
              
          %Loop through to see who is infected and update will_take_PCR   
          [InfNextday, will_take_PCR] = Interactingyeargroups_Infection_mex(squeeze(history.Infection(:,:,day)), history.TotInfstudent,  will_take_PCR, Symptomatic, day, day_of_PCR);
          history.Infection(:,:,day+1) = InfNextday;
          
    else
        %still need to do updates
        histInftoday =  history.Infection(:,:,day);
        histInfnextday = history.Infection(:,:,day+1);          
        %update day of infection for infected/recovered individuals
        histInfnextday(histInftoday > 0) = histInftoday(histInftoday > 0) +1;
        history.Infection(:,:,day+1) = histInfnextday;
    end
    
else
    %still need to do updates
    histInftoday =  history.Infection(:,:,day);
    histInfnextday = history.Infection(:,:,day+1);
    %update day of infection for infected/recovered individuals
    histInfnextday(histInftoday > 0) = histInftoday(histInftoday > 0) +1;
    history.Infection(:,:,day+1) = histInfnextday;
end

%% ------ infections from external sources ----- 


if sum(sum(ExtInfectLogical(:,:,day))) > 0 
    %only do this if anybody is infected from external sources that day
   
    %find who becomes infected
    extinfected = find(ExtInfectLogical(:,:,day));
    
    histInftoday = history.Infection(:,:,day);
    histIsoltoday = history.Isolation(:,:,day);
    histInfnextday = history.Infection(:,:,day+1);
    
    for i = extinfected'
              
       if histInftoday(i) == 0 & histIsoltoday(i) == 0 
       %pupil becomes infected if not already infected and not isolating
           histInfnextday(i) = 1;
           
           if Symptomatic(i)
               %set day they will take PCR
                will_take_PCR(i) = day + day_of_PCR(i);
           end
       end
    end
 %update infection   
 history.Infection(:,:,day+1) = histInfnextday;   
end


    end
    
end

%Remove final days for history
history.Infection = history.Infection(:,:,1:(Weeks*7));
history.Isolation = history.Isolation(:,:,1:(Weeks*7));
history.Symptomatic = Symptomatic;
history.gettingatest = history.gettingatest(:,:, 1:(Weeks*7));
history.IsolatingthroughCovid = history.IsolatingthroughCovid(:,:,(1:Weeks*7));
history.pos_test_day = history.pos_test_day(:,1:(Weeks*7));
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
   
   
  %  Now.PCR  =   (Symptomatic.*Prob_profiles.PCR_symday(history.Infection(:,:,Day) +1) + (1-Symptomatic).*Prob_profiles.PCR_asymday(history.Infection(:,:,Day) + 1));
  %  Now.lat  =   (Symptomatic.*Prob_profiles.lat_symday(history.Infection(:,:,Day) +1) + (1-Symptomatic).*Prob_profiles.lat_asymday(history.Infection(:,:,Day) + 1));
    
    
end


function Now = Update4(Now,histInf, Prob_profiles, Symptomatic, Day, K, histIsol)

Now.Infect = K*((1-histIsol).*(Prob_profiles.newvar(Day)*(Symptomatic.*Prob_profiles.Inf_symday(histInf +1) + (1-Symptomatic).*Prob_profiles.Inf_asymday(histInf +1))));

end

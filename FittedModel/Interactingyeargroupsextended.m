function [history] = Interactingyeargroups(school_pop, Infection, Testing, Strategy, Adherence, Prob_profiles, randnum, school_pop2)
%This function simulates the spread of infection between pupils within a
%secondary school formed of close contact groups in within year groups.

%This version of the model code explicitly tracks R_school, i.e. it must 
%explicitly track who infects who. Interactingyeargroupsmuchquicker is used
%for model fitting. More quantities are also tracked in this version 

%% --- Inputs --- %%
if nargin == 0
% - school_pop and school_pop2 are structs containing the close contact
    [school_pop] = SchoolPopulation(200, 5, 25, 0);
    school_pop2 = school_pop;
%   structure of the whole school (.school_matrix) , a vector containing
%   the year group each row/column of the matrix belong to (.yeargroup_vec)
%   the number of year groups (.YearGroups) and the size of each
%   year group (.YearSize). 

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
    %masstestinguptake (0-1) %proportion of individuals agreed to mass
                             %testing - assumed to be 1 in paper but lower uptake 
                             %is now captured through the Adherence vector 
                             %probtakelatflow
    Strategy.masstestinguptake = 1;
    %SCTuptake (0-1) %%Proportion of individuals agreed to participate 
                     %in serial contact testing - assumed to be 1 in paper
    Strategy.SCTuptake = 1;
    %initialuptake (0-1)  %uptake to initial testing - no initial testing
                       %assumed
    Strategy.initialuptake = 1;
    
%- Adherence is a struct containing parameters relevant to adherence

    %used in paper
    %probtakelatflow % vector of daily values of the expected proportion of
                     % the school population who take an LFT each dsay
    Adherence.probtakelatflow = 0.5*ones(1,Infection.Weeks*7);
    
    %not used in paper
    % takenproperly (0-1) % of LFTs taken properly (scaling sensitivity)
    Adherence.takenproperly = 0;
    % serialnonadherence (0-1)  %proportion of individuals who agree to be 
                                %tested but never actually take home tests
    Adherence.serialnonadherence = 0;
    % C_isolate (0-1) %relative force of external infection on isolating individual
    Adherence.C_isolate = 0;
    
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
  
  
  % Set random seed
  randnum = randi(1e6);
  rng(randnum);
end
            
%% ----Output ----
%   - history is a structure, containing the infection, testing, and
%   isolation details throughout the course of the simulation. Useful
%   outputs from this can be obtained using the files Modeloutputs.m and
%   Moremodeloutputs.m

%Authors: Trystan Leng and Edward M. Hill 
%Last update 11/11/2021.

%{
if nargin == 0 %To update
    %random number   
    randnum = 5;
    rng(randnum);
    
    %School population
    school_pop = SchoolPopulation(200, 5, 200, 0);
    
    
    %Infection_parameters
    %Infection.K = 3*rand+1; %secondary infections expected from symptomatic individual (without depletion of susceptibles)

    
    
    Infection.Ext = 0.0013133; %probability of external infection to each non-individual on daya not isolating
    Infection.Inf_0 = 0.02; %Initial population prevalence
    Infection.Rec_0 = 0.2;  %Initial population immunity
    Infection.K_asym = 0.4*rand + 0.3; % Relative infectiousness of asymptomatic individuals
    Infection.Sym_proportion = 0.12 + 0.19*rand; %Proportion of school population who are symptomatic
    Infection.alpha_withinyear = 0; %Interaction with other members of the year group (not close contacts)
    Infection.alpha_betweenyears = 0; %Interaction between year groups
    Infection.leak_infect = 0; %0 or 1, whether infection occurs on test days
    
    
    Infection.Weeks = 6; %Number of weeks (including week before term)
    
    Infection.HolidayWeek = []; % Weeks of school holidays
    
    %Testing_parameters
    Testing.sens_PCR = 1; % parameter deciding which column of the PCR sensitivity profile to read
    %(1) - baseline sensitivity, (2) - low sensitivity, (3) - high sensitivity 
    Testing.sens_lat = 1; %parameter deciding which column of the LFT sensitivity profile to read.
    %(1) - baseline sensitivity, (2) - low sensitivity, (3) - high sensitivity 
    Testing.spec_PCR = 0.997 + 0.003*rand; %specificity of PCR tests
    Testing.spec_lat = 0.997 + 0.003*rand; %specificity of LFT tests
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
   % Adherence.probtakelatflow = 0.05*ones(1, Infection.Weeks*7); %probability of taking lat flow (0 - 1), larger integers could specify the inclusion of correlations into this
   Adherence.probtakelatflow = [1,1,1,1,1,1,1,0.0160368054295554,0.242141672776696,0.221278350684790,0.173843395227479,0.230068323684321,0.222998568751878,0.00797233690200095,0.0482577242359299,0.242081620407951,0.189148536967751,0.192661343905288,0.175212024639949,0.0847528832194434,0.0148324218119442,0.272193520554720,0.125323134353830,0.0576525837018385,0.223569836157121,0.115413466974662,0.0363096125621343,0.0146869103030617,0.282932629197924,0.0183326536808124,0.0252409857929954,0.134351263755906,0.0681450670187002,0.0251627124063831,0.0164964370211048,0.0984050450148963,0.0698763203501310,0.0290835674905219,0.0752266271176480,0.0752266271176480,0.0752266271176480,0.0752266271176480];
   
   
    Adherence.serialnonadherence = 0; %proportion of individuals who agree to be tested but never actually take home tests
    
    
    Adherence.takenproperly = 1; %proportion of LFTs that are not taken properly
    
    
    
   % Adherence.isolatingproperly = [0,0,0,0,0,0,0,0.867924528301887,0.844243792325056,0.802469135802469,0.750000000000000,0.701863354037267,0,0,0.601503759398496,0.577519379844961,0.552123552123552,0.496078431372549,0.442386831275720,0,0,0.393762183235867,0.378947368421053,0.393442622950820,0.354838709677419,0.356250000000000,0,0];
    
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
%}

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

Prob_profiles.lat_symday = Adherence.takenproperly*Prob_profiles.lat_symday;
Prob_profiles.lat_asymday = Adherence.takenproperly*Prob_profiles.lat_asymday;


YearSize =  school_pop.yearsize;
YearGroup = school_pop.yeargroups;
Ks = Infection.K*ones(YearGroup,YearSize);
Ks2 = Infection.K2*ones(YearGroup, YearSize);
Ks3 = Infection.K3*ones(YearGroup, YearSize);

%%Infectivity profile%%
Prob_profiles.Inf_symday = [0;Prob_profiles.Infectivity; zeros(max(51,Weeks*7-15 + 50),1)];
Prob_profiles.Inf_asymday = Infection.K_asym*Prob_profiles.Inf_symday;

%% ---- storing history ----
history.Infection = zeros(YearGroup, YearSize, Weeks*7+10); %Infection status of individuals through time
history.Isolation = zeros(YearGroup, YearSize, Weeks*7+10); %Whether individuals are isolated or not through time
history.Rday = zeros(2, Weeks*7+10); % First row -> number of individuals 
                                     % infected on day d
                                     % Second row -> number of secondary
                                     % infections from individuals infected
                                     % on day d
                                     
history.DayInf = zeros(YearGroup,YearSize); %Day individual is infected
history.TotInfstudent = zeros(YearGroup, YearSize, Weeks*7+10); %total infection to a student, stored through time
history.ext_or_int = zeros(YearGroup, YearSize,Weeks*7+10); %Infected externally or from pupil-to-pupil tansmission
history.Isolated_Infecteds = zeros(YearGroup, Weeks*7 + 10); %Number of isolated infecteds each day
history.pos_test_day = zeros(YearGroup, Weeks*7+10); %Number of positive tests that day

history.taken_test = zeros(YearGroup, YearSize, Weeks*7+10); %storing whether an individual has take an LFT test on a day 
history.posLFT = zeros(YearGroup, YearSize, Weeks*7+10); %storing whether an individual tests +ve to an LFT on a day
history.posLFTandPCR = zeros(YearGroup, YearSize, Weeks*7+10); %storing whether an inividuasl tests +ve to an LFT and PCR on a day

history.InfectedTwice = zeros(YearGroup, Weeks*7+10); %storing whether an individual is infected through pupil-to-pupil
                                                      %transmission and community transmission on the same time-step

history.IsolatingthroughCovid = zeros(YearGroup, YearSize,Weeks*7 +10); %whether an individual is isolating because they are confirmed covid positive
history.gettingatest = zeros(YearGroup, YearSize, Weeks*7+10); %2 for in-school test, 1 for home test


%% Things that depend on randomness
rng(randnum);
%generate random numbers beforehand
%for PCR tests and for infection       
r_lat = rand(YearGroup, YearSize, Weeks*7);
r_PCR = rand(YearGroup, YearSize, Weeks*7);
%for confirmatory PCR tests
r_conf = rand(YearGroup, YearSize, Weeks*7);


%for testing day
r_testingday = randi(5, YearGroup, YearSize) - 1;
r_taketest = rand(YearGroup, YearSize, Weeks*7);

%Set who is symptomatic
Symptomatic = zeros(YearGroup, YearSize);
%who is symptomatic?
for yr = 1:YearGroup
    Symptomatic(yr, randsample(YearSize, binornd(YearSize, Infection.Sym_proportion))) = 1;
end

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





%Find groups for mass testin, SCT and initial testing
%assume minimum agree to all
%N.B. In paper, assumed all pupils agree to mass testing and SCT, so not
%used
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
        
    % tests occur on random days
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
%prop_tested = squeeze(sum(sum(history.gettingatest > 0)));


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


%N.B Assume same close contact structure in every year group, so just use
%submatrix from first year group
YearGroupMatrix1 = full(subMatrix(school_pop.school_matrix, yr, YearSize));
YearGroupMatrix2 = full(subMatrix(school_pop2.school_matrix, yr, YearSize));

%Chosen isolation, mass testing, and SCT strategies
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
                  
       Ext = Infection.Ext*Infection.commext(day);
       Ext_isolate = Adherence.C_isolate*Ext;
                    
     %% --- Testing students via PCR (self-seeking) --- %%
        for yr = 1:YearGroup
                       
            for student = all_students(will_take_PCR(yr,:) == day)
                    if Now.had_PCR(yr, student) == 0 %student does not seek PCR if already tested positive to a PCR
                        %student seeks test
                        if r_PCR(yr, student, day) < Now.PCR(yr, student)
                            %student tests positive
                            Isolation_period = 10; %Isolates for 10 days from symptom onset                            
                            notisolmarker = 0; %not already isolating 
                            Now.had_PCR(yr, student) = 1; %has had a +ve PCR so set this to 1
                           
                            if history.Isolation(yr,student,day) == 0
                                notisolmarker = 1; %update this to 1 if already isolating
                            end
                            
                            history.Isolation(yr, student, day:(day+Isolation_period)) = 1; %isolates for 10 days
                            history.IsolatingthroughCovid(yr,student,(day+PCR_delay):(day + Isolation_period)) = 1;  %isolating because confirmed covid positive
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
        
         %% --- mass testing using LFTs -----
            
              for yr = 1:YearGroup
                 
               % find group who are mass testing -> those who actually take
              % a test that day but have not already taken a PCR
                 masstestinggroup =  all_students(history.gettingatest(yr,:,day) > 0 & will_take_PCR(yr,:) ~= day);
                 
                 for student = masstestinggroup
                   
                   %To match Testing uptake
                   if history.gettingatest(yr,student,day) == 1
                         P_taketest = min(Adherence.probtakelatflow(day)/prop_tested(day), 1); %tested at home (via LFT)
                      else
                         P_taketest = 1; %SCT testing
                   end                    
                   %For % participation plots
                   %{
                   if history.gettingatest(yr,student,day) == 1
                         P_taketest =  Adherence.probtakelatflow(day);
                      else
                         P_taketest = 1; %SCT testing
                   end    
                   %}
                      
                      if r_taketest(yr,student,day) < P_taketest %Find out if pupil takes test
                          
                         history.taken_test(yr,student,day) = 1; %pupil takes test

                            if r_lat(yr, student, day) < (1 - history.Isolation(yr,student,day))*Now.lat(yr,student) %Find out if pupil tests +ve (isolating pupils assumed not to take test so not test +ve)
                                                               
                             history.posLFT(yr,student,day) = 1;    %record student test positive to an LFT that day 
                      
                               
                                  history.ever_test_pos(yr,student) = 1; %record if student ever has tested positive

                                  %confirmatory PCR
                                 if history.Infection(yr, student, day) == 0 || history.Infection(yr, student, day) > 15
                                      %false positive, but recovered, so may
                                      %test positive by chance                                 
                                     if  r_conf(yr, student, day) < Now.PCR(yr, student)
                                          %tests positive - so close
                                          %contacts will isolate and those
                                          %agreeing to be tested will test
                                          Isolation_period = 10;
                                          Testing_period = 7;                                 
                                          %Assume they will not seek another PCR test
                                          Now.had_PCR(yr,student) = 1; 

                                           
                                     history.posLFTandPCR(yr,student,day) = 1; %record student test positive to both tests that day
                                     history.IsolatingthroughCovid(yr,student,(day+PCR_delay):(day + Isolation_period)) = 1;  %isolating because confirmed covid positive

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
                                      history.posLFTandPCR(yr,student,day) = 1; %record student test positive to both tests that day
                                      history.IsolatingthroughCovid(yr,student,(day+PCR_delay):(day + Isolation_period)) = 1; %isolating because confirmed covid positive

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
         
      
      %Not isolating during holiday weeks
      if any(Week == Infection.HolidayWeek)
            history.Isolation(:,:,day) = 0;
      end
      
      %% ---- Transmission ----
            
      %tracking R version
      
       Fromotheryears = zeros(1,YearGroup); %storing the number of individuals infected in each year from other years 
      
      %% ---- infections within school --- 

      if ~(Week == 1 ||day_of_week == 6 || day_of_week == 7 || any(Week == Infection.HolidayWeek))
         %if pupils are in school
          
          for yr = 1:YearGroup
            
              %find infectious pupils attending school
              j = find(Now.Infect(yr,:) > 0 & history.Isolation(yr,:,day) == 0);
                           
              for student = j
                 
                  %close contact infections            
                  closecontactprob = Now.Infect(yr,student)/(NumCloseContacts + alpha_withinyear*(YearSize-NumCloseContacts - 1) + alpha_betweenyears*(YearSize*(YearGroup-1)));
                  contactprob = closecontactprob;
                  
                  %find group who are close contacts of student are not yet infected, and are not
                  %isolating
                  for contact = all_students((YearGroupMatrix(student,:) == 1 & history.Infection(yr,:,day) == 0 & history.Isolation(yr,:,day) == 0 ) & history.Infection(yr,:,day+1) == 0)
                            if rand < contactprob
                                
                                %contact  is infected
                                history.Infection(yr,contact,day+1) = 1; %update contact infection status
                                history.DayInf(yr,contact) = day; %set day of infection of contact to be day
                                history.Rday(1,day) = history.Rday(1,day)+1; %add 1 to the number of pupils infected that day                              
                                history.ext_or_int(yr, contact, day+1) = 2; %2 representing infected from pupil-to-pupil transmission

                                if Symptomatic(yr, contact)
                                    %set day will seek and take PCR
                                    will_take_PCR(yr, contact) = day + day_of_PCR(yr, contact);                                                  
                                end
                                
                                history.Rday(2, history.DayInf(yr,student)) = history.Rday(2, history.DayInf(yr,student)) + 1; %add 1 to the number of secondary infections
                                                                                                                               %from indivs infected on history.DayInf(yr,student)
                            end                      
                  end
                                 
                  
                  %within year infections
                  withinyearprob = alpha_withinyear*closecontactprob;
                  contactprob = withinyearprob;
                  
                  %find group who are not close contacts of student, are in
                  %the same year, are not yet infected, and are not
                  %isolating
                  for contact  = all_students(((YearGroupMatrix(student,:) == 0 & history.Infection(yr,:,day) == 0 & history.Isolation(yr,:,day) == 0 )) & history.Infection(yr,:,day+1) == 0)
                            if rand < contactprob
                                
                                %contact  is infected
                                history.Infection(yr,contact,day+1) = 1; %update contact infection status
                                history.DayInf(yr,contact) = day; %set day of infection of contact to be day
                                history.Rday(1,day) = history.Rday(1,day)+1; %add 1 to the number of pupils infected that day                                 
                                history.ext_or_int(yr, contact, day+1) = 2; %2 representing infected from pupil-to-pupil transmission

                                if Symptomatic(yr, contact)
                                    %set day will seek and take PCR
                                    will_take_PCR(yr, contact) = day + day_of_PCR(yr, contact);                                                  
                                end

                                history.Rday(2, history.DayInf(yr,student)) = history.Rday(2, history.DayInf(yr,student)) + 1; %add 1 to the number of secondary infections
                                                                                                                               %from indivs infected on history.DayInf(yr,student)
                            end                      
                  end
                                                  
          %between year infections
             
                 if alpha_betweenyears > 0   %only do this step if there is infection between years
  
                      closecontactprob = Now.Infect(yr,student)/(NumCloseContacts + alpha_withinyear*(YearSize-NumCloseContacts - 1) + alpha_betweenyears*(YearSize*(YearGroup-1)));  

                      %rest of school infections
                      restofschoolprob = alpha_betweenyears*closecontactprob;

                      other_years = 1:YearGroup; other_years(yr) = []; %only loop through other years

                     for whichyear = other_years

                         %find non-isolating susceptible population of
                         %whichyear
                         suspop = length(all_students((history.Infection(whichyear,:,day) == 0 & history.Infection(whichyear,:, day+1) == 0) & history.Isolation(whichyear,:,day) == 0));

                         %find out how many are infected by student (not who
                         %yet)
                         binominfect = mybinornd(suspop, restofschoolprob);                                      
                         Fromotheryears(whichyear) = Fromotheryears(whichyear) + binominfect; %add num infected by student to the num infected Fromotheryears

                         history.Rday(2, history.DayInf(yr,student)) = history.Rday(2, history.DayInf(yr,student)) + binominfect; % add num infected by student to the no.

                     end
                 
                 end  
                  
             end
               
          end
          
          %now loop through and find who is infected from other years 
          for whichyear = 1:YearGroup
                    
                     %find non-isolating susceptible population of
                     %whichyear
                     susgroup = all_students((history.Infection(whichyear,:,day) == 0 & history.Infection(whichyear,:, day+1) == 0) & history.Isolation(whichyear,:,day) == 0);                      
                     actualinfecteds = randsample(susgroup, min(Fromotheryears(whichyear), length(susgroup))); %choose who is infected
                     
                    for contact = actualinfecteds

                             history.Infection(whichyear,contact,day+1) = 1; %update contact infection status
                             history.DayInf(whichyear,contact) = day; %set day of infection of contact to be day
                             history.Rday(1,day) = history.Rday(1,day) + 1; %add 1 to the number of pupils infected that day 
                             history.ext_or_int(whichyear,contact,day+1) = 2; %2 representing infected from pupil-to-pupil transmission

                            if Symptomatic(whichyear,contact)
                                %set day will seek and take PCR
                                will_take_PCR(whichyear,contact) = day + day_of_PCR(whichyear,contact);                                                  
                            end                                         
                    end         
          end
                      
      end
                
        %% ------ infections from external sources ----- 
      
      if any(Week == Infection.HolidayWeek)
          %higher force of external infection during holidays
          TodayExt  = Infection.HolidayExt*Ext;
      else
          TodayExt = Ext;
      end
            
      for yr = 1:YearGroup
          
        
           susgroup = all_students(history.Infection(yr,:,day) == 0); %find susceptible population of yr
           %update infections
           
           %mexed version      
           [histInfnextday, history.DayInf(yr,:), history.Rday, histExtorIntnextday, will_take_PCR(yr,:)] = Interactingyeargroups_externalinfection_mex(susgroup,  history.DayInf(yr,:),  Symptomatic(yr,:), will_take_PCR(yr,:), day_of_PCR(yr,:), history.Isolation(yr,:,day), history.Rday, rand(1,length(susgroup)), TodayExt, Ext_isolate, day);                   
           %non-mexed version (equivalent but slower)
           %[histInfnextday, history.DayInf(yr,:), history.Rday, histExtorIntnextday, will_take_PCR(yr,:)] = Interactingyeargroups_externalinfection(susgroup,  history.DayInf(yr,:),  Symptomatic(yr,:), will_take_PCR(yr,:), day_of_PCR(yr,:), history.Isolation(yr,:,day), history.Rday, rand(1,length(susgroup)), TodayExt, Ext_isolate, day);         

           
           %see if individuals have been recorded as infected twice (from
           %pupil-to pupil transmission and external transmission)
           history.InfectedTwice(yr,day+1) = history.InfectedTwice(yr,day+1) + sum((histInfnextday == 1) & history.Infection(yr,:,day+1) == 1);
           
           
           %if infected twice need to adjust things
           if history.InfectedTwice(yr,day+1) > 0
           
               j = find(histInfnextday >0 & history.Infection(yr,:,day+1) >0); %find individuals infected twice 
               otherstudents = all_students;
               otherstudents(j) = []; %individuals not infected twice 

               tempiore = randi(2,1,length(j)); %temp internally (2) or externally (1) infected - randomly assigned because infected both           

               infgroup = find(history.Infection(:,:,day) > 0 & history.Infection(:,:,day) < 16 & history.Isolation(:,:,day) == 0); %find infected individidlas
               infgroupinfectiousness = Now.Infect(infgroup); %find inffectiousness of idnidviuals

               if sum(tempiore == 1) > 0

                   %have to find out who was assigned to having infected l, to
                   %adjust Rday(2,:)

                  for l = j(tempiore ==1) 

                     if length(infgroup) > 1                 
                         %figure out who infected l

                         %scale infectivity for rest of year non-close contacts

                         sameyear = infgroup((mod(infgroup, YearGroup) == yr)); %find out who infecteds in same year
                         notclosecontact = YearGroup*(find(YearGroupMatrix(l,:) == 0) - 1) + yr; %find out who is not close contact
                         [~, sameyearnotclose] = ismembertol(intersect(sameyear, notclosecontact), infgroup); %find intersection and their position in infgroup
                         infgroupinfectiousness(sameyearnotclose) = alpha_withinyear*infgroupinfectiousness(sameyearnotclose); %adjust infectiousness to individual
                         diffyearpos = find((mod(infgroup, YearGroup) ~= yr));                  %find out infecteds in different years and position in infgroup
                         infgroupinfectiousness(diffyearpos) = alpha_betweenyears*infgroupinfectiousness(diffyearpos);  %adjust infectiousness to individual

                         infector = randsample(infgroup, 1, true, infgroupinfectiousness);
                         history.Rday(2,history.DayInf(infector)) = history.Rday(2,history.DayInf(infector)) - 1; %adjust Rday(2,;)
                     else
                        infector  = infgroup;
                        history.Rday(2,history.DayInf(infector)) = history.Rday(2,history.DayInf(infector)) - 1; %adjust Rday(2,:)

                     end
                  end
               end

               history.ext_or_int(yr,j,day+1) = tempiore; %assign as either internally or externally infected
               history.ext_or_int(yr,otherstudents,day+1) = max(histExtorIntnextday(otherstudents), history.ext_or_int(yr,otherstudents,day+1)); %assign other students 
               history.Infection(yr,:,day+1) = 1*(histInfnextday|history.Infection(yr,:,day+1)); %add external infections to history.Infection  
               history.Rday(1,day) = history.Rday(1,day) -history.InfectedTwice(yr,day+1); %Infected twice individuals have been counted twice, so remove one instance of them from Rday          
           
           else
 
              history.Infection(yr,:,day+1) = 1*(histInfnextday|history.Infection(yr,:,day+1));   %add external infections to history.Infection  
              history.ext_or_int(yr,:,day+1) = max(histExtorIntnextday, history.ext_or_int(yr,:,day+1)); %assign all students                
          end
           
           
      end
      
           
      for yr = 1:YearGroup
          %update infection status of infected/recovered pupils
          infrecstudents = history.Infection(yr,:,day) > 0;
          history.Infection(yr, infrecstudents, day+1) = history.Infection(yr, infrecstudents, day) +1;
      end
                  
        %Updates  
        if ~isempty(Infection.HolidayWeek) && Week < Infection.HolidayWeek(1)            
            Now = Update3(Now, history.Infection(:,:,(day+1)), Prob_profiles, Symptomatic, day+1, Ks);            
        elseif isempty(Infection.HolidayWeek)
             Now = Update3(Now, history.Infection(:,:,(day+1)), Prob_profiles, Symptomatic, day+1, Ks);
        elseif ~isempty(Infection.HolidayWeek) && Week < Infection.HolidayWeek(end-1)
             Now = Update3(Now, history.Infection(:,:,(day+1)), Prob_profiles, Symptomatic, day+1, Ks2);
        else
             Now = Update3(Now, history.Infection(:,:,(day+1)), Prob_profiles, Symptomatic, day+1, Ks3);
        end
                
        %Quantities through time
        for yr = 1:YearGroup
            %Isolated infecteds       
            history.Isolated_Infecteds(yr,day) = sum(history.Isolation(yr, history.Infection(yr,:,day) > 0 & history.Infection(yr,:,day) < 16, day));
            history.Tot_Isolated(yr,day) = sum(history.Isolation(yr,:,day));
        end    
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
    Now.PCR  =   (Symptomatic.*Prob_profiles.PCR_symday(history.Infection(:,:,Day) +1) + (1-Symptomatic).*Prob_profiles.PCR_asymday(history.Infection(:,:,Day) + 1));
    Now.lat  =   (Symptomatic.*Prob_profiles.lat_symday(history.Infection(:,:,Day) +1) + (1-Symptomatic).*Prob_profiles.lat_asymday(history.Infection(:,:,Day) + 1));
    
    
end

function Now = Update3(Now, histInf, Prob_profiles, Symptomatic, Day, Ks)
    %Symptomatics need to be transposed when 1 year group 
    %{
    Now.Infect = (Symptomatic'.*Probs.Inf_sym(history.Infection(:,:,Day) +1) + (1-Symptomatic)'.*Probs.Inf_asym(history.Infection(:,:,Day) +1))';
    Now.PCR  =   (Symptomatic'.*Probs.PCR_sym(history.Infection(:,:,Day) +1) + (1-Symptomatic)'.*Probs.PCR_asym(history.Infection(:,:,Day) + 1))';
    Now.lat  =   (Symptomatic'.*Probs.lat_sym(history.Infection(:,:,Day) +1) + (1-Symptomatic)'.*Probs.lat_asym(history.Infection(:,:,Day) + 1))';
    %}
    %But not with larger year groups? figure out
    Now.Infect = Ks.*(Prob_profiles.newvar(Day)*(Symptomatic.*Prob_profiles.Inf_symday(histInf +1) + (1-Symptomatic).*Prob_profiles.Inf_asymday(histInf +1)));
    Now.PCR  =   (Symptomatic.*Prob_profiles.PCR_symday(histInf +1) + (1-Symptomatic).*Prob_profiles.PCR_asymday(histInf + 1));
    Now.lat  =   (Symptomatic.*Prob_profiles.lat_symday(histInf +1) + (1-Symptomatic).*Prob_profiles.lat_asymday(histInf + 1));
    
    
end

function ret = mybinornd( n, p )
  ret = sum(rand(1,n)<p);
end

function subM = subMatrix(M, whichyear, YearSize)
            subM = M(((whichyear-1)*YearSize + 1): ((whichyear-1)*YearSize + YearSize), ((whichyear-1)*YearSize + 1): ((whichyear-1)*YearSize + YearSize));
end

function [history] = Interactingyeargroups(params, PCR_test_sym, PCR_test_asym, lat_test_sym, lat_test_asym, randnum)
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
%Last update 05/02/2021.

if nargin == 0
    PCR_test_sym = readtable('PCR_Curve_summary.csv');
    PCR_test_sym = table2array(PCR_test_sym(:, 2:4));
    PCR_test_asym = csvread('PCR_Curve_asym.csv');

    lat_test_sym = readtable('lat_Curve_summary.csv');
    lat_test_sym = table2array(lat_test_sym(:, 2:4));
    lat_test_asym = csvread('lat_Curve_asym.csv');
    
    randnum = 5;
    rng(randnum);
    params = [3*rand+1, 0.00135, 1, 1, 1, 0.02, 0.2, 0, 0.4*rand + 0.3, 0.12 + 0.19*rand, 0, 200, 5, 0, 0];
    
    
end

%params are K, Ext, isol_or_test, sens_lat, sens_PCR, Inf_0, Rec_0, C_isolate, K_asym, Sym_proportion, alpha
K = params(1); % 'school R0' of symptomatic indviduals
Ext = params(2); % constant force of infection on non-isolating individuals
isol_or_test = params(3); % proportion of population not agreeing to be tested
sens_lat = params(4); %parameter deciding which column of the LFT sensitivity profile to read.
%This parameter currently also contains the options for no testing, regular mass testing and serial contact testing
%The details of each value are included below:
%(1) - baseline sensitivity, (2) - low sensitivity, (3) - high sensitivity %These assume serial contact testing
%(4) - No school testing
%(5)- baseline sensitivity, (9) - low sensitivity (10) - high sensitivity %These assume weekly mass testing
%(7) - baseline sensitivity, (11) - low sensitivity, (12) - high sensitivity %These assume combined serial contact testing and weekly mass testing 
%(6) - baseline sensitivity, biweekly mass testing
%(8) - baseline sensitivity, biweekly mass testing and serial contact

%(13) - weekly mass testing and isolation of close contacts (no serial
%contact testing)
%testing

%(14) - twice weekly mass testing and isolation of close contacts (no serial
%contact testing)
%testing (15) - low sensitivity, (16) - high sensitivity

sens_PCR = params(5); % parameter deciding which column of the PCR sensitivity profile to read
%(1) - baseline sensitivity, (2) - low sensitivity, (3) - high sensitivity 

Inf_0 = params(6); %initial population prevelance
Rec_0 = params(7); %initial population immunity
C_isolate = params(8); %external force of infection on isolating individuals 
K_asym = params(9); % Relative infectiousness of asymptomatic individuals
Sym_proportion = params(10); %Proportion of school population who are symptomatic
alpha = params(11); %Interaction between year groups

YearSize = params(12);  %Size of year
YearGroup = params(13); %How many year groups (or how many close contact groups)
leak_infect = params(14); %whether infection happens on day of tests
freqordens = params(15); % whether frequency (0) or density dependent

if freqordens
    const =  200-1;   %hardcoded to be same as below for YearSize = 200, alpha = 0. (YearGroup = 5).   
else
    const =  (YearSize + ((YearGroup-1)*alpha*YearSize - 1));
end

    

Not_testing_prop = isol_or_test; %Individuals not getting tested
Not_initial_prop = isol_or_test; %Individuals initially not getting tested 

%1 = baseline, 2 = lower bound, 3 = upper bound
%Probability of PCR tests
Probs.PCR_sym = [0;  PCR_test_sym(10*(1:30) + 1, sens_PCR); zeros(69,1)];
%Probs.PCR_asym = Probs.PCR_sym;
Probs.PCR_asym = [0; PCR_test_asym(10*(1:30) + 1, sens_PCR); zeros(69,1)];


%Probability of lateral flow tests
if sens_lat == 4
    Probs.lat_sym = zeros(100,1);
   Probs.lat_asym = zeros(100,1);
   
elseif sens_lat == 9 || sens_lat == 11 || sens_lat == 15
    %low sensitivity weekly mass testing
    Probs.lat_sym = [0; lat_test_sym(10*(1:30) + 1, 2); zeros(69,1)];
    %Probs.lat_asym = Probs.lat_sym;
    Probs.lat_asym = [0; lat_test_asym(10*(1:30) + 1, 2); zeros(69,1)];
    
elseif sens_lat == 10 || sens_lat == 12 || sens_lat == 16
    %high sensitivity weekly mass testing
    Probs.lat_sym = [0; lat_test_sym(10*(1:30) + 1, 3); zeros(69,1)];
    %Probs.lat_asym = Probs.lat_sym;
    Probs.lat_asym = [0; lat_test_asym(10*(1:30) + 1, 3); zeros(69,1)];
   
elseif sens_lat > 4 
    %baseline params in this instance
    Probs.lat_sym = [0; lat_test_sym(10*(1:30) + 1, 1); zeros(69,1)];
    Probs.lat_asym = [0; lat_test_asym(10*(1:30) + 1, 1); zeros(69,1)];
    
else
    Probs.lat_sym = [0; lat_test_sym(10*(1:30) + 1, sens_lat); zeros(69,1)];
    Probs.lat_asym = [0; lat_test_asym(10*(1:30) + 1, sens_lat); zeros(69,1)];
    
end


%specificity here for now
spec = 0.9997;

if sens_lat ~= 4
%don't do this step for sens_lat == 4, as mass testing doesn't happen
    Probs.lat_sym(Probs.lat_sym < (1- spec)) = (1-spec);
    Probs.lat_asym(Probs.lat_asym < (1-spec)) = (1-spec);
end


%Infectivity profile since day of infection
Infectivity_since_infection = [0.0063 0.0563 0.1320 0.1798 0.1817 0.1521 0.1117 0.0746 0.0464 0.0272 0.0152 0.0082 0.0043 0.0022 0.0011 0.0005 0.0002 0.0001 0.0001 0.0000];
Infectivity_since_infection = Infectivity_since_infection'/sum(Infectivity_since_infection);
Infectivity_since_infection(15) = sum(Infectivity_since_infection(15:end)); Infectivity_since_infection = Infectivity_since_infection(1:15);

%Symptom onset (i.e. incubation period)
Symptom_onset =  [0.0055 0.0534 0.1307 0.1814 0.1847 0.1545 0.1129 0.0747 0.0458 0.0265 0.0146 0.0077 0.0039 0.0020 0.0010 0.0005 0.0002 0.0001 0 0];
Symptom_onset = Symptom_onset/sum(Symptom_onset);
Symptom_onset(15) = sum(Symptom_onset(15:end)); Symptom_onset = Symptom_onset(1:15);


%%Parameters%%
Weeks = 8; %Weeks
PCR_delay = 2;    %PCR test delay

Ext_isolate = C_isolate*Ext;

%%Infectivity profile%%
Probs.Inf_sym = K*[0;Infectivity_since_infection; zeros(85,1)];
Probs.Inf_asym = K_asym*Probs.Inf_sym;

%Keep entire history
history.Infection = zeros(YearGroup, YearSize, Weeks*7+10); %Infection status of individuals through time
history.Isolation = zeros(YearGroup, YearSize, Weeks*7+10); %Whether individuals are isolated or not through time
history.TotInf = zeros(YearGroup, Weeks*7+10); %Total number of infected individuals each day
history.ext_or_int = zeros(YearGroup, YearSize,Weeks*7+10); %Infected externally or from school
history.Isolated_Infecteds = zeros(YearGroup, Weeks*7 + 10); %Number of isolated infecteds each day

history.Known_Infecteds = zeros(YearGroup, Weeks*7+10); %Number of known infecteds each day
history.Known_asymptomatics = zeros(YearGroup, Weeks*7+10); %Number of known infected asymptomatics each day
history.pos_test_day = zeros(YearGroup, Weeks*7+10); %Number of positive tests that day
history.ever_test_pos = zeros(YearGroup, YearSize); %logical if ever test positive
history.presyms = zeros(YearGroup, YearSize); %Number Presymptomatic individuals captured


history.day_caught = zeros(YearGroup, YearSize, (Weeks-1)*7); %Day of mass testing that individuals are caught
Rollingtestdays = zeros(1,YearGroup); %which 'day' of mass testing is it?



%not fixed yet!
history.Rday = zeros(2, Weeks*7+10); % First row -> number of individuals 
                                     % infected on day d
                                     % Second row -> number of secondary
                                     % infections from individuals infected
                                     % on day d
                                     
history.DayInf = zeros(YearGroup,YearSize); %Day individual is infected



%Test days for those getting tested
history.Test_days = zeros(YearGroup, Weeks*7);
history.Test_days(:,1) = 1; %test individuals on the monday before term
history.Test_days(:,4) = 1; %test individuals on the Thursday before term





%Individual's infection probability and testing probabilities
Now.Infect = zeros(YearGroup, YearSize); 
Now.PCR = zeros(YearGroup, YearSize);
Now.lat = zeros(YearGroup, YearSize);
Now.had_PCR = zeros(YearGroup, YearSize); %have they had a +ve PCR?


Symptomatic = zeros(YearGroup, YearSize); %Proportion who are symptomatic




if sens_lat == 6 || sens_lat > 7 
   %3 years tested monday and thursday, 2 years tested tuesday and friday
    for day = 8:Weeks*7
      for yr = 1:YearGroup
          if yr < 1+(YearGroup/2)
              if mod(day, 7) == 1 || mod(day,7) == 4
                  history.Test_days(yr,day) = 1;
                  history.Test_days(yr,day) = 1;
              end
          else
              if mod(day, 7) == 2 || mod(day,7) == 5
                  history.Test_days(yr,day) = 1;
                  history.Test_days(yr,day) = 1;
              end
          end
       end
    end
   
elseif sens_lat == 5 || sens_lat == 7 || sens_lat > 8    
    %different year tested every day
    for day = 8:Weeks*7
        for dow = 1:5
            if mod(day, 7) == dow
                history.Test_days(dow, day) = 1;
            end
        end
        
    end

    
end








%% Things that depend on randomness
rng(randnum);
%generate random numbers beforehand
%for PCR tests, for LF tests, and for infection       
r_lat = rand(YearGroup, YearSize, Weeks*7);
%r_inf = rand(YearGroup, YearSize, Weeks*7);
%r_inf = rand(YearGroup,YearSize, YearGroup, YearSize, Weeks*7); %Internal infections
r_infext = rand(YearGroup, YearSize, Weeks*7);
%r_inf = r_infext;
r_PCR = rand(YearGroup, YearSize, Weeks*7);

%for confirmatory PCR tests
r_conf = rand(YearGroup, YearSize, Weeks*7);



%who is symptomatic?
for j = 1:YearGroup
    Symptomatic(j, randsample(YearSize, round(Sym_proportion*YearSize))) = 1;
end

%day an individual will seek PCR test if symptomatic i.e. during of
%presymptomatic period    
day_of_PCR = zeros(YearGroup, YearSize);
PCR_days = randsample(1:15, sum(sum(Symptomatic)), true, Symptom_onset); 
day_of_PCR(Symptomatic == 1) = PCR_days;

will_take_PCR = zeros(YearGroup, YearSize);

%who is initially infected and recovered?
for j = 1:YearGroup
    init = mnrnd(YearSize, [Inf_0; Rec_0; 1 - Inf_0 - Rec_0]);
    init_not_sus = randsample(YearSize, init(1)+init(2));
    init_Inf = init_not_sus(1:init(1));   
    init_Rec = init_not_sus((init(1)+1):(init(1)+init(2)));
               
    %randomly assign which day of infection infecteds are at
    history.Infection(j, init_Inf,1) = randi(15, 1, init(1));
    %assume recovereds do not test positive to any tests
    history.Infection(j, init_Rec,1) = 32;
   
    % testing for initially infecteds
    for student = init_Inf'
        history.DayInf(j,student) = 1;
        history.Rday(1,1) = history.Rday(1,1)+1;

        if Symptomatic(j, student)
             will_take_PCR(j,student) = day_of_PCR(j, student) - history.Infection(j, student,1 ) +1; 
             
             if will_take_PCR(j, student) < 1 
                 if rand <  Probs.PCR_asym(day_of_PCR(j, student))
                     %tests positive before simulation begins
                     history.Isolation(j, student, 1: (1 + 9 + will_take_PCR(j,student))) = 1; %if tests is on day 0, then student should isolate  from day 1 to day 10                                                                                            %if test is on day -1, then student should isolate from day 1 to day 8 
                     history.ever_test_pos(j, student) = 1;
                 end
             end

            
        end        
    end
end



%Initial test and infection probabilities 
Now = Update2(Now, history, Probs, Symptomatic, 1);


%set population of who gets tested
getting_tested_1 = 1:YearSize;
%all_students = [getting_tested_1; getting_tested_1; getting_tested_1; getting_tested_1; getting_tested_1];
all_students = repmat(getting_tested_1, YearGroup, 1);
Not_getting_tested = zeros(YearGroup, round(Not_testing_prop*YearSize));
getting_tested = zeros(YearGroup, YearSize - round(Not_testing_prop*YearSize));
Not_initially_tested = zeros(YearGroup, round(Not_initial_prop*YearSize));
initially_tested = zeros(YearGroup, YearSize - round(Not_initial_prop*YearSize));

for j = 1:YearGroup
    
    %Set getting tested throughout the term
    Not_getting_tested(j,:) = randsample(YearSize, round(Not_testing_prop*YearSize));
    getting_tested_temp = getting_tested_1;
    getting_tested_temp(Not_getting_tested(j,:)) = [];
    getting_tested(j,:) = getting_tested_temp;

    %set initially tested
    Not_initially_tested(j,:) = randsample(YearSize, round(Not_initial_prop*YearSize));
    getting_tested_temp = getting_tested_1;
    getting_tested_temp(Not_initially_tested(j,:)) = [];
    initially_tested(j,:) = getting_tested_temp;   
    
end
    
rng(randnum);

%Simulation
for Week = 1:Weeks
    for day_of_week = 1:7
            day = (Week-1)*7 + day_of_week;
            
        %Update which 'test day' it is if there is SCT    
        for j = 1:YearGroup
           if history.Test_days(j, day) == 1
               Rollingtestdays(j) = Rollingtestdays(j) + 1;
           else
               Rollingtestdays(j) = 0;
           end
        end
        
       
        
            
        %Test students via PCR (self seeking) with probability       
        for j = 1:YearGroup
            for student = all_students(j, will_take_PCR(j,:) == day)
                    if Now.had_PCR(j, student) == 0 
                        %student seeks test
                        Now.had_PCR(j, student) = 1;
                        if r_PCR(j, student, day) < Now.PCR(j, student)
                            %student tests positive
                            notisolmarker = 0;
                            
                            

                                
                            
                            
                            if history.Isolation(j,student,day) == 0
                                notisolmarker = 1;
                            end
                            
                            history.Isolation(j, student, day:(day+10)) = 1; %isolates for 10 days from day following symptom onset
                            history.ever_test_pos(j, student) = 2;
                            
                            
                            history.pos_test_day(j, day) = history.pos_test_day(j, day) +1;                      

                            if Week > 1 && notisolmarker
                                %school only isolates and tests  after children are back in
                                %school
                                if sens_lat ~= 5 && sens_lat ~= 6 && sens_lat ~= 9 && sens_lat ~=10 && sens_lat ~=4
                                    
                                 if sens_lat < 13
                                    
                                    nottestingandnotisol = intersect(Not_getting_tested(j,:), find(history.Isolation(j,:,day) == 0)); 
                                     
                                    history.Isolation(j, nottestingandnotisol, (day+PCR_delay+1):(day+10)) = 1; %isolates indivs not getting tested                    
                                    history.Test_days(j, (day+PCR_delay+1):(day+7)) = 1; %school tests 
                                else
                                    %sens_lat = 13 is isolation + mass
                                    %testing
                                    notisol = find(history.Isolation(j,:,day) == 0); 
                                    history.Isolation(j, notisol, (day+PCR_delay+1):(day+10)) = 1; %isolates indivs  getting tested                    
                                end
                                    
                                    
                                end
                            end
                                                      
                            
                        else
                            %student tests negative
                            history.Isolation(j, student, day:(day+PCR_delay)) = 1; %isolates until PCR returns
                        end
                    end
            end
        end
        
        
         %Mass testing students by lateral flow
         if day_of_week ~= 6 && day_of_week~=7
             for j = 1:YearGroup
                 if history.Test_days(j, day) == 1
                     
                     if day == 1 || day == 4
                         testing_group = initially_tested(j,:);
                     else   
                         testing_group = getting_tested(j,:);
                     end
                     
                      for student = testing_group %only test students agreed to be tested
                          if r_lat(j, student, day) < (1 - history.Isolation(j,student,day))*Now.lat(j,student)
                              %lateral flow tests positive
                              history.ever_test_pos(j,student) = 1;
                              
                              if history.Infection(j, student, day) == 0
                                  %false positive so confirmatory PCR comes
                                  %back negative
                                  Isolation_period = PCR_delay;
                                  Testing_period = PCR_delay;
                                  
                              elseif history.Infection(j, student, day) > 15
                                  %false positive, but recovered, so may
                                  %test positive by chance                                 
                                 if  r_conf(j, student, day) < Now.PCR(j, student)
                                      %tests positive
                                      Isolation_period = 10;
                                      Testing_period = 7;                                 
                                      %Assume they will not seek another PCR test
                                      Now.had_PCR(j,student) = 1;                                    
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
                                  Now.had_PCR(j,student) = 1;                                
                              end

                              history.pos_test_day(j,day) = history.pos_test_day(j,day) + 1;

                              history.Isolation(j, student, (day):(day+Isolation_period)) = 1; %isolates +ve individual
                             if Week > 1 && sens_lat ~= 5 && sens_lat ~= 6 && sens_lat ~= 9 && sens_lat ~= 10
                                %year only isolates and tests  after children are back in
                                %school   
                                %senslat 5,6,9 and 10 represents regular mass tests,
                                %so positive tests don't trigger isolation
                                %or test days
                                
                                if sens_lat < 13
                                    nottestingandnotisol = intersect(Not_getting_tested(j,:), find(history.Isolation(j,:,day) == 0)); 
                                    history.Isolation(j, nottestingandnotisol, (day+1):(day+Isolation_period)) = 1; %isolates indivs not getting tested                    
                                    history.Test_days(j, (day+1):(day+Testing_period)) = 1; %school tests 
                                else
                                    %sens_lat = 13 is isolation + mass
                                    %testing
                                    notisol = find(history.Isolation(j,:,day) == 0); 
                                    history.Isolation(j, notisol, (day+1):(day+Isolation_period)) = 1; %isolates indivs  getting tested                    
                                end
                                
                                
                                
                                
                                
                              end
                              
                              if Week > 1
                                  %only count this from when the school
                                  %starts
                                 history.day_caught(j, student, Rollingtestdays(j)) = 1; 
                              end
                                                           
                             if history.Infection(j, student, day) <  day_of_PCR(j,student)
                                 %count individuals caught in
                                 %presymptomatic stage
                                history.presyms(student,day) = 1;
                             end
                             
                                                                 %do not isolate individual on that day if
                            %leak_infect
                           if leak_infect
                               history.Isolation(j, student, day) = 0;  %does not isolate on that day
                            end
                             
                            
                          end                         
                          
                      end
                 end
             end            
         end
         
         
      %Isolating at weekends if serial contact testing
       %if sens_lat equals 4, then no isolation at weekends as noone is being tested    
       if sens_lat ~= 4         
           for j = 1:YearGroup
                if (day_of_week == 6 || day_of_week == 7) && history.Test_days(j, day) == 1
                    %if this is a test day isolate testing students because weekend
                    history.Isolation(j, getting_tested, day) = 1;        
                end 
           end
       end
       
       
       
       %Infection 
       %{
       
       if Week == 1 ||day_of_week == 6 || day_of_week == 7
            %weekend so no school infection or before kids are back to
            %school
            history.TotInf(:,day) = Ext;             
       else

       
       %Infectivity of individuals who are attending school
       Infnow = Now.Infect.*(1-history.Isolation(:,:,day));     
       yg_inf = prod(1 - Infnow'/const)'; %probability not infected by other members of year
       rs_inf = ((prod(prod(1 - (alpha/const)*Infnow')))./prod(1 - (alpha/const)*Infnow'))'; %probability not infected by other years
       rs_inf(isnan(rs_inf)) = 0; %set isnans to 0       
       history.TotInf(:,day) = 1 - (yg_inf.*rs_inf).*(1-Ext);
       
       end
       
        for j = 1:YearGroup
            for student = 1:YearSize
                if history.Infection(j,student,day) == 0  
                %only susceptible individuals can be infected
                    if r_inf(j,student, day) < (1 - history.Isolation(j,student,day))*history.TotInf(j,day) + history.Isolation(j,student,day)*Ext_isolate

                       %find out if internally or externally infected
                       if r_inf(j, student,day) < (1 - history.Isolation(j,student,day))*Ext + history.Isolation(j,student,day)*Ext_isolate
                           %1 for external
                           history.ext_or_int(j, student, day) = 1;
                       else
                           %2 for school
                           history.ext_or_int(j, student, day) = 2;

                       end

                       history.Infection(j, student,day+1) = 1;

                       if Symptomatic(j, student)
                          %set day will seek and take PCR
                          will_take_PCR(j, student) = day + day_of_PCR(j, student);
                       end                   

                    end

                else
                %update infection status of all other individuals
                history.Infection(j, student, day+1) = history.Infection(j, student, day) +1;
                end
            end

        end
       
       %}
       
       
    
       %Within-school infections
       if ~(Week == 1 ||day_of_week == 6 || day_of_week == 7)
         %if pupils are in school
         
         %because there are no other random things going on I can probably
         %do this just every day?
         
         r_inf = rand(YearGroup,YearSize, YearGroup, YearSize); %Internal infections
       % r_inftemp1(day) = sum(sum(sum(sum(r_inf))));
       % r_inftemp2(day) = r_inf(1,20,3,50);
         
         if sum(sum(Now.Infect > 0 & history.Isolation(:,:,day) == 0)) > 0
         

             
          for yr = 1:YearGroup
            
              %find infectious pupils attending school
              j = find(Now.Infect(yr,:) > 0 & history.Isolation(yr,:,day) == 0);
                           
              for student = j
                 
                  %year group infections           
                  contactprob = Now.Infect(yr,student)/const;
                  
                  %find group who are close contacts of student are not yet infected, and are not
                  %isolating
                  for contact = all_students(yr,(history.Infection(yr,:,day) == 0 & history.Isolation(yr,:,day) == 0 ) & history.Infection(yr,:,day+1) == 0)
                            if r_inf(yr,student, yr, contact) < contactprob
                            %if rand < contactprob
                                
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
                  
                  
                  if alpha > 0
                      
                      contactprob = alpha*Now.Infect(yr,student)/const;

                      
                      otheryears = 1:YearGroup;
                      otheryears(yr) = [];
                      
                      for otheryear = otheryears
                          for contact = all_students(otheryear,(history.Infection(yr,:,day) == 0 & history.Isolation(otheryear,:,day) == 0 ) & history.Infection(otheryear,:,day+1) == 0)
                              if r_inf(yr,student, otheryear, contact) < contactprob
                              %if rand < contactprob
                                
                                %contact  is infected
                                history.Infection(otheryear,contact,day+1) = 1; %update contact infection status
                                history.DayInf(otheryear,contact) = day; %set day of infection of contact to be day
                                
                                
                                
                                history.Rday(1,day) = history.Rday(1,day)+1; %add 1 to the number of pupils infected that day                              
                                history.ext_or_int(otheryear, contact, day+1) = 2; %2 representing infected from pupil-to-pupil transmission

                                if Symptomatic(otheryear, contact)
                                    %set day will seek and take PCR
                                    will_take_PCR(otheryear, contact) = day + day_of_PCR(otheryear, contact);                                                  
                                end
                                
                                history.Rday(2, history.DayInf(yr,student)) = history.Rday(2, history.DayInf(yr,student)) + 1; %add 1 to the number of secondary infections
                                                                                                                               %from indivs infected on history.DayInf(yr,student)
                              end    
                              
                          end
                      
                      end
                  
                  end
              end
          end
          
         end
        
       end
       
       
       
       
       %External infections
        for yr = 1:YearGroup
          
        
           susgroup = all_students(yr,history.Infection(yr,:,day) == 0 & history.Infection(yr,:,day+1) == 0); %find susceptible population of yr
           %update infections
           
           %mexed version      
           [histInfnextday, history.DayInf(yr,:), history.Rday, histExtorIntnextday, will_take_PCR(yr,:)] = Interactingyeargroups_externalinfection(susgroup,  history.DayInf(yr,:),  Symptomatic(yr,:), will_take_PCR(yr,:), day_of_PCR(yr,:), history.Isolation(yr,:,day), history.Rday, r_infext(yr,susgroup,day), Ext, Ext_isolate, day);
           
           history.Infection(yr,:,day+1) = 1*(histInfnextday|history.Infection(yr,:,day+1));   %add external infections to history.Infection  
           history.ext_or_int(yr,:,day+1) = max(histExtorIntnextday, history.ext_or_int(yr,:,day+1)); %assign all students
        end
        
      for yr = 1:YearGroup
          %update infection status of infected/recovered pupils
          infrecstudents = history.Infection(yr,:,day) > 0;
          history.Infection(yr, infrecstudents, day+1) = history.Infection(yr, infrecstudents, day) +1;
      end
      %}  
        %Updates
        Now = Update2(Now, history, Probs, Symptomatic, day+1);
        
        %Quantities through time
        for j = 1:YearGroup
            %Isolated infecteds       
            history.Isolated_Infecteds(j,day) = sum(history.Isolation(j, history.Infection(j,:,day) > 0 & history.Infection(j,:,day) < 16, day));
            history.Known_Infecteds(j,day) = sum(history.ever_test_pos(j, history.Infection(j,:,day) > 0 & history.Infection(j,:,day) < 16) ~= 0);
            history.Known_asymptomatics(j,day) = sum(history.ever_test_pos(j, history.Infection(j,:,day) > 0 & history.Infection(j,:,day) < 16 & (1-Symptomatic(j,:))) ~= 0);
            history.Tot_Isolated(j,day) = sum(history.Isolation(j,:,day));
        end   
            
    end
    
end


history.Infection = history.Infection(:,:,1:(Weeks*7));
history.Isolation = history.Isolation(:,:,1:(Weeks*7));
history.TotInf = history.TotInf(:,1:(Weeks*7));
history.ext_or_int = history.ext_or_int(:,:, 1:(Weeks*7));
history.Isolated_Infecteds = history.Isolated_Infecteds(:, 1:(Weeks*7));
history.Known_Infecteds = history.Known_Infecteds(:, 1:(Weeks*7));
history.Known_asymptomatics = history.Known_asymptomatics(:, 1:(Weeks*7));
history.pos_test_day = history.pos_test_day(:,1:(Weeks*7));
history.Symptomatic = Symptomatic;

%history.r_inf1 = r_inftemp1;
%history.r_inf2 = r_inftemp2;


end
    
function Now = Update2(Now, history, Probs, Symptomatic, Day)
    %Symptomatics need to be transposed when 1 year group 
    %{
    Now.Infect = (Symptomatic'.*Probs.Inf_sym(history.Infection(:,:,Day) +1) + (1-Symptomatic)'.*Probs.Inf_asym(history.Infection(:,:,Day) +1))';
    Now.PCR  =   (Symptomatic'.*Probs.PCR_sym(history.Infection(:,:,Day) +1) + (1-Symptomatic)'.*Probs.PCR_asym(history.Infection(:,:,Day) + 1))';
    Now.lat  =   (Symptomatic'.*Probs.lat_sym(history.Infection(:,:,Day) +1) + (1-Symptomatic)'.*Probs.lat_asym(history.Infection(:,:,Day) + 1))';
    %}
    %But not with larger year groups? figure out
    Now.Infect = (Symptomatic.*Probs.Inf_sym(history.Infection(:,:,Day) +1) + (1-Symptomatic).*Probs.Inf_asym(history.Infection(:,:,Day) +1));
    Now.PCR  =   (Symptomatic.*Probs.PCR_sym(history.Infection(:,:,Day) +1) + (1-Symptomatic).*Probs.PCR_asym(history.Infection(:,:,Day) + 1));
    Now.lat  =   (Symptomatic.*Probs.lat_sym(history.Infection(:,:,Day) +1) + (1-Symptomatic).*Probs.lat_asym(history.Infection(:,:,Day) + 1));
end
    
    


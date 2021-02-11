function [Prev_true, Prev_asym, Known_asym, pos_tests, Absences, Tot_Infected, Schooldays_missed, Infected_within_school, Infected_during_term, AsymsCaptured, AsymsTotal, Peak_infect, Tot_Isolating, Days_tested, Presyms_Captured, Day_caught_vec] =   Modeloutputs(history)
%This function extracts useful quantities from the history struct outputted
%by interactingyeargroups.m

%Input:
%   - history is the structure outputted from running interactingyeargroups.m

%Output:
%   - Outputs are described next to their introduction

%Authors: Trystan Leng and Edward M. Hill 
%Last update 05/02/2021.


Weeks = 8;
YearSize = 200;


for day = 1:(Weeks*7)
    %True prevalence through time
    Prev_true(:, day) = sum(history.Infection(:,:,day) > 0 & history.Infection(:,:,day) < 16, 2);
    
    %prevalence of asymptomatics through time
    Prev_asym(:,day) = sum(history.Infection(:,:,day) > 0 & history.Infection(:,:,day) < 16 & (1-history.Symptomatic), 2);
end

%Known asymptomatics (those that have tested positive)
Known_asym = history.Known_asymptomatics;

%Positive tests on a given day
pos_tests = history.pos_test_day;

%pupils absent from school 
Schooldays = [8:12, 15:19, 22:26, 29:33, 36:40, 43:47, 50:54];
Weekends = 1:56;
Weekends(Schooldays) = [];
Absences = history.Tot_Isolated;
Absences(:, Weekends) = 0;

%total infected over the course of the simulation
Tot_Infected = sum(history.Infection(:,:,56) > 0, 2) - sum(history.Infection(:, :,1) == 32, 2);

%schooldays missed for each student
Schooldays_missed = sum(history.Isolation(:,:, Schooldays), 3);

%Infected within school
Infected_within_school = sum(sum(history.ext_or_int(:,:, 8:end) == 2, 3), 2);

%Infected during the half term
Infected_during_term = sum(sum(history.ext_or_int(:,:, 8:end) ~= 0, 3), 2);

%Total infected asymptomatics
AsymsCaptured = sum(((1-history.Symptomatic).*(history.Infection(:,:,Weeks*7)  > 0) & history.Infection(:,:,Weeks*7) < 87).*(history.ever_test_pos ~= 0), 2); 

%Total  asymptomatics
AsymsTotal = sum((1-history.Symptomatic).*(history.Infection(:,:,Weeks*7)  > 0) & history.Infection(:,:,Weeks*7) < 87, 2);

%Peak infected
Peak_infect = max(sum(Prev_true));

%Total days isolating
Tot_Isolating = sum(history.Isolation, 3);

%how many school years are isolating each day?
Days_tested = zeros(1, 56);
Days_tested(Schooldays) = sum(history.Test_days(:, Schooldays));

%Symptomatics captured before developing symptoms
Presyms_Captured=  sum(sum(history.presyms));

%'Day of serial contact testing' caught
Day_caught_vec = sum(sum(history.day_caught),2);
Day_caught_vec = Day_caught_vec(:);



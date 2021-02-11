function [Incidence, Num_lat, Num_PCR, Num_conf, School_incidence] = Moremodeloutputs(history)
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
YearGroup = 5;

for day = 1:(Weeks*7)
    %True prevalence through time
    Incidence(:, day) = sum(history.Infection(:,:,day) == 1, 2);
end

%Number of tests
Schooldaysandinitial = [1, 4, 8:12, 15:19, 22:26, 29:33, 36:40, 43:47, 50:54];
Not_isolating = sum( (1 - history.Isolation(:,:, Schooldaysandinitial)), 2);
Not_isolating = reshape(Not_isolating, YearGroup, length(Schooldaysandinitial));
%Number of LFTs N.B. this only applies for strategies involving mass
%testing and for all individuals opting into testing
Num_lat = sum(sum(history.Test_days(:, Schooldaysandinitial).*Not_isolating));
%Number of PCR tests
Num_PCR = sum(sum(history.ever_test_pos == 2));
%Number of confirmatory PCR tests (the same as the number of positive LFTs)
Num_conf = sum(sum(history.ever_test_pos == 1));


%Within school incidence on last day of term
School_incidence = sum(sum(history.ext_or_int(:,:, 54) == 2));

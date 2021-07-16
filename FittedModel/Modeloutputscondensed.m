function [pos_PCRsday, peakcovabs, Absences, Prev_true, num_LFTs, pos_LFTs, pos_LFTsandPCR, Incidence] = Modeloutputscondensed(history, Weeks)

%Useful outcome measures obtained from Interacting year groups



pos_PCRsday = sum(history.pos_test_day(:,8:(Weeks*7))); %number of positive PCR tests each day for each year
num_LFTs = squeeze(sum(sum(history.taken_test(:,:,8:(Weeks*7))))); %number of LFTs taken each day for each year 
pos_LFTs = squeeze(sum(sum(history.posLFT(:,:,8:(Weeks*7))))); %number of positive LFTs taken each day for each year
pos_LFTsandPCR = squeeze(sum(sum(history.posLFTandPCR(:,:,8:(Weeks*7))))); %number of positive LFTs that result in positive confirmatory PCR each day for each year



%pupils absent from school 
Schooldays = 8:(Weeks*7);
Weekends = zeros(1, length(8:Weeks*7));
Weekends(mod((8:Weeks*7), 7) == 6) = 1;
Weekends(mod((8:Weeks*7), 7) == 0) = 1;
Weekends = Schooldays(Weekends == 1);
Absences = history.Tot_Isolated; 
Absences(:, Weekends) = 0; %Don't count isolations as absences at weekends


%pupils absent because they have tested positive to COVID
CovAbsences = history.IsolatingthroughCovid;
CovAbsences(:, Weekends) = 0; %Don't count isolations as absences at weekends


if ~isempty(history.HolidayWeek)
   %Don't count isolations as absences at weekends
   Absences(:, (((history.HolidayWeek-1)*7) + 1): (((history.HolidayWeek-1)*7) + 7)) = 0; 
   CovAbsences(:, (((history.HolidayWeek-1)*7) + 1): (((history.HolidayWeek-1)*7) + 7)) = 0; 

end


peakcovabs = squeeze(sum(sum(CovAbsences))); %records total number of absences in the school each day - can use max over the relevant time-frame to find peak confirmed covid absences



Prev_true = squeeze(sum(history.Infection > 0 & history.Infection < 16, 2)); %Model prevalence
Incidence = squeeze(sum(history.Infection == 1, 2));  %Model incidence


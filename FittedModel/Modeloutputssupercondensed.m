function [pos_PCRsday, peakcovabs, pos_LFTs] = Modeloutputssupercondensed(history, Weeks)


pos_PCRsday = sum(history.pos_test_day(:,8:(Weeks*7)));

%pupils absent from school 
Schooldays = 8:(Weeks*7);

Weekends = zeros(1, length(8:Weeks*7));

Weekends(mod((8:Weeks*7), 7) == 6) = 1;
Weekends(mod((8:Weeks*7), 7) == 0) = 1;

Weekends = Schooldays(Weekends == 1);



%Schooldays = [8:12, 15:19, 22:26, 29:33, 36:40, 43:47, 50:54];
%Weekends = 1:56;
%Weekends(Schooldays) = [];


%{
Absences = history.Tot_Isolated;
Absences(:, Weekends) = 0;
%}

CovAbsences = history.IsolatingthroughCovid;

%CovAbsences = history.IsolatingthroughCovid_day;
CovAbsences(:, Weekends) = 0;


if ~isempty(history.HolidayWeek)
   %Absences(:, (((history.HolidayWeek-1)*7) + 1): (((history.HolidayWeek-1)*7) + 7)) = 0; 
   CovAbsences(:, (((history.HolidayWeek-1)*7) + 1): (((history.HolidayWeek-1)*7) + 7)) = 0; 

end

%peakcovabs = squeeze(sum(sum(history.IsolatingthroughCovid)));

peakcovabs = squeeze(sum(sum(CovAbsences)));
%peakcovabs = sum(CovAbsences);


pos_LFTs = squeeze(sum(sum(history.posLFT(:,:,8:(Weeks*7)))));
%pos_LFTs = sum(history.posLFT_day(:,(8:Weeks*7)));

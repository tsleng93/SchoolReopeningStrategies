function [histInfnextday, histDayInf, histRday, histExtorIntnextday, will_take_PCR] = Interactingyeargroups_externalinfection(SusIndiv,  histDayInf,  Symptomatic, will_take_PCR, day_of_PCR, histIsol, histRday, r_infext, TodayExt, Ext_isolate, day)
      

histInfnextday = zeros(size(histIsol));
histExtorIntnextday = zeros(size(histIsol));



for ii = 1:length(SusIndiv)
    
    
    student = SusIndiv(ii);
    
     if histIsol(student) == 0
         if r_infext(ii) < TodayExt
             histInfnextday(student) = 1;
             histDayInf(student) = day;
             histRday(1,day) = histRday(1,day)+1;

             histExtorIntnextday(student) = 1;


             if Symptomatic(student)
                  %set day will seek and take PCR
                  will_take_PCR(student) = day + day_of_PCR(student);
             end   



         end
     elseif Ext_isolate > 0
        if r_infext(ii) < Ext_isolate
        % if rand  < Ext_isolate
            histInfnextday(student) = 1;
            histDayInf(student) = day;
            histRday(1,day) = histRday(1,day)+1;

            histExtorIntnextday(student) = 1;

           if Symptomatic(student)
              %set day will seek and take PCR
              will_take_PCR(student) = day + day_of_PCR(student);
           end                          
         end
     end
end
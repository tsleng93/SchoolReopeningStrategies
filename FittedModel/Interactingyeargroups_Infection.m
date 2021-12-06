function [histInfNextday, will_take_PCR] = Interactingyeargroups_Infection(histInf, histTotInf,   will_take_PCR, Symptomatic, day,  day_of_PCR)

%Try seeing if making this 2D speeds things up

histInfNextday = zeros(size(histInf));

YearGroup = size(histInf,1);
YearSize = size(histInf,2);

r_inf = rand(YearGroup,YearSize);

    for yr = 1:YearGroup
            for student = 1:YearSize
                

                if histInf(yr,student) == 0  
                %only susceptible individuals can be infected                
                    if r_inf(yr,student) < histTotInf(yr,student)
                                             

                       histInfNextday(yr, student) = 1;

                       if Symptomatic(yr, student)
                          %set day will seek and take PCR
                          will_take_PCR(yr, student) = day + day_of_PCR(yr, student);
                       end                   

                    end

                else
                %update infection status of all other individuals
                histInfNextday(yr, student) = histInf(yr, student) +1;
                end
            end

    end

% This function implements the CUSUM test. The theory and pseudocode is
% obatined from 'Adaptive Filtering and Change Detection',
% F.Gustafsson, Pg 66

function [totAlarms, avgRes, alarms, percentGIsZero] = testCUSUM(residual)

    g = 0;
  
    noOfTimeInstants = size(residual,1);
    alarms = zeros(noOfTimeInstants,1);
    
    gSeries = zeros(noOfTimeInstants,1);
    percentGIsZero = 0;
    
    totAlarms = 0;   
    
    threshold = 500;   
    forgetFact = 2.5;
    
    avgRes = mean(residual);

    for i = 1:noOfTimeInstants 
    
        if(residual(i,1) > 0)
            g = g + residual(i,1) - forgetFact;  

            if(g < 0)
               g = 0;
            end

            if(g > threshold)
               g = 0;
               alarms(i,1) = 1;    
               totAlarms = totAlarms + 1;
            else
               alarms(i,1) = 0;
            end

            gSeries(i,1) = g;
        end
        

    end
    
    for j = 1:noOfTimeInstants
       if(gSeries(j,1) == 0)
           percentGIsZero = percentGIsZero + 1;
       end
    end
    
    percentGIsZero = percentGIsZero / noOfTimeInstants;
   
end
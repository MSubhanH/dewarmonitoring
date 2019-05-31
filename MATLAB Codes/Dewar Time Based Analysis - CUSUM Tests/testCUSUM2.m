% This function implements the LS CUSUM test. The theory and pseudocode is
% obatined from 'Adaptive Filtering and Change Detection',
% F.Gustafsson, Pg 68

function [totAlarms, avgRes, alarms, percentGIsZero] = testCUSUM2(residual)

    g1 = 0;
    g2 = 0;
     
    noOfTimeInstants = size(residual,1);
    alarms = zeros(noOfTimeInstants,1);
    
    g1Series = zeros(noOfTimeInstants,1);
    g2Series = zeros(noOfTimeInstants,1);
    percentGIsZero = 0;
    
    totAlarms = 0;   
    
    threshold = 25;   %18;
    forgetFact = 0.12;
    
    avgRes = mean(residual);

    t = 0:15:size(residual,1)*15;
    
    for i = 2:noOfTimeInstants 
           
        if(residual(i,1) > 0)
            
            res = 1/(t(i) - t(1))*sum(residual(2:i, 1));
            
            s1 = res;
            s2 = -res;
            
            g1 = max([g1 + s1 - forgetFact, 0]);  
            g2 = max([g2 + s2 - forgetFact, 0]);
            
            if(g1 > threshold || g2 > threshold)
               g1 = 0;
               g2 = 0;
               alarms(i,1) = 1;    
               totAlarms = totAlarms + 1;               
            else    
               alarms(i,1) = 0;
            end

            g1Series(i,1) = g1;
            g2Series(i,1) = g2;
        end
        

    end
    
    for j = 1:noOfTimeInstants
       if(g1Series(j,1) == 0)
           percentGIsZero = percentGIsZero + 1;
       end
    end
    
    percentGIsZero = percentGIsZero / noOfTimeInstants;
   
end
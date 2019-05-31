function [f1, f2, f3,f4,f5, y0] = normFeatures(dateAligned, sigDataAligned, sigETAligned, dayOfYearAligned)

    startRangeData = 1;
    endRangeData =   size(sigDataAligned,1); 

    f1 = sigDataAligned(startRangeData:endRangeData,4);
    f1 = (f1(:,1) - min(f1(:,1))) / (max(f1(:,1)) - min(f1(:,1)));
    f1 = log(f1) ;

    for i = 1:size(f1,1)
       if(f1(i,1) == -Inf)
          f1(i,1) = f1(i-1,1); 
       end
    end

    f1 = (f1(:,1) - min(f1(:,1))) / (max(f1(:,1)) - min(f1(:,1)));

    
    
    f2 = sigETAligned(startRangeData:endRangeData,2);
    f2 = (f2(:,1) - min(f2(:,1))) / (max(f2(:,1)) - min(f2(:,1)));

    
    
    for i = 1:size(dayOfYearAligned(startRangeData:endRangeData,1),1)
        if(dayOfYearAligned(i,1) > 305)
            cycleDay(i,1) = dayOfYearAligned(i,1) - 305 + 1;
        else
            cycleDay(i,1) = 365 - ((305 - dayOfYearAligned(i,1)) - 1);
        end
    end


    f3 = (cycleDay(:,1) - 1)  / (365 - 1);

    sigma = 0.1;
    coldIndex = (1/(sigma*sqrt(2*pi))) * exp(-0.5*((f3(:,1) - 0.2527)/sigma).^2);
    coldIndex = ((coldIndex(:,1) - 0) / (8 - 0))*(1 - 0.8) + 0.8;

    f3 = coldIndex;
    
    
    f4 = 25 - sigETAligned(startRangeData:endRangeData,2);
    f4 = (f4(:,1) - min(f4(:,1))) / (max(f4(:,1)) - min(f4(:,1)));

    x = 1:size(sigETAligned,1);
    coeff = polyfit(x', sigETAligned(:,2), 3);
    eqn_et = coeff(1,1)*(x.^3) + coeff(1,2)*(x.^2) + coeff(1,3)*x + coeff(1,4);

    f5 = sigETAligned(:,2) - eqn_et';
    f5 = (f5(:,1) - min(f5(:,1))) / (max(f5(:,1)) - min(f5(:,1)));

    
    y0 = sigDataAligned(startRangeData:endRangeData,3);
    
    for i = 1:size(y0,1)
        if(y0(i,1) == 0)
           y0(i,1) = y0(i-1,1); 
        end
    end

end
function [totAlarms, avgRes, alarms, testQuantityRes] = testQuantile(y0, eqn)

    C = [2.30046776476642,2.04421545918589,1.90359863595407,1.76796015098872,1.61777484665221,1.46956740874188,1.30567087032978,1.15977351692190,0.953300304551104,0.598386737013665,0.0732726375841887;
        -1.13787358751092,-1.66375010567253,-1.94719643813235,-2.14675064107972,-2.32067574041200,-2.47394482486541,-2.63286273846755,-2.82909459802170,-3.05419449375825,-3.26389945613897,-3.49732521768526;
        -0.0828521335859643,-0.576107230721373,-0.834787206967767,-1.02999971721861,-1.19194016344773,-1.38266549783672,-1.55536747447899,-1.77055170658599,-2.01912846760733,-2.30097842377618,-2.75019427084803;
        0.768636509980674,0.454651641983256,0.235561493882798,0.0465315518351459,-0.165677076625891,-0.391642627173204,-0.615682227945935,-0.825430457179592,-1.07125080489373,-1.31511461687061,-1.71209433354796;
        1.90022766850592,1.56895899400586,1.29051197917614,0.990982347441981,0.745794365112853,0.537457518376578,0.285824061004898,0.0685003574505601,-0.202924862060755,-0.529224790570599,-1.09581654537788];

    indicesRes = find(y0 >  eqn);
    resModelAndStage2 = y0(indicesRes) - eqn(indicesRes);  
    
    testSignalRes = resModelAndStage2;
    avgRes = mean(resModelAndStage2);
    
    tau = [0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99];
    noOfQuants = size(tau,2);
    N = 100;
    seqQuantEstTestSig = zeros(size(testSignalRes,1),size(tau,2));

    for k = 1:noOfQuants

        n = 0;
        m = 0;
        seqQuantEst = 0;

        for i = 1:size(testSignalRes,1)
            [seqQuantEst,n,m] = seqQuantile(testSignalRes(i,1),seqQuantEst,n,m,tau(1,k),N); 
            seqQuantEstTestSig(i,k) = seqQuantEst;
        end
    end

    testQuantity = zeros(size(seqQuantEstTestSig,1),1);
    K = 5;

    % Method 1 - Distance = distance from nearest distribution

    for j = 1:size(seqQuantEstTestSig,1)   
        minDistFromRef = 1000;
        for i = 1:K
            normOfDiff = norm((seqQuantEstTestSig(j,3:11) - C(i,3:11)));       

            if(normOfDiff < minDistFromRef)
                minDistFromRef = normOfDiff;
            end
        end

        testQuantity(j,1) = sqrt(minDistFromRef);
    end
    
    testQuantityRes = zeros(size(y0,1),1);
    testQuantityRes(indicesRes) = testQuantity;

    thresh = 1.5;

    alarms(1:size(y0,1))= 0;
    indicesAlarms = find(testQuantityRes > thresh);
    alarms(indicesAlarms) = 1;
    
    totAlarms = size(find(alarms == 1),2);

end
% Project: Dewar Monitoring ==========================================%
% Code : TimeBasedAnalysis_QuantileBased_ModelTesting.mat 

% Description : This code was written to test the trained reference
% distributionns on Stage 2 testing data for change detection based on
% quantile estimation.

% The algorithm works by computing the probability distribution of the test
% residual signal at each time instant and then finding its normalized distance 
% from the nearest reference distribution. If the distance exceeds beyond a
% threshold, then a change is said to be detected.

% The theory and the pseudocode for this work has been taken from the
% following paper:

% 'Residual change detection using low-complexity sequential quantile
% estimation', D. Jung, E. Frisk and M. Krysander


clc
clear all

%% Loading Trained CDFs (Reference probability distributions)

C = [2.30046776476642,2.04421545918589,1.90359863595407,1.76796015098872,1.61777484665221,1.46956740874188,1.30567087032978,1.15977351692190,0.953300304551104,0.598386737013665,0.0732726375841887;
    -1.13787358751092,-1.66375010567253,-1.94719643813235,-2.14675064107972,-2.32067574041200,-2.47394482486541,-2.63286273846755,-2.82909459802170,-3.05419449375825,-3.26389945613897,-3.49732521768526;
    -0.0828521335859643,-0.576107230721373,-0.834787206967767,-1.02999971721861,-1.19194016344773,-1.38266549783672,-1.55536747447899,-1.77055170658599,-2.01912846760733,-2.30097842377618,-2.75019427084803;
    0.768636509980674,0.454651641983256,0.235561493882798,0.0465315518351459,-0.165677076625891,-0.391642627173204,-0.615682227945935,-0.825430457179592,-1.07125080489373,-1.31511461687061,-1.71209433354796;
    1.90022766850592,1.56895899400586,1.29051197917614,0.990982347441981,0.745794365112853,0.537457518376578,0.285824061004898,0.0685003574505601,-0.202924862060755,-0.529224790570599,-1.09581654537788];

%% Concatenating Stage 2 test data sets to one vector

fileID = ['1_nm.mat';
          '2_nm.mat';
          '3_nm.mat';
          '4_bm.mat';
          '5_am.mat';
          '6_nm.mat';
          '7_nm.mat';
          '8_bm.mat';
          '9_am.mat';
          '10bm.mat';
          '11am.mat';
          '12nm.mat';  ];

dateAligned = [];      
sigDataAligned = [];
sigETAligned = [];
dayOfYearAligned = [];

for n = 1:12
    dateAlignedHandler = load(['Aligned Data Matrices\dateAligned_', fileID(n,:)]);
    sigDataAlignedHandler = load(['Aligned Data Matrices\dataAligned_', fileID(n,:)]);
    sigETAlignedHandler = load(['Aligned Data Matrices\etAligned_', fileID(n,:)]);
    dayOfYearAlignedHandler = load(['Aligned Data Matrices\dayOfYearAligned_', fileID(n,:)]);

    dateAligned = [dateAligned; dateAlignedHandler.dateAligned];
    sigDataAligned = [sigDataAligned; sigDataAlignedHandler.sigDataAligned];
    sigETAligned = [sigETAligned; sigETAlignedHandler.sigETAligned];
    dayOfYearAligned = [dayOfYearAligned; dayOfYearAlignedHandler.dayOfYearAligned];    
end

%% Computing datamodel (predicted Stage 2)
  
[f1,f2, f3,f4,f5, y0] = normFeatures(dateAligned, sigDataAligned, sigETAligned, dayOfYearAligned);
  
theta_m_sol = [5.15634525286613,0,16.4699636273839,5.52409305859874,7.36962642092408];
eqn = theta_m_sol(1,1)*f1(:,1).*(1+f4) + theta_m_sol(1,5)*f5(:,1).*(1-f4) + + (theta_m_sol(1,4)*f4) + theta_m_sol(1,3);

%% Computing residual for the Stage 2 test dataset (Stage 2 - data model)

indicesRes = find(y0 >  eqn);   % Time samples where the model overestimates, i.e when model > Stage 2, are ignored.
resModelAndStage2 = y0(indicesRes) - eqn(indicesRes);       

%% Computing quantiles estimates (probability distribution) of residual for each time step

testSignalRes = resModelAndStage2(1:end,1);

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
    
    subplot(2,1,1);
    plot(seqQuantEstTestSig(1:end, k),'-m');
    hold on;   
end

subplot(2,1,1);
plot(testSignalRes);
title('Estimated quantiles superimposed on residual plot');
xlabel('Sample no.');
ylabel('Residual');

subplot(2,1,2);
plot(fliplr(seqQuantEstTestSig(1:1000:end,2:11)), tau(2:11),'-b');
title('Blue : Estimated CDFs , Pink: Ref CDFs');
xlabel('Residual');
ylabel('Probability');
xlim([-5 8]);
hold on;

plot(fliplr(C(:,2:11)), tau(2:11),'-m','LineWidth',3);

%% Deriving test quantity - distance from ref CDFs for change detection

testQuantity = zeros(size(seqQuantEstTestSig,1),1);
K = 5;

% Method 1 - Distance = distance from nearest distribution

for j = 1:size(seqQuantEstTestSig,1)   
    minDistFromRef = 1000;
    for i = 1:K
        normOfDiff = norm((seqQuantEstTestSig(j,3:11) - C(i,3:11)));       

        if(normOfDiff < minDistFromRef)
            minDistFromRef = normOfDiff;
            refDist(j,1) = i;
        end
        
    end
    
    testQuantity(j,1) = sqrt(minDistFromRef);
end

testQuantityRes = zeros(size(y0,1),1);
testQuantityRes(indicesRes) = testQuantity;


%% Defining Alarms when distance from closest ref distribution is greater than threshold i.e (testQuantityRes > thresh)

thresh = 1.5;

alarms(1:size(y0,1))= -10;
indicesAlarms = find(testQuantityRes > thresh);
alarms(indicesAlarms) = -5;

plot(y0 - eqn); hold on;plot(testQuantityRes); hold on; plot(y0); hold on; plot([1 size(testQuantityRes,1)], [thresh thresh],'-bla'); hold on; plot(alarms,'LineWidth',2);
legend('Residual', 'Dist from closest ref CDF', 'Stage 2', 'Thresh', 'Alarms');


%% Stability Analysis

% Here we define three regions for the estimated probability distributions:
% Stable, SemiStable and Unstable

% The idea is to color-code the Stage 2 temperature in real-time for being
% in one of the three (Stable, SemiStable and Unstable) regimes for the
% ease of the operator.

indicesStable = find(testQuantity < thresh);
indicesUnstable = find(testQuantity > thresh);

plot(fliplr(seqQuantEstTestSig(indicesStable(1:1000:end),2:11)), tau(2:11),'-b');
hold on;
plot(fliplr(seqQuantEstTestSig(indicesUnstable(1:1000:end),2:11)), tau(2:11),'-r');
hold on;

estMeanQuantStable = mean(fliplr(seqQuantEstTestSig(indicesStable,2:11)));
estMeanQuantUnstable = mean(fliplr(seqQuantEstTestSig(indicesUnstable,2:11)));
estMeanQuantSemiStable = mean([estMeanQuantStable; estMeanQuantUnstable]);  % The SemiStable CDF is the mean of CDF of Stable and Unstable regime

% Plotiing the CDFs of the stability regimes

plot(estMeanQuantStable, tau(2:11),'-bla','LineWidth',5);
hold on;
plot(estMeanQuantUnstable, tau(2:11),'-bla','LineWidth',5);
hold on;
plot(estMeanQuantSemiStable, tau(2:11),'-bla','LineWidth',5);
xlim([-1 20]);

title('CDFs of Stable and Unstable regimes for test residual data. Blue: Stable, Pink: Unstable');
xlabel('Residual');
ylabel('Probability');
%% Evaluating stability of Stage 2 w.r.t defined stability distributions

% The distance between the CDF of Stage 2 and each stability regime is
% calculated for every time step. The sample is assigned to the regime
% which it is closest to (having least disatnce). 

% All of the samples are then color coded w.r.t to the assigned regime.

stabilityLevel = zeros(size(seqQuantEstTestSig,1),1);

for i = 1:size(seqQuantEstTestSig,1)
    
   normDistFromStable = norm(fliplr(seqQuantEstTestSig(i,2:11)) -  estMeanQuantStable);
   normDistFromUnstable = norm(fliplr(seqQuantEstTestSig(i,2:11)) -  estMeanQuantUnstable);
   normDistFromSemiStable = norm(fliplr(seqQuantEstTestSig(i,2:11)) -  estMeanQuantSemiStable);
   
   dist = 1000;
   
   if(normDistFromStable < dist)
       stabilityLevel(i,1) =  10;
       dist = normDistFromStable;
   end
   
   if(normDistFromUnstable < dist)
       stabilityLevel(i,1) =  1;
       dist = normDistFromUnstable ;
   end
   
   
   if(normDistFromSemiStable < dist)
       stabilityLevel(i,1) =  5;
       dist = normDistFromSemiStable;
   end
   
end

stabilitySeries = zeros(size(resModelAndStage2,1),1);
stabilitySeries(indicesRes) = stabilityLevel; 

plot(find(stabilitySeries == 0), y0(find(stabilitySeries == 0)),'-b');
hold on;
plot(find(stabilitySeries == 10), y0(find(stabilitySeries == 10)),'-g');
hold on;
plot(find(stabilitySeries == 5), y0(find(stabilitySeries == 5)),'-y');
hold on;
plot(find(stabilitySeries == 1), y0(find(stabilitySeries == 1)),'-r');
title('Stage 2 color coded with stability regimes');
xlabel('Sample no.');
ylabel('Temp (K)');
legend('Ignored', 'Stable', 'SemiStable', 'Unstable');

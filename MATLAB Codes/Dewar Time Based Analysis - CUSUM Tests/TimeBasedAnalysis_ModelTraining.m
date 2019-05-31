% Project: Dewar Monitoring 
% Code : TimeBasedAnalysis_ModelTraining

% Description : 
% The following code was written to train the Stage 2 temperature prediction data
% model on a dataset of Stage 2 that contains values of Stage 2 falling
% below the threshold of 25K.

% The features of the data model are derived from the Pressure and
% Environment Temperature (ET). Prior to this training step, the ET
% measurements were aligned with the Pressure and Stage 2 dataset. The
% aligned datasets were stored as .mat files.

clc
clear all;

%% Reading Aligned Datasets

fileID = '25K.mat'; 

dateAligned = load(['Aligned Data Matrices\dateAligned_', fileID]);  % Read mat file containing aligned timestamps
sigDataAligned = load(['Aligned Data Matrices\dataAligned_', fileID]); % Read mat file containing Stage 2 measurements falling below 25K
sigETAligned = load(['Aligned Data Matrices\etAligned_', fileID]); % Read mat file containing aligned ET measurements 

dateAligned = dateAligned.dateAligned25K;
sigDataAligned = sigDataAligned.sigDataAligned25K;
sigETAligned = sigETAligned.sigETAligned25K;

for i = 1:size(dateAligned, 1)
    dayOfYearAligned(i,1) = doy(dateAligned(i,1));  % Converting timestamp to date for all epochs
end

%% Defininf and normalizing features

% The following features were derived for the data model:

% f1 : Pressure 
% f2 : Environment Temperature
% f4 : 25K - f2
% f5 : Filtered f2 (containing mainly daily variations with seasonal
% variations filtered)

% The features are normalized for model training step.

startRangeData = 1;
endRangeData =   size(sigDataAligned,1); 


% f1 : Pressure
f1 = sigDataAligned(startRangeData:endRangeData,4);
f1 = (f1(:,1) - min(f1(:,1))) / (max(f1(:,1)) - min(f1(:,1)));
f1 = log(f1) ;
for i = 1:size(f1,1)
   if(f1(i,1) == -Inf)
      f1(i,1) = f1(i-5,1); 
   end
end
f1 = (f1(:,1) - min(f1(:,1))) / (max(f1(:,1)) - min(f1(:,1)));


% f2 : Environment Temperature
f2 = sigETAligned(startRangeData:endRangeData,2);
f2 = (f2(:,1) - min(f2(:,1))) / (max(f2(:,1)) - min(f2(:,1)));


% f4 : 25K - f2
f4 = 25 - sigETAligned(startRangeData:endRangeData,2);
f4 = (f4(:,1) - min(f4(:,1))) / (max(f4(:,1)) - min(f4(:,1)));


% f5 : Filtered f2
x = 1:size(sigETAligned,1);
coeff = polyfit(x', sigETAligned(:,2), 3);
eqn_et = coeff(1,1)*(x.^3) + coeff(1,2)*(x.^2) + coeff(1,3)*x + coeff(1,4);
f5 = sigETAligned(:,2) - eqn_et';
f5 = (f5(:,1) - min(f5(:,1))) / (max(f5(:,1)) - min(f5(:,1)));


% y0 : Training data -> Stage 2 temperature below 25K 
y0 = sigDataAligned(startRangeData:endRangeData,3);
y = (y0(:,1) - min(y0(:,1))) / (max(y0(:,1)) - min(y0(:,1)));
for i = 1:size(y0,1)
    if(y0(i,1) == 0)
       y0(i,1) = y0(i-5,1); 
    end
end

%% Least Squares Estimation Training

theta_m0 = [0,0,0,0,0];  % initializing weights
feature = [f1,f2,f4,f5];
          
fun = @(theta_m,feature) (theta_m(1,1)*f1.*(1+f4)) + (theta_m(1,5)*f5.*(1-f4)) + (theta_m(1,4)*f4) + (theta_m(1,3));

options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
[theta_m_sol, resNorm_m] = lsqcurvefit(fun,theta_m0,[f1,f2,f4,f5],y0);

eqn_m = theta_m_sol(1,1)*f1(:,1).*(1+f4) + theta_m_sol(1,5)*f5(:,1).*(1-f4) + + (theta_m_sol(1,4)*f4) + theta_m_sol(1,3);  % Final equation of data model

resStage2Eqn = (y0 - eqn_m);  % Overall residual of model and training data

plot(dateAligned,y0(:,1)); hold on; plot(dateAligned,eqn_m); hold on; plot(dateAligned,f1+theta_m_sol(1,3)); hold on;plot(dateAligned,f4+theta_m_sol(1,3)); hold on;plot(dateAligned,f5+theta_m_sol(1,3)); 
legend('Stage 2', 'Model','F1','F4', 'F5');
xlabel('Date'), ylabel('Temperature (K)');
title('Stage 2 training data vs trained data model');

preFactOfP_m = theta_m_sol(1,1)*(1+f4);
preFactOfET_m = theta_m_sol(1,5)*(1-f4);

% The data model equation 'eqn_m' was adjusted with different combinations
% of the derived features to obtain a minimal resStage2Eqn.

% The plot indicates a signficant deviation of the data model from Stage 2 
% for the period of Nov-Dec where it over-estiamtes Stage 2. This happens
% because of the seasonal variation in Stage 2 but can be adjusted for by
% omitting all over-estimations of Stage 2 during the testing stage. The
% final eqn_m is kept to be suitable as general case and to avoid
% overfitting the training dataset.
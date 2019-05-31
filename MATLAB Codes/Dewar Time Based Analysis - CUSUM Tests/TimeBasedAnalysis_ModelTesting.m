% Project: Dewar Monitoring 
% Code : TimeBasedAnalysis_ModelTesting 

% Description : 
% The following code was written to test the residuals between the trained
% data model and actual Stage 2 temperature for change detection. 

% To analyze the results, the testing data set was manually divided into
% three classes : 1) No Maintenance (nm)  2) Before Maintenance (bm) 3)
% After Maintenance (am)

% The testing datasets were assigned each class by visually examining the
% datasets for example a data set the showed a maintenance was cut in two
% new datasets : one before maintenance and one after maintenance. This
% way, 12 datasets were derived and were assigned with any one of the suffixes 
% (nm, bm or am) to indicate the assigned class. 

% The CUSUM and the LS CUSUM tests are applied on the computed
% residual for change detection. The theory and pseudocode of these tests
% are taken from the book 'Adaptive Filtering and Change Detection',
% F.Gustafsson, Pg 66,68

% The result of these tests are alarms when a potential change is detected
% in the time series. As the datasets are already classified, if an alarm
% is detected in a dataset labelled with suffix (nm or bm) then it is
% considered as a False Alarm. A rough measure of the probabilty of false
% alarm pFA for a given test is obtained by dividing total false alarms by 
% total detected alarms for all the datasets.

clc
clear all;

%% Read Aligned Matrices

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
          '12nm.mat';  ];   % Storing file names with assigned suffixes
      
totCorrectAlarms = 0;
totFalseAlarms = 0;
totAlarms = 0;
           
resTimeSeries = [];
dateTimeSeries = [];
allIndices = [];
y0TimeSeries = [];

for n = 1:12
    
    % Storing aligned testing data
    
    dateAlignedHandler = load(['Aligned Data Matrices\dateAligned_', fileID(n,:)]);
    sigDataAlignedHandler = load(['Aligned Data Matrices\dataAligned_', fileID(n,:)]);
    sigETAlignedHandler = load(['Aligned Data Matrices\etAligned_', fileID(n,:)]);
    dayOfYearAlignedHandler = load(['Aligned Data Matrices\dayOfYearAligned_', fileID(n,:)]);

    dateAligned = dateAlignedHandler.dateAligned;
    sigDataAligned = sigDataAlignedHandler.sigDataAligned;
    sigETAligned = sigETAlignedHandler.sigETAligned;
    dayOfYearAligned = dayOfYearAlignedHandler.dayOfYearAligned;
    
    [f1,f2, f3,f4,f5, y0] = normFeatures(dateAligned, sigDataAligned, sigETAligned, dayOfYearAligned);  % Normalizing features 
   
    theta_m_sol = [5.15634525286613,0,16.4699636273839,5.52409305859874,7.36962642092408];  % Assigning values of weights obtained from training stage
    eqn = theta_m_sol(1,1)*f1(:,1).*(1+f4) + theta_m_sol(1,5)*f5(:,1).*(1-f4) + + (theta_m_sol(1,4)*f4) + theta_m_sol(1,3); % Calculation of data model (Prediction of Stage 2)
  
    indicesRes = find(eqn >  22.5); 
    resModelAndStage2 = y0(indicesRes) - eqn(indicesRes);
 
    % Applying either CUSUM or LS CUSUM (uncomment the one to be used)
    
    [alarmsInSample, avgRes, alarms] = testCUSUM(resModelAndStage2);
    %[alarmsInSample, avgRes, alarms] = testCUSUM2(resModelAndStage2);
    
    
    % Plotting Alarms over Stage 2 and Model 
    % Blue : Stage 2 
    % Red : Model
    % Yellow : Alarm
    
    figure(1)
    subplot(6,2,n)
    plot(dateAligned,y0(:,1)); hold on; plot(dateAligned,eqn); hold on; plot(dateAligned(indicesRes), max(y0)*alarms); hold off;
    title(['Dataset : ', num2str(n)]);
    xlabel('Date'), ylabel('Temp(K)');
     
    % Plotting Alarms only
    figure(2)
    subplot(6,2,n);
    plot(dateAligned(indicesRes), alarms)
    xlabel('Date'), ylabel('Alarm');
    title(['Dataset : ', num2str(n)]);
    
    suffix = fileID(n,3:4);
    
    if(strcmp(suffix, 'bm'))
        totCorrectAlarms = totCorrectAlarms + alarmsInSample;
    else
        totFalseAlarms = totFalseAlarms + alarmsInSample;
    end
    
    totAlarms = totAlarms + alarmsInSample;
    
end

pFA = 0;
if(totAlarms ~= 0)
    pFA = totFalseAlarms / totAlarms;  % Rough estimate of probabilty of false alarm
end
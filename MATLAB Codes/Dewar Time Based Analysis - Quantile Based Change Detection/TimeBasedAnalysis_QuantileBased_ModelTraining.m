% Project: Dewar Monitoring ==========================================%
% Code : TimeBasedAnalysis_QuantileBased_ModelTesting.mat 

% Description : 

% The following code was written to obtain probabilty distributions of 
% nominal data which then shall be taken as reference distributions for 
% perfoming change detection using sequential quantile estimation process.

% The quantile based change detection works by computing the deviation of
% the probabilty distribution of the residual(difference between Stage 2 
% and data model)from a set of reference distributions at each time instant. 
% If the this deviation exceeds a predetermined threshold, then one can 
% conclude that a change has been detected. Therefore, the method does not
% depend on a specific statistcal quantity such as mean but rather depends
% on change in the pdf of residual.


% The theory and the pseudocode for this work has been taken from the
% following paper:

% 'Residual change detection using low-complexity sequential quantile
% estimation', D. Jung, E. Frisk and M. Krysander

clc
clear all

%% Training reference distributions on training residual data 

% To obtain the set of reference probabilty distributions, the residual
% between the datamodel and the training Stage 2 data (Stage 2 below 25 K)
% has been chosen. The prior assumption is that this residual data set 
% contains the residual readings when the system is stable. 
% Therefore, there shall be changes in the probability distributions within
% this data set as well but all these changes correspond to pdfs that would
% exist in the stable system regime.

load('resStage2Model.mat');
resSignal = resStage2Eqn(:,1);
plot(resSignal);
title('Residual between Stage 2(below 25 K) and trained data model');
ylabel('Residual');
xlabel('Sample No.');

%% Sequential Quantile Estimation of Residual Signal

% The method works by computing the pdf of the residual at every time instant. 
% The algorithm used is called Sequential Quantile Estimation and is
% explained in the paper refered in the description.


% In the following loop, the distribution of residual is estimated for each
% time step and is stored in seqQuantEstDC.

tau = [0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99];
noOfQuants = size(tau,2);

seqQuantEstDC = zeros(size(resSignal,1),size(tau,2));
N = 100;

for k = 1:noOfQuants
    
    seqQuantEst = 0;
    n = 0;
    m = 0;

    for i = 1:size(resSignal,1)
        [seqQuantEst,n,m] = seqQuantile(resSignal(i,1),seqQuantEst,n,m,tau(1,k),N); 

        seqQuantEstDC(i,k) = seqQuantEst;
    end
end

estQuant = zeros(1, noOfQuants);
for i = 1:noOfQuants
    estQuant(1,i) = mean(seqQuantEstDC(:,i));
end

plot(resSignal); hold on; plot(seqQuantEstDC);
title('Estimated quantiles superimposed on residual plot');
ylabel('Residual');
xlabel('Sample No.');
%% K Means Classification of Sequential Quantile Estimates

% After obtaining the distribution for each time step, K-means
% classification is applied to obtain 5 clusters and hence 5 reference
% distributions.

% The final outcome of the training is the vector of centroids C which
% contains 5 reference distributions. A distribution is defined by a set of
% quantiles for the data set. In this case, 11 quantiles have been chosen
% to define the quantile set.

trainDataSeqQuantEst = seqQuantEstDC(1:size(seqQuantEstDC,1),:);

K = 5;
[id, C] = kmeans(trainDataSeqQuantEst, K);

for i = 1:size(id,1)
   distFromCentroid(i,1) = norm(seqQuantEstDC(i,:) - C(id(i,1),:)); 
end

% Plot trained Cumulative Distribution Functions

for i = 1:K
    plot(fliplr(C(i,2:11)), tau(2:11));
    hold on
end

title('Trained reference Cumulative Distribution Functions');
xlabel('Residual');
ylabel('Probability');
legend('1','2','3','4','5');

%% Plot different data samples of distributions
for i = 1:K
    ind = find(id == i);
    plot(resSignal(ind));
    hold on
end
legend('1','2','3','4','5');
title('Residual data as per reference CDF');
xlabel('Sample No.');
ylabel('Residual');
% temps_read.m subroutine
% Brings in string of monthly temperatures to calculate distribution,
% growth equation parameters, etc to feed into COMBO.

% UPDATED 3/28/11 to utilize temperature .txt files provided by Bob B. and
% stored in C:\projects\COMBO\Matlab\2011.03.28.subroutines\Temps_hist

% NOTE: THIS VERSION OF CODE WAS USED FOR DEBUGGING BLEACHING ISSUES
% ENCOUNTERED IN RUNNING FULLY MODULARIZED VERSION OF CODE WITH VALUATION.
% TAKEN DIRECTLY FROM 3/28 FOLDER, BUT FUNCTION CALL WAS MODIFIED TO
% REFLECT NEW INPUTS IN MOST RECENT VERSION OF COMBOV3.
function [baselineT,growthparams,monthlytemps,base_3month] = temps_readv5(temps,startDate,endDate,plotsave)

%monthlyT = xlsread('monthlyT.xls');

% ************************************************************************
% ******************** BASELINE TEMPERATURE READ-IN **********************
% ************************************************************************
% Format requirements: Data should be complete years, monthly data. Column
% 1 should be the monthly temperature data from the time period of
% interest. Column 2 should be the start and end year for the data 

% USER INPUT: select number of years to use from record to construct
% statistics on historical data and monthly averages.

%monthlength = input('Enter number of years to use for monthly averages (default = 10) ');
monthlength = 20;

% CHANGE 4/5/11: numyears becomes avgyears for calculating monthly
% averages, and distyears for calculating distribution
recordyears = [startDate:endDate];
numyears = length(recordyears);
% % Distyears = 30: Distribution is calculated for 1950-1979.
%distyears = 30;          

% Calculate monthly averages from subset of all data. Reshape first.
tempyears1 = reshape(temps,12,[]); % reshape so rows = months; columns = yrs
monthlytemps = mean(tempyears1(:,numyears-monthlength:numyears),2);
monthlytemps_max = max(tempyears1(:,numyears-monthlength:numyears),[],2);
monthlytemps_min = min(tempyears1(:,numyears-monthlength:numyears),[],2);
% monthlytemps = mean(tempyears1(:,numyears-monthlength:numyears),2);
% monthlytemps_max = max(tempyears1(:,numyears-monthlength:numyears),[],2);
% monthlytemps_min = min(tempyears1(:,numyears-monthlength:numyears),[],2);

% ************************************************************************
% ********************* BASELINE DISTRIBUTION CALCULATOR *****************
% ************************************************************************
% Calculate 3 month moving average T, and baselineT distribution for use in
% calculating threshold bleaching probabilities
% Note: "filter" algorithm does not treat endpoints properly. Need to shift
% vector up by one, then replace first and last elements with raw temps for
% endpoints (see calculation of thrmo_avgT below)
thrmo_avgT1 = filter(ones(1,3)/3,1,temps);
thrmo_avgT = [temps(1);thrmo_avgT1(3:length(temps));temps(length(temps))];
thrmo_tempyears1 = reshape(thrmo_avgT,12,[]); % Reshape so rows = months
% CHANGE from earlier versions: take only most recent data from 3 month avg.
thrmo_tempyears = thrmo_tempyears1(:,1:numyears);
%thrmo_tempyears = thrmo_tempyears1(:,1:distyears);
thrmo_means = mean(thrmo_tempyears);
thrmo_maxs = max(thrmo_tempyears);
% spit out the mean three month maximum temperature
base_3month = mean(thrmo_maxs);
baselineT = thrmo_maxs-mean(thrmo_maxs);

% ************************************************************************
% *********************** GROWTH COEFFICIENT CALCULATOR ******************
% ************************************************************************
% Note: calculations are on raw monthly temperatures (not 3 month averages) 
% first, clip reshaped raw data to period of interest
tempyears = tempyears1(:,1:numyears); 

% Column vectors of annual mean, min, and max
annmeans = mean(tempyears);
annmins = min(tempyears);
annmaxs = max(tempyears);
annstdev = std(annmaxs);    % standard deviation of annual max T

% Calculate parameters for growth equations
x0 = mean(annmins)-5;
x1 = mean(annmaxs)-2*annstdev;
x2 = mean(annmaxs)+5;

% Populate coefficient matrix, A
A = ones(4);
A(:,1) = [x0^3 x1^3 x2^3 3*x1^2]';
A(:,2) = [x0^2 x1^2 x2^2 2*x1]';
A(:,3) = [x0 x1 x2 1]';
A(:,4) = [1 1 1 0]';

% Populate right hand side (solution vector)
B = [0 1 0 0]';

% Calculate coefficients
growthparams = A\B;

 % visual display of temperatures: DISABLED
 figure(5)
 subplot(2,1,1)
 clf
 plot(monthlytemps,'r')
 hold on
 plot(monthlytemps_max,'r-.')
 plot(monthlytemps_min,'r-.')
 xlabel('Month')
 ylabel('Monthly average T')
 title('Monthly mean T for recent historical record')
 axis([1 12 min(monthlytemps_min) max(monthlytemps_max)])
 grid on
 
 subplot(2,1,2)
 clf
 plot(recordyears,annmeans,'r.')
 xlabel('Year')
 ylabel('Annual Average T')
 title('Annual average T for complete historical record')
 axis([startDate endDate min(annmeans) max(annmeans)])
 grid on



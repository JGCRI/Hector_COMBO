% COMBO model implementation in Matlab
% Designed to mimic original .xls version of COMBO but to allow more flexibility
%  to ingest alternative location information, climate scenarios, etc.
% First coded by C. Wobus, Stratus Consulting, November 2010.

% Version 2 of COMBO code: integrating all subroutines into a master code
% that allows for maximum flexibility and user input. 

% Major changes to Version 3 of COMBO code (completed 2/2/11): 
% 1. Includes subroutine temps_readv3.m to calculate growth coefficients,
%    calculate stats on historical data, and construct monthly patterns
% 2. Modified User_inputs to remove growth equation constants and other
%   unneeded variables that are now calculated instead of input.

% Additional changes to Version 4 of COMBO code (completed 2/2012):
% 1. Add three additional mortality events with TTd, TTe, TTf to User
% inputs_v4.xls file

% COMBO code renamed to go with combo_loopv5.m and valuation2.m. Minimal
% changes to COMBO code itself.

function [DateVec, Cov, Cov1, Cov2, PbleachA, PbleachB, PbleachC, ...
    PbleachD, PbleachE, PbleachF, bleach_index] = ...
    combov5(startDate,endDate,latlong,temps,bigloop,...
    path1,path2,path4,path4txt,deltaTTa,plotsave,flexTTa,placename);

% User input for type of growth (scaled or flat):
%growthtype = input('Enter type of growth to use (f = flat or s = scaled): ','s');
growthtype = 's';

% *********************************************************************
% Load Excel file with user input parameters. Query this file later to get
% parameters of interest in locations of interest.
% *********************************************************************
if placename == 'HI'
    thresh_inputs = xlsread('User_inputsv5.xls','HI_thresh');
elseif placename == 'PR'
    thresh_inputs = xlsread('User_inputsv5.xls','PR_thresh');
elseif placename == 'FL'
    thresh_inputs = xlsread('User_inputsv5.xls','FL_thresh');
end

% Get latitude, longitude, salinity and alkalinity from variable "latlong"
lat = latlong(bigloop,1);
long = latlong(bigloop,2);
sal = latlong(bigloop,3); % Salinity of the sample (Default = 35 ppt)
par1 = latlong(bigloop,4);

% *********************************************************************
% Call function temps_readv5 to process monthly temperature data for
% calculating growth and mortality parameters, generating monthly
% temperatures, and creating statistics on past temperature variability.
% NOTE: monthly temperatures now come from raw R1 and R2 datasets,
% available by cell, and called in shell code combo_loop.m at lines 24, 74 
% *********************************************************************
[baselineT,growthparams,monthlytemps,base_3month] = temps_readv5(temps,startDate,endDate,plotsave);
baseT = mean(monthlytemps);


% *********************************************************************
% Call function scengen_extractv3 to cull relevant data from Scengen *.OUT
% files and corresponding MAGICC CO2 data. These data feed into the
% detailed_tempsv5 subroutine below.
% *********************************************************************
[temp_CO2,scenarios] = scengen_extractv5(baseT,path1,path2,path4,...
    path4txt,lat,long);


% *********************************************************************
% Call function detailed_temps.m to calculate full temperature scenario
% *********************************************************************
[DateVec,detailedtemps] = detailed_tempsv5(monthlytemps,scenarios,baseT,plotsave);
FutureT = detailedtemps;

% Calculate statistics based on historical data
% NOTE: These stats are based on the deviations from mean.
meanT = mean(baselineT);
stdevT = std(baselineT);
minT = min(baselineT);
maxT = max(baselineT);

% ************************************************************
% **Mortality and Growth parameters (user interface Line 21)**
% ************************************************************
Gmax = thresh_inputs(1);
Gb = Gmax/12;
Mmax = thresh_inputs(2);
Mb = Mmax/12;

% ************************************************************
% ***Bleaching Mortality Module (User Interface Line 26-35)***
% ************************************************************
% Threshold temperatures (User Interface A29, A31, A33)
if flexTTa == 'y'
    % These lines are for user-input threshold temperatures
    TTa = base_3month + deltaTTa;
    TTb = TTa + 0.2;
    TTc = TTb + 0.2;
    % NEW 2/2012: Add three additional mortality events
    TTd = TTc + 0.2;
    TTe = TTd + 0.2;
    TTf = TTe + 0.2;
else
    % These lines are for default threshold temperatures
    TTa = thresh_inputs(3);
    TTb = thresh_inputs(4);
    TTc = thresh_inputs(5);
    % NEW 2/2012: Add three additional mortality events
    TTd = thresh_inputs(18);
    TTe = thresh_inputs(19);
    TTf = thresh_inputs(20);
end
% % DeltaT adaptation options (converted to degrees/month)
% NEW 2/2012: Add three additional mortality events dTTd, dTTe, dTTf
dTTa = thresh_inputs(6)/12;
dTTb = thresh_inputs(7)/12;
dTTc = thresh_inputs(8)/12;
dTTd = thresh_inputs(21)/12;
dTTe = thresh_inputs(22)/12;
dTTf = thresh_inputs(23)/12;
% Bleaching Factors (User Interface C29 - E29)
% NEW 2/2012: Add three additional mortality events
BFa = thresh_inputs(9);
BFb = thresh_inputs(10);
BFc = thresh_inputs(11);
BFd = thresh_inputs(24);
BFe = thresh_inputs(25);
BFf = thresh_inputs(26);
% Trigger Probabilities
TPa = thresh_inputs(12);
TPb = thresh_inputs(13);
TPc = thresh_inputs(14);
TPd = thresh_inputs(27);
TPe = thresh_inputs(28);
TPf = thresh_inputs(29);
% Mortality Factors
MFa = thresh_inputs(15);
MFb = thresh_inputs(16);
MFc = thresh_inputs(17);
MFd = thresh_inputs(30);
MFe = thresh_inputs(31);
MFf = thresh_inputs(32);

% Threshold Exceedances (Columns E-G on "3 Threshold Calculator")
% NOTE: DOES NOT INCORPORATE DELTATTa, Tb and Tc here to include adaptation
% option. See Combo model 3 Threshold Calculator, Columns B-D
Tdiffa = TTa - FutureT;
Tdiffb = TTb - FutureT;
Tdiffc = TTc - FutureT;
Tdiffd = TTd - FutureT;
Tdiffe = TTe - FutureT;
Tdifff = TTf - FutureT;


% ************************************************************
% ****************** CO2 concentration inputs ****************
% ************************************************************
% Omega sensitivity values (from A16:B16 of sheet User Interface)
% also calculates sensitivity values tabulated in Combo Calculator X4:Y12.
% This may need to be modified...
omega_sens1 = 0.2;
omega_sens2 = 0.4;
sens1 = omega_sens1*(4.6/3.73);
sens2 = omega_sens2*(4.6/3.73);

% ************************************************************
% ************ Execute "3 Threshold Calculator" **************
% ************************************************************
% Create an empirical cdf based on historical temperature data
% f,x are the percentiles and values associated with baselineT data, where
% baselineT is populated with temperature deviations from average.
% NOTE: f,x will not necessarily have the same dimensions as baselineT
% since there may be repeat values. Currently use a loop to populate the
% probability vectors, but this may be worth improving for speed.
[f,x]=ecdf(baselineT);
% debugging: look at ecdf to see how it jives with previous versions.
if plotsave=='y'
    figure(6)
    hold on
    plot(x,f,'g','LineWidth',2);
    grid on
    title('Empirical CDF for baselineT data')
end

threshAprob = zeros(1,length(Tdiffa));
threshBprob = zeros(1,length(Tdiffb));
threshCprob = zeros(1,length(Tdiffc));
threshDprob = zeros(1,length(Tdiffd));
threshEprob = zeros(1,length(Tdiffe));
threshFprob = zeros(1,length(Tdifff));
for i=1:length(Tdiffa);
    %loop to fill Threshold A vector (col H); leave at 0 if Tdiffa<minT
    if Tdiffa(i)>maxT
        threshAprob(i) = 1;
    elseif minT<Tdiffa(i)<maxT
        % Find temperature in Tdiffa closest to baselineT cdf
        [~,prob_index1]=min(abs(Tdiffa(i)-x));
        threshAprob(i) = f(prob_index1);
    end
    
    %loop to fill Threshold B vector (col K); leave at 0 if Tdiffb<minT
    if Tdiffb(i)>maxT
        threshBprob(i) = 1;
    elseif minT<Tdiffb(i)<maxT
        % debugging: threshAprob(i)=99;
        [~,prob_index2]=min(abs(Tdiffb(i)-x));
        threshBprob(i) = f(prob_index2);
    end
    
    %loop to fill Threshold C vector (col N); leave at 0 if Tdiffc<minT
    if Tdiffc(i)>maxT
        threshCprob(i) = 1;
    elseif minT<Tdiffc(i)<maxT
        % debugging: threshAprob(i)=99;
        [~,prob_index3]=min(abs(Tdiffc(i)-x));
        threshCprob(i) = f(prob_index3);
    end
    
    %loop to fill Threshold D vector (col N); leave at 0 if Tdiffd<minT
    if Tdiffd(i)>maxT
        threshDprob(i) = 1;
    elseif minT<Tdiffd(i)<maxT
        % debugging: threshAprob(i)=99;
        [~,prob_index4]=min(abs(Tdiffd(i)-x));
        threshDprob(i) = f(prob_index4);
    end
    
    %loop to fill Threshold E vector (col N); leave at 0 if Tdiffe<minT
    if Tdiffe(i)>maxT
        threshEprob(i) = 1;
    elseif minT<Tdiffe(i)<maxT
        % debugging: threshAprob(i)=99;
        [~,prob_index5]=min(abs(Tdiffe(i)-x));
        threshEprob(i) = f(prob_index5);
    end
    
    %loop to fill Threshold F vector (col N); leave at 0 if Tdifff<minT
    if Tdifff(i)>maxT
        threshFprob(i) = 1;
    elseif minT<Tdifff(i)<maxT
        % debugging: threshAprob(i)=99;
        [~,prob_index6]=min(abs(Tdifff(i)-x));
        threshFprob(i) = f(prob_index6);
    end
    
end

% Populate vectors of non-occurrence prob. w/ bleaching factors
% (columns J, M, and P in sheet "3 threshold calculator")
threshA_bleach = 1-(1-threshAprob)*BFa;
threshB_bleach = 1-(1-threshBprob)*BFb;
threshC_bleach = 1-(1-threshCprob)*BFc;
threshD_bleach = 1-(1-threshDprob)*BFd;
threshE_bleach = 1-(1-threshEprob)*BFe;
threshF_bleach = 1-(1-threshFprob)*BFf;

% Calculate cumulative exceedance probability vectors (columns T-V in sheet
% "3 Threshold Calculator")
PbleachA = 1-cumprod(threshA_bleach);
PbleachB = 1-cumprod(threshB_bleach);
PbleachC = 1-cumprod(threshC_bleach);
PbleachD = 1-cumprod(threshD_bleach);
PbleachE = 1-cumprod(threshE_bleach);
PbleachF = 1-cumprod(threshF_bleach);
%****************** END "3 Threshold Calculator" *************

% ************************************************************
% ************ Execute "Trigger Calculator" **************
% ************************************************************
triggerPullA = PbleachA>=TPa;
triggerPullB = PbleachB>=TPb;
triggerPullC = PbleachC>=TPc;
triggerPullD = PbleachD>=TPd;
triggerPullE = PbleachE>=TPe;
triggerPullF = PbleachF>=TPf;

% find index of first nonzero element in each triggerPull vector (note that
% indA, indB and indC correspond to the indices in Trigger Calculator
% columns D, G, and J, respectively, where "Mortality Addition Comb Cover
% Surv" is populated
indA = find(triggerPullA,1);
indB = find(triggerPullB,1);
indC = find(triggerPullC,1);
indD = find(triggerPullD,1);
indE = find(triggerPullE,1);
indF = find(triggerPullF,1);

% Save a vector of all bleaching dates to pass to shell code
bleach_index = [indA indB indC indD indE indF];

% Create Survival Factor vector (Column F of Combo Calculator)
SurvFact = ones(size(FutureT));
SurvFact(indA) = 1 - MFa;
SurvFact(indB) = 1 - MFb;
SurvFact(indC) = 1 - MFc;
SurvFact(indD) = 1 - MFd;
SurvFact(indE) = 1 - MFe;
SurvFact(indF) = 1 - MFf;
%****************** END "Trigger Calculator" *************

% ************************************************************
% *************** Execute Combo Calculator *******************
% ************************************************************

% Create a column vector of mortality based on user-defined Mb value
Mtotal = Mb.*ones(size(DateVec));

%Growth equation (equivalent to column G in sheet Combo Calculator)
Geqn = growthparams(1).*FutureT.^3 + growthparams(2).*FutureT.^2 + ...
     growthparams(3).*FutureT + growthparams(4);

%SST multiplier: populate column H of Combo Calculator with seasonally
%variable SST-dependent growth.
% Set initial vectors of SST, Gnet and Cov to fill in loop below.
SST_Gmonth = zeros(size(FutureT));
Gnet = zeros(size(FutureT));
Cov = ones(size(FutureT));
if growthtype == 'f'
    SST_Gmonth(1:12) = Gb;
elseif growthtype == 's'
    SST_Gmonth(1:12) = Geqn(1:12).*Gmax/sum(Geqn(1:12));
end

% NOTE: Try to eliminate this loop using cumsum...
for i = 13:length(SST_Gmonth)
    SST_Gmonth(i) = SST_Gmonth(i-12)*(Geqn(i)/Geqn(i-12));
    Gnet(i) = SST_Gmonth(i) - Mtotal(i);
    Cov(i) = SurvFact(i).*Cov(i-1)*(1+Gnet(i));
end

% ********************** CO2-1 mortality ************************
% Calculate columns L-P in Combo Calculator. NOTE: this assumes for the
% moment that everything is input as 5-year inputs. This should still work
% if inputs are in 1-year increments...interpolation will then be by month.
% Interpolate between 5-year increments to fill in omega decline 1 (Column
% L of Combo Calculator). NOTE: original spreadsheet does not interpolate
% but rather copies and pastes values. Output will be slightly different.
% Call function omega_slope to fill the omegaslope matrix:
[omegaslope] = omega_slopev5(temp_CO2,par1,sal);
omega_dec1 = interp1(omegaslope(:,1),omegaslope(:,5),DateVec,'linear','extrap');
monthly_dec1 = sens1*(omega_dec1/12);
% zero out first 12 entries of omega_dec1 and monthly_dec1
omega_dec1(1:12,1) = 0;
monthly_dec1(1:12,1) = 0;
CO21_Gnet = SST_Gmonth.*(1+cumsum(monthly_dec1)) - Mtotal;
Cov1 = ones(size(FutureT));

% ********************** CO2-2 mortality ************************
% Calculate columns Q-U in Combo Calculator. Note that omega decline does
% not change between the two CO2 mortality calculations; only the
% sensitivity parameter (cells A16 and B16 of User Interface, which
% translate to cells X13 and X14 in Combo Calculator)
monthly_dec2 = sens2*(omega_dec1/12); 
% zero out first 12 entries of monthly_dec2
monthly_dec2(1:12,1) = 0;
CO22_Gnet = SST_Gmonth.*(1+cumsum(monthly_dec2)) - Mtotal;
Cov2 = ones(size(FutureT));

% Loop to calculate CO2-dependent growth
% NOTE: Try to eliminate this loop using cumsum...or embed both Cov1 and
% Cov into the same loop for speed.
for i = 13:length(SST_Gmonth)
    Cov1(i) = SurvFact(i).*Cov1(i-1)*(1+CO21_Gnet(i));
    Cov2(i) = SurvFact(i).*Cov2(i-1)*(1+CO22_Gnet(i));
end


% End function
end







%function [detailedtemps] = detailed_temps(monthlytemps,scenarios)

% Creates a vector of future temperatures (detailedtemps) based on
% user-input for monthly temperature patterns and model output from future
% temperature scenarios. 

function [DateVec,detailedtemps] = detailed_tempsv5(monthlytemps, ...
    scenarios, baseT, plotsave)

% Create vector of future dates in decimal years 
DateVec = [2000:(1/12):(2100+(11/12))]';

% Create vector of monthly deviations from average temperature:
diffs = monthlytemps-baseT;

% Calculate the annual increment of temperature change. NOTE: this requires
% that all dates are expressed in decimal years.
% 11/30/10: this is bunk. No need to calculate increments as in .xls
% version. Can simply interpolate 5-yr changes into monthly values.
temp_incr1 = diff(scenarios(:,2))./diff(scenarios(:,1));

% interpolate so that temperature increases are applied across all months. 
% NOTE: this is slightly different than the original COMBO model in that it
% spreads projected temperature increases out by month, rather than using a
% "stepped" increase in temp every x years. But the method used here is
% also likely to be more realistic.
scenarios_ann = interp1(scenarios(:,1),scenarios(:,2),DateVec,'linear','extrap');

% Add monthly variability to incremental temperature increases 
% First need to make a vector with repeating monthly temperature variations
% that's the same length as DateVec...
monthlypattern = repmat(diffs(:,1),length(DateVec)/12,1);
temp_incr3 = scenarios_ann + monthlypattern;

% need to add increments to baseline temperature
detailedtemps = baseT + temp_incr3;

if plotsave=='y'
figure(4)
%clf
hold on
plot(DateVec,detailedtemps,'g')
xlabel('Year')
ylabel('Projected Temperature (C)')
grid on
end

end




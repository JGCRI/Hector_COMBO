% Function omegaA_lookup
% To be used in concert with code "combo.m" for determining omegaslope
% values given projected timeseries of temperature and CO2 concentrations.

% NEW VERSION 3/1/11: Calculates Omega values using J. Kleypas CO2SYS code,
% rather than using lookup table of Omega values. Requires as input values
% for alkalinity and salinity, in addition to the CO2 and temperature data
% provided by MAGICC/SCENGEN output

function [omegaslope] = omega_slopev5(temp_CO2,par1,sal)

omegaA = zeros(length(temp_CO2),1);

% Fill in values for salinity, alkalinity and parameter types for CO2SYS
par1type=    1;    % The first parameter supplied is of type "1", which means "Total Alkalinity"
par2type=    4;    % The second parameter supplied is of type "4", which means "pCO2"
% par2 is CO2 concentration, which is an output from MAGICC
% tempin  is temperature, which is an output from SCENGEN
presin  =    1;    % Pressure at input conditions (presin=1; upper meter of water column)
presout =    1;    % Pressure at output conditions (presout=presin)
sil     =    0;    % Concentration of silicate  in the sample (in umol/kg)
po4     =    0;    % Concentration of phosphate in the sample (in umol/kg)
pHscale =    1;    % pH scale at which the input pH is reported ("1" means "Total Scale")
k1k2c   =    4;    % Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
kso4c   =    1;    % Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")


% Loop to fill vector omegaA using lookup table. Use temperature and CO2
% values from lookup table that are closest to those from output.
for k = 1:length(temp_CO2);
    %temp_co2 = three columns: Year CO2 Temperature 
    par2 = temp_CO2(k,2);
    tempin = temp_CO2(k,3);
    tempout = tempin;
    % This is the line that calls CO2SYS. 
    % CO2SYS() calculates an array of 52 values
    % COMBO needs OmegaArout to determine calcification rates
    % which is the 31st value of the array
    CO2SYSOutput=CO2SYS(par1,par2,par1type,par2type,sal,tempin,tempout,...
        presin,presout,sil,po4,pHscale,k1k2c,kso4c);
    omegaA(k) = CO2SYSOutput(31);
end

temp_CO21 = [temp_CO2 omegaA];

% Calculate omegaslope from omega values and dates. NOTE: this algorithm
% uses a forward difference, whereas the original COMBO .xls sheet uses a
% centered difference. May slightly alter output but the FD is
% substantially less expensive computationally. Note that length of
% omegaslope vector is one shorter than length of temp_CO2 vector.
omegaslope1 = diff(omegaA)./diff(temp_CO2(:,1));

% fill omegaslope matrix for export back to combo.m
omegaslope = [temp_CO21(1:length(omegaslope1),:) omegaslope1];
end







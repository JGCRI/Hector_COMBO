% combo_loopv5.m

% Shell code to loop through individual columns of user-input baseline
% temperature data, and to run the COMBO code on each individual cell.
% Outputs from this code are strings of projected decreases in coral cover,
% for each cell from which baseline temperature data were extracted.

% Version 4/18/11: Add in functionality to append file "latlong" with %
% coral cover under different scenarios in different years (choose 2050 and
% 2100) for different scenarios

% Modifications 5/5/11: Automate as many steps as possible to minimize
% confusion with other users.

% Modifications 5/26/11: Add valuation subroutine to calculate value of
% coral reef under different scenarios. Update calculation of TTa, TTb,
% and TTc so they're 0.9, 1.1, and 1.3 degrees above baseline (needs to be
% hard-wired in code combov4.m: comment out lines 87-89, uncomment 83-85)
% Variable with values is called "discvalues_all" and is organized by cell
% (columns) and by year (rows).

% Modifications 7/12/11: Add in plots to compare bleaching scenarios with
% and without policy options, and total value of difference.

% Modifications 2/2012: Add additional three mortality events

% Clear all variables 
clear

% initialize variables to save value and cover data at end.
discvalues_all = [];
covvalues_all = [];
cover_all = [];
cover_Tonly_all = [];
cover4val = [];

% tempseries is a variable that points to the correct folder for extracting
% temperature data, depending on which dataset is being used.
%temptype = input('Enter temperature dataset to use (R2, R1R2): ','s');
temptype = 'R1R2';

% User input for location:
placename = input('Enter location (HI, PR, or FL; ALL CAPS): ','s');
placefile = strcat(placename,'_locs.xls');
tempseries = strcat(placename,'_',temptype);
latlong = xlsread(placefile);

% Thresholds for bleaching
disp('Model is sensitive to default deltaT for first threshold')
disp('We recommend first running with default thresholds in user input file')
disp('Alternatively, thresholds can be entered manually here for sensitivity analysis')
flexTTa = input('Enter thresholds manually? (y/n)','s');
if flexTTa == 'y'
    deltaTTa = input('Enter deltaT for first threshold [Start with 0.75(FL) or 0.9(HI, PR) ]');
else
    deltaTTa = -9999;
end

% Set file paths for codes, SCENGEN output, and temperatures
path1 = pwd;
% DEBUGGING: use path to old SCENGEN output
path2 = strcat(path1,'\Scengen_EPA');
path3 = strcat(path1,'\Temps_hist\',tempseries);

% Fix path to relevant SCENGEN files to extract model output (BAU vs POL)
disp('')
disp('NOTE: for valuation, you must run policy scenario first as baseline')
disp('After policy scenario is run, re-run with no policy for valuation')
disp('')
blyesno = input('Save policy output as baseline for value comparison? (y/n) ','s');
path4txt = input('Enter Scenario to run (NoParis, Paris, ParisPlus, ParisRef): ','s');

% User inputs for recreational use and nonuse values:
% PUT THESE INTO USER INPUTS SPREADSHEETS
disp('')
[valueinputs,types,~] = xlsread('User_inputsv5','values_all');
useRUdef = input('Use default recreational use values [HI - $1741M; FL - $2039M; PR - $0.24M]? (y/n) ','s');
if useRUdef == 'y'
    if placename == 'FL'
        recval_min = valueinputs(1,1);
        recval = valueinputs(1,2);
        recval_max = valueinputs(1,3);
    elseif placename == 'HI'
        recval_min = valueinputs(3,1);
        recval = valueinputs(3,2);
        recval_max = valueinputs(3,3);
    elseif placename == 'PR'
        recval_min = valueinputs(5,1);
        recval = valueinputs(5,2);
        recval_max = valueinputs(5,3);
    end
elseif useRUdef == 'n'
    recval = input('Enter Recreational Use value for reef site ($M): ');
end

useNUVdef = input('Use default Non-Use Values [HI: $32,444M; FL: $458M; PR: $305M] (y/n) ','s');
if useNUVdef == 'y'
    if placename == 'FL'
        nuval_min = valueinputs(2,1);
        nuval = valueinputs(2,2);
        nuval_max = valueinputs(2,3);
    elseif placename == 'HI'
        nuval_min = valueinputs(4,1);
        nuval = valueinputs(4,2);
        nuval_max = valueinputs(4,3);
    elseif placename == 'PR'
        nuval_min = valueinputs(6,1);
        nuval = valueinputs(6,2);
        nuval_min = valueinputs(6,3);
    end
elseif useNUVdef == 'n'
    nuval = input('Enter Non-Use Value for reef site ($M): ');
end

valtype = input('What values should be used for valuation routine (Rec, Nonuse, or Total)? ','s');

% Discounting parameters go here:
base2015 = input('Use 2015 as base year (default)? ','s');
if base2015 == 'y'
    baseyr = 2015;                          % Base year for calcs
elseif base2015 == 'n'
    baseyr = input('Enter base year for discounting: ');
end
disc05 = input('Use discount rate of 5% (default)? ','s');
if disc05 == 'y'
    disc = 0.05;                            % discount rate (%)
elseif disc05 == 'n'
    disc = input('Enter discount rate as decimal (e.g., 0.05): ');
end

% Allow suppression of plot outputs to speed analysis
disp('Show all plots (e.g., cover through time, bleaching probabilities, maps?')
plotsave = input('(y/n): ','s');

path4 = strcat(path2,'\',path4txt);

% CHANGE from previous versions: Start and end dates are all the same at
% 1900 to 2010. Longer timeseries are used for calculating all relative
% values, and short (10-20 year) timeseries are used for calculating
% monthly patterns.
% Currently Hawaii has only 1950-2010. Others have 1900-2010
startDate = 1950;
endDate = 2010;
% hitest = strcmp(placename,'HI');
% if hitest
%     startDate = 1950;
%     endDate = 2010;
% else
%     startDate = 1900;
%     endDate = 2010;
% end

% Look through directory to find all temperature output files
tempfiles = strcat(path3,'\*.txt');
tfiles = dir(tempfiles);
numtempfiles = size(tfiles,1);
cellnums = [];

% % Initialize output file
% HI_output = zeros(numtempfiles,11);

% ************************************************************************
% ************* BEGIN THE BIG LOOP THRU ALL LOCATIONS ********************
% ************************************************************************
for bigloop=1:numtempfiles
    %bigloop = 1;
    
    % change to file path with relevant temperature data and load it.
    cd(path3)
    
    temps = [];
    tempfname = tfiles(bigloop).name
    % create strings for data and plot filenames
    namespace = find(isspace(tempfname));
    tlocn = tempfname(namespace+1:length(tempfname)-4);    % extract cell name
    % Keep track of cell numbers in order of analysis
    cellnums = [cellnums str2num(tlocn(3:length(tlocn)))];
    fig2 = [tempseries,'_',path4txt,'_prob_',tlocn,'.pdf'];
    fig3 = [tempseries,'_',path4txt,'_cov_',tlocn,'.pdf'];
    % Open relevant temperature file and save data to variable "temps"
    fid = fopen(tfiles(bigloop).name);
    temps = [temps; fscanf(fid, '%*s %f')];
    fclose(fid);
    
    % switch back to root file path
    cd(path1)
    
    % NEED TO EXPORT STARTDATE, ENDDATE, AND TEMPS TO REST OF CODE...
    [DateVec, Cov, Cov1, Cov2, PbleachA, PbleachB, PbleachC, ...
        PbleachD, PbleachE, PbleachF, bleach_index] = ...
        combov5(startDate,endDate,latlong,temps,bigloop,...
        path1,path2,path4,path4txt,deltaTTa,plotsave,flexTTa,placename);
    
    % ************************************************************
    % ******************  Plotting Subroutines *******************
    % ************************************************************
    % NOTE: these were previously contained within combov4 code. Moved here to
    % gather all relevant plots and variables in the shell code
    % Plot cumulative exceedance probabilities for three thresholds
    if plotsave=='y'
             figure(2)
             clf
             plot(DateVec,PbleachA,'LineWidth',2)
             hold on
             grid on
             plot(DateVec,PbleachB,'r','LineWidth',2);
             plot(DateVec,PbleachC,'g','LineWidth',2);
             axis([2000 2100 0 1])
             xlabel('Year')
             ylabel('Probability')
             title({['Cell: ',tlocn,'; Temp Series: ',tempseries,'; Scenario: ',path4txt,]; ...
                 'Model Calculated Cumulative Probability of Bleaching'});
             legend('Threshold A','Threshold B','Threshold C')
             eval (['print -dpdfwrite ',fig2])
        
        % Plot reef declines for all three calculation methods.
        figure(3)
        clf
        plot(DateVec,Cov,'b','LineWidth',2)
        hold on
        plot(DateVec,Cov1,'k','LineWidth',2)
        plot(DateVec,Cov2,'r','LineWidth',2)
        grid on
        axis([2000 2100 0 1.2])
        xlabel('Year')
        ylabel('Percent Initial Cover')
        title({['Cell: ',tlocn,'; Scenario: ',path4txt,]; ...
            'Percent Initial Cover Over Time'});
        legend('SST Only','CO2 Sensitivity 1','CO2 Sensitivity 2')
        eval (['print -dpdfwrite ',fig3])
    end
    
    % compile and save cover output for all three sensitivities
    cov_output = [DateVec Cov Cov1 Cov2];
    out_fname = strcat('PctInitCov_',tlocn,'_',path4txt,'.txt');
    save(out_fname,'cov_output','-ascii');
    
    % ******************* INITIAL COVER DATA GO HERE **********************
    % Source for initial cover data: see sheet "Cover" in each of the three
    % workbooks used as inputs: HI_locs.xls, PR_locs.xls, FL_locs.xls
    % ************************ Data for Hawaii ****************************
    % These data are no longer needed for valuation; however, they are needed
    % for making composite plots of percent cover by cell

    if placename == 'HI'
        [hinitcovers,hcellnames,~] = xlsread('User_inputsv5','HI_initcov');
        % Grab the initial cover data corresponding to cell of interest
        [rownum,~] = find(strcmp(tlocn,hcellnames));
        cov_init = hinitcovers(rownum);
        % NEW 5/08: save column vector of cover through time for CO2 sens1
        if ~isempty(rownum)
            cover_CO21 = Cov1.*cov_init;
            cover_all = [cover_all cover_CO21];
            cover_Tonly = Cov.*cov_init;
            cover_Tonly_all = [cover_Tonly_all cover_Tonly];
        end
        % ********************** Data for Puerto Rico ************************
    elseif placename == 'PR'
        [pinitcovers,pcellnames,~] = xlsread('User_inputsv5','PR_initcov');
        % Grab the initial cover data corresponding to cell of interest
        [rownum,~] = find(strcmp(tlocn,pcellnames));
        cov_init = pinitcovers(rownum);
        % NEW 5/08: save column vector of cover through time for CO2 sens1
        if ~isempty(rownum)
            cover_CO21 = Cov1.*cov_init;
            cover_all = [cover_all cover_CO21];
            cover_Tonly = Cov.*cov_init;
            cover_Tonly_all = [cover_Tonly_all cover_Tonly];
        end
        % ********************** Data for Florida ******************************
    elseif placename == 'FL'
        [finitcovers,fcellnames,~] = xlsread('User_inputsv5','FL_initcov');
        % Grab the initial cover data corresponding to cell of interest
        [rownum,~] = find(strcmp(tlocn,fcellnames));
        cov_init = finitcovers(rownum);
        % NEW 5/08: save matrix of cover through time for CO2 sens1
        if ~isempty(rownum)
            cover_CO21 = Cov1.*cov_init;
            cover_all = [cover_all cover_CO21];
            cover_Tonly = Cov.*cov_init;
            cover_Tonly_all = [cover_Tonly_all cover_Tonly];
        end
    end
    
    % ADD HERE THE TEMPERATURE SENSITIVITY 1 AND 2 TO PASS TO VALUATION
    % ALSO NEED TO PASS ALONG COVER_AVG
    % Save matrix of all Cov1 values for valuation routine
    cover4val = [cover4val Cov1];
    
    % Find the years of all bleaching events. If statement below is to
    % ensure no logical indexing errors, since some cells may have NaN if
    % not all bleaching events are achieved.
    nobleach = find(isnan(bleach_index));
    if ~isempty(nobleach)
        nobleach_ind = min(find(isnan(bleach_index)));
        bleachdates = DateVec(bleach_index(1:nobleach_ind-1),1);
    else
        bleachdates = DateVec(bleach_index,1);
    end
    
    out_fname2 = strcat(tlocn,'_bleachdates_','_',path4txt,'.txt');
    save(out_fname2,'bleachdates','-ascii');
    

end

% Calculate average cover (linear average of all cells with initial cover
% data)
cover_avg = mean(cover_all,2);
cover_avg_Tonly = mean(cover_Tonly_all,2);
cover_output = [DateVec cover_avg cover_avg_Tonly];

% ************************************************************************
% **************** END BIG LOOP THROUGH TEMPERATURE FILES ****************
% ************************************************************************


% ************************************************************************
% ********************* CALL VALUATION SUBROUTINE ***********************
% ************************************************************************

% Input hard bottom initial cover values and calculate percentages
hardbottom = xlsread('User_inputsv5','val_cover');

% Figure out what value to export to valuation subroutine
if strcmp(valtype,'Total')==1
    tvalue_min = recval_min + nuval_min;
    tvalue = recval + nuval;
    tvalue_max = recval_max + nuval_max;
    tvalues = [tvalue_min tvalue tvalue_max];
elseif strcmp(valtype,'Rec')==1
    tvalue_min = recval_min;
    tvalue = recval;
    tvalue_max = recval_max;
    tvalues = [tvalue_min tvalue tvalue_max];
elseif strcmp(valtype,'Nonuse')==1
    tvalue_min = nuval_min;
    tvalue = nuval;
    tvalue_max = nuval_max;
    tvalues = [tvalue_min tvalue tvalue_max];
end

[disc_vals,ndisc_vals] = valuationv5(DateVec,cover4val,cellnums,hardbottom,placename,tvalues,baseyr,disc);

%
%     % Save policy data for valuation comparison to no policy options
%
% Save matrices discvalues_all
% Plot results for valuation comparisons
if strcmp(path4txt, 'ParisPlus')==1
    % filenames for output
    val_fname = strcat(path1, '\output\', placename,'_',path4txt,'_val.txt');
    cov_fname = strcat(path1, '\output\', placename,'_',path4txt,'_cov.txt');

    DateVec_yr = [2000:2100];
    
    figure(7)
    clf
    plot(DateVec,cover_avg,'b','LineWidth',2);
    title('Cover decline through time')
    xlabel('Year')
    ylabel('Average % Cover')
    title(['Reduced emissions scenario cover through time for ',placename])
    grid on
    axis([2000 2100 0 100])
    axis 'auto y'
    
    figure(9)
    clf
    plot(DateVec,cover_avg_Tonly,'r','LineWidth',2)
    hold on
    plot(DateVec,cover_avg,'b','LineWidth',2);
    legend('Temperature only','CO2 Sensitivity 1')
    xlabel('Year')
    ylabel('Average % Cover')
    title(['Reduced emissions scenario cover through time: CO2 vs Temp Only for ',placename])
    grid on
    axis([2000 2100 0 100])
    axis 'auto y'
    
    figure(10)
    clf
    plot(DateVec_yr,disc_vals(:,2),'b','LineWidth',2)
    hold on
    plot(DateVec_yr,disc_vals(:,1),'b-.')
    plot(DateVec_yr,disc_vals(:,3),'b-.')
    legend('Reduced emissions')
    grid on
    xlabel('Year')
    ylabel('Discounted Value (M 2007$)')
    title(['Reduced emissions scenario discounted value through time for ',placename])
    
    figure(11)
    clf
    plot(DateVec_yr,ndisc_vals(:,2),'b','LineWidth',2)
    hold on
    plot(DateVec_yr,ndisc_vals(:,3),'b-.')
    plot(DateVec_yr,ndisc_vals(:,1),'b-.')
    legend('Reduced emissions')
    grid on
    xlabel('Year')
    ylabel('Nominal Value (M 2007$)')
    title(['Reduced emissions scenario nondiscounted value through time for ',placename])

    % If policy data already saved, open saved data and save variables for
    % no policy scenarios.
else
    try
        % NEED TO RE-WRITE FOR NEW VARIABLE NAMES
        val_basename = strcat(path1, '\output\', placename,'_ParisPlus_val.txt');
        cov_basename = strcat(path1, '\output\', placename,'_ParisPlus_cov.txt');

        base_val = dlmread(val_basename);
        base_cov = dlmread(cov_basename);
        DateVec_yr = [2000:2100];
        
        % Grab discounted and nondiscounted values from policy run for
        % comparison with BAU run
        polval_disc_min = base_val(:,2);
        polval_disc = base_val(:,3);
        polval_disc_max = base_val(:,4);
        polval_ndisc_min = base_val(:,5);
        polval_ndisc = base_val(:,6);
        polval_ndisc_max = base_val(:,7); 
                
        BAUValMin = sum(disc_vals(:,1));
        BAUValMean = sum(disc_vals(:,2));
        BAUValMax = sum(disc_vals(:,3));
        
        POLValMin = sum(polval_disc_min);
        POLValMean = sum(polval_disc);
        POLValMax = sum(polval_disc_max);
        
        lostbenefit_min = POLValMin-BAUValMin;
        lostbenefit = POLValMean - BAUValMean;
        lostbenefit_max = POLValMax - BAUValMax;
        
        
        
        %         lostbenefit_min = sum(polval_disc_min-disc_vals(:,1));
        %         lostbenefit = sum(polval_disc-disc_vals(:,2));
        %         lostbenefit_max = sum(polval_disc_max-disc_vals(:,3));
        
        % filenames for output
        val_fname = strcat(path1, '\output\', placename,'_',path4txt,'_val.txt');
        cov_fname = strcat(path1, '\output\', placename,'_',path4txt,'_cov.txt');
        valdata_fname = strcat(path1, '\output\', placename,'_',valtype,'_','val_out.txt');
        
        % Output a file with data summarizing values for each scenario
        out_value_data = [BAUValMin BAUValMean BAUValMax; ...
            POLValMin POLValMean POLValMax; lostbenefit_min lostbenefit lostbenefit_max]
        save(valdata_fname,'out_value_data','-ascii');
        
        figure(12)
        plot(DateVec,base_cov(:,2),'b','LineWidth',2);
        hold on
        plot(DateVec,cover_avg,'r-.','LineWidth',2);
        legend('Reduced Emissions','BAU')
        title('Cover decline through time')
        xlabel('Year')
        ylabel('Average % Cover')
        title(['Reduced emissions vs BAU cover comparison for ',placename]);
        grid on
        axis([2000 2100 0 100])
        axis 'auto y'
        cover_plotname = strcat(path1, '\output\', placename,'_cover.pdf');
        eval (['print -dpdfwrite ',cover_plotname])
        
        figure(13)
        clf
        plot(DateVec,cover_avg_Tonly,'r','LineWidth',2)
        hold on
        plot(DateVec,cover_avg,'b','LineWidth',2);
        legend('Temperature Only','CO2 Sensitivity 1')
        xlabel('Year')
        ylabel('Average % Cover')
        title(['BAU scenario cover through time: CO2 vs Temp Only for ',placename])
        grid on
        axis([2000 2100 0 100])
        axis 'auto y'
        TvsCo2_plotname = strcat(path1, '\output\', placename,'_cover_TvsCO2.pdf');
        eval (['print -dpdfwrite ',TvsCo2_plotname])
        
        figure(14)
        clf
        plot(DateVec_yr,polval_disc,'b','LineWidth',2)
        hold on
        plot(DateVec_yr,disc_vals(:,2),'r-.','LineWidth',2)
        legend('Reduced emissions','BAU')
        % Plot POL error bounds
        plot(DateVec_yr,polval_disc_min,'b-.')
        plot(DateVec_yr,polval_disc_max,'b-.')
        % Plot BAU error bounds
        plot(DateVec_yr,disc_vals(:,1),'r-.')
        plot(DateVec_yr,disc_vals(:,3),'r-.')
        grid on
        xlabel('Year')
        ylabel('Discounted Value (M 2007$)')
        title({['Reduced emissions vs BAU valuation comparison for ',placename],...
            ['Sum of Lost Annual Benefits = $',num2str(round(lostbenefit)),'M']})
        discval_plotname = strcat(path1, '\output\', placename,'_discvaluation.pdf');
        eval (['print -dpdfwrite ',discval_plotname])
        
        figure(15)
        clf
        % Plot POL results with error bounds
        plot(DateVec_yr,polval_ndisc,'b','LineWidth',2)
        hold on
        plot(DateVec_yr,ndisc_vals(:,2),'r-.','LineWidth',2)
        legend('Reduced emissions','BAU')
        % Plot POL error bounds
        plot(DateVec_yr,polval_ndisc_min,'b-.')
        plot(DateVec_yr,polval_ndisc_max,'b-.')
        % Plot BAU error bounds
        plot(DateVec_yr,ndisc_vals(:,1),'r-.')
        plot(DateVec_yr,ndisc_vals(:,3),'r-.')
        grid on
        xlabel('Year')
        ylabel('Nominal Value (M 2007$)')
        title(['Reduced emissions vs BAU valuation comparison for ',placename]);
        ndiscval_plotname = strcat(path1, '\output\', placename,'_ndiscvaluation.pdf');
        eval (['print -dpdfwrite ',ndiscval_plotname])


    catch
        s1 = 'No policy data found for ';
        s2 = '. Please run policy scenario for ';
        s3 = ' before running valuation.';
        msg = [s1 placename s2 placename s3];
        disp(msg)
    end

end
% End "if" loop

% Save value and average cover files: Value file now has both discounted 
% and non-discounted value, by year; cover is by month.
% UPDATE 5/18/12: disc_vals and ndisc_vals are now 101x3. Middle column of 
% each is the old mean value.
values_time = [DateVec_yr' disc_vals ndisc_vals];
save(val_fname,'values_time','-ascii');
save(cov_fname,'cover_output','-ascii'); 
% script valuation2.m
% Revised 5/9/12 to calculate valuation as weighted average of all cells
% with initial hard bottom data
% Output is a single column vector of weighted average discounted values by
% year, from 2000-2100. variable name = 'disc_val'

function [disc_vals,ndisc_vals] = valuationv5(DateVec,cover4val,cellnums,hardbottom,...
    placename,tvalues,baseyr,disc);

% Find the correct hard bottom percentage array for valuation
if placename == 'HI'
    hard_pct = hardbottom(:,2)./sum(hardbottom((~isnan(hardbottom(:,2))),2));
elseif placename == 'FL'
    hard_pct = hardbottom(:,3)./sum(hardbottom((~isnan(hardbottom(:,3))),3));
elseif placename == 'PR'
    hard_pct = hardbottom(:,4)./sum(hardbottom((~isnan(hardbottom(:,4))),4));
end

% Parameters that don't change
%baseyr = 2007;                          % Base year for calcs
%disc = 0.03;                            % discount rate (%)
PVfact = 1./(1+disc).^(DateVec-baseyr); % Vector of PV factors


% Sum the discounted, weighted value for all cells in valuation routine.
% NOTE: The individual cell columns are not preserved. Just the weighted
% sum.
discvalue_min_allcells = 0;
discvalue_allcells = 0;
discvalue_max_allcells = 0;

ndiscvalue_min_allcells = 0;
ndiscvalue_allcells = 0;
ndiscvalue_max_allcells = 0;

for valloop = 1:length(cellnums)
    % Calculate cell-weighted discounted value through time
    discvalue_cell_min = cover4val(:,valloop).*PVfact.*hard_pct(cellnums(valloop)).*tvalues(1);
    discvalue_cell = cover4val(:,valloop).*PVfact.*hard_pct(cellnums(valloop)).*tvalues(2);
    discvalue_cell_max = cover4val(:,valloop).*PVfact.*hard_pct(cellnums(valloop)).*tvalues(3);
    
    discvalue_min_allcells = discvalue_min_allcells + discvalue_cell_min;
    discvalue_allcells = discvalue_allcells + discvalue_cell;
    discvalue_max_allcells = discvalue_max_allcells + discvalue_cell_max;
    
    % Calculate cell-weighted nondiscounted value through time
    ndiscvalue_cell_min = cover4val(:,valloop).*hard_pct(cellnums(valloop)).*tvalues(1);
    ndiscvalue_cell = cover4val(:,valloop).*hard_pct(cellnums(valloop)).*tvalues(2);
    ndiscvalue_cell_max = cover4val(:,valloop).*hard_pct(cellnums(valloop)).*tvalues(3);  
    
    ndiscvalue_min_allcells = ndiscvalue_min_allcells + ndiscvalue_cell_min;
    ndiscvalue_allcells = ndiscvalue_allcells + ndiscvalue_cell;
    ndiscvalue_max_allcells = ndiscvalue_max_allcells + ndiscvalue_cell_max;   
end

% Reshape value to extract just annual values.
nyears = length(DateVec)/12;

% Calculate min, max and mean discounted and non-discounted values
% Minimum discounted value
disc_val_min_mat1 = reshape(discvalue_min_allcells,12,nyears);      % discounted (matrix, rows=months)
disc_val_min_mat = disc_val_min_mat1';              % discounted (matrix, rows = years)
disc_val_min = disc_val_min_mat(:,1);                       % discounted (by year)
% Mean discounted value
disc_val_mat1 = reshape(discvalue_allcells,12,nyears);      % discounted (matrix, rows=months)
disc_val_mat = disc_val_mat1';              % discounted (matrix, rows = years)
disc_val = disc_val_mat(:,1);                       % discounted (by year)
% Maximum discounted value
disc_val_max_mat1 = reshape(discvalue_max_allcells,12,nyears);      % discounted (matrix, rows=months)
disc_val_max_mat = disc_val_max_mat1';              % discounted (matrix, rows = years)
disc_val_max = disc_val_max_mat(:,1);                       % discounted (by year)

% Minimum discounted value
ndisc_val_min_mat1 = reshape(ndiscvalue_min_allcells,12,nyears);      % discounted (matrix, rows=months)
ndisc_val_min_mat = ndisc_val_min_mat1';              % discounted (matrix, rows = years)
ndisc_val_min = ndisc_val_min_mat(:,1);                       % discounted (by year)
% Mean non-discounted value
ndisc_val_mat1 = reshape(ndiscvalue_allcells,12,nyears);      % discounted (matrix, rows=months)
ndisc_val_mat = ndisc_val_mat1';              % discounted (matrix, rows = years)
ndisc_val = ndisc_val_mat(:,1);                       % discounted (by year)
% Maximum discounted value
ndisc_val_max_mat1 = reshape(ndiscvalue_max_allcells,12,nyears);      % discounted (matrix, rows=months)
ndisc_val_max_mat = ndisc_val_max_mat1';              % discounted (matrix, rows = years)
ndisc_val_max = ndisc_val_max_mat(:,1);                       % discounted (by year)

disc_vals = [disc_val_min disc_val disc_val_max];
ndisc_vals = [ndisc_val_min ndisc_val ndisc_val_max];

end
% End of function call





% open temperature baseline files
%
% Assumes files are simply two columns: month and temperature (no headers)

fid = fopen('R1 hi5b.txt');
data = [data; fscanf(fid, '%*s %f')];
fclose(fid);

temps = baseline;

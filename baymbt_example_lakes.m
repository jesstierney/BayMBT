% The following code provides an example of how to use the BayMBT lakes
% calibration to calibrate MBT5Me data. The example data is from:
%
% Miller et al. (2018) A 900-year New England temperature reconstruction
% from in situ seasonally produced branched glycerol dialkyl glycerol
% tetraethers (brGDGTs), Climate of the Past 
% https://doi.org/10.5194/cp-14-1653-2018
%
% These data also appear in Figure 9 of the BayMBT paper:
%
% Martínez-Sosa, Pablo, Tierney, J. E., Stefanescu, I. C., Dearing
% Crampton-Flood, E., Shuman, B. N. & Routson, C (2021). BayMBT: A global
% Bayesian temperature calibration for lacustrine brGDGTs, Geochimica et
% Cosmochimica Acta, in revision.

% read in data
lakeData = readtable('BasinPondMiller.csv');
% pull out mbt5me for ease
mbt5me = lakeData.MBT5Me;
% set prior mean, standard deviation, and model type
prior_mean = 10;
prior_std = 10;
model = "T0";

% calibrate both uncorrected mbt and corrected mbt:
calibratedData = baymbt_predict(mbt5me,prior_mean,prior_std,model,"lake");

% plot the median values and 1-sigma confidence intervals
figure(1); clf;
p1=plot(lakeData.Year,calibratedData.T(:,2)); hold on;
%lower 1-sigma
p2=plot(lakeData.Year,calibratedData.T(:,2)-(calibratedData.T(:,2)-calibratedData.T(:,1))/2);
%upper 1-sigma
p3=plot(lakeData.Year,calibratedData.T(:,2)+(calibratedData.T(:,3)-calibratedData.T(:,2))/2);
legend([p1 p2 p3],'Median','Lower 1-sigma','Upper 1-sigma');
xlabel('Year');
ylabel('MAT above zero');
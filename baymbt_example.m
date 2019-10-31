% The following code provides an example of how to read in fractional
% brGDGT abundances and apply BayMBT. The example data is from:
%
% Crampton-Flood, E. D., Peterse, F., Munsterman, D., & Damsté, J. S. S.
% (2018). Using tetraether lipids archived in North Sea Basin sediments to
% extract North Western European Pliocene continental air temperatures.
% Earth and Planetary Science Letters, 490, 193-205.
%
% These are Pliocene GDGT data from the Hank core. Crampton-Flood describe
% a way to correct brGDGT data for marine contributions, an example of that
% correction is also given here.
%
% These data also appear in Figure 9 of the BayMBT paper:
%
% Dearing Crampton-Flood, E., Tierney, J. E., Peterse, F., Kirkels, F. M.,
% & Sinninghe Damsté, J. S. (2020). BayMBT: A Bayesian calibration model
% for branched glycerol dialkyl glycerol tetraethers in soils and peats.
% Geochimica et Cosmochimica Acta, 268, 142-159.

% read in data
hank = readtable('Hank_core.txt');

% calculate mbt5me
mbt5me = (hank.Ia + hank.Ib + hank.Ic) ./ (hank.Ia + hank.Ib + hank.Ic ...
    + hank.IIa5 + hank.IIb5 + hank.IIc5 + hank.IIIa5 + hank.IIIb5 + hank.IIIc5);

% adjust for marine contributions according to Crampton-Flood et al., (2018)
bwt = 14.5;
mbtmarine = (bwt + 23.7)./59.5;
mbtterr = (mbt5me - mbtmarine.*(1 - hank.percentTerrestrial./100)) ./ ...
    (hank.percentTerrestrial./100);

% set prior mean, standard deviation, and model type
prior_mean = 10;
prior_std = 15;
model = "T0";

% calibrate both uncorrected mbt and corrected mbt:
uncorrected = baymbt_predict(mbt5me,prior_mean,prior_std,model);
corrected = baymbt_predict(mbtterr,prior_mean,prior_std,model);

% plot the median values
figure(1); clf;
p1=plot(hank.Depth,uncorrected.T(:,2)); hold on;
p2=plot(hank.Depth,corrected.T(:,2));
legend([p1 p2],'Uncorrected','Corrected for marine contribution');
xlabel('Depth (m)');
ylabel('MAT above zero');
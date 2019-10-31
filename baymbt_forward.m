function mbt5me = baymbt_forward(t,model)
% function output = baymbt_forward(t,model)
%
% BAYMBT forward model for MBT5Me measured in soils and peats.
% Predicts MBT5Me based on mean annual air temperature or mean temperatures
% above zero.
% ----- Inputs -----
% t: A scalar or vector of temperature values (1 x N) or (N x 1)
%
% model: A string corresponding the model you want to use. Options are:
% "T" = assumes mean annual air temperature (BayMBT)
% "T0" = assumes mean annual temperatures above zero (BayMBT0)
%
% ----- Outputs -----
%
% mbt5me: A 1000-member ensemble of mbt5me values (N x 1000)
%
% ----- Citation -----
% Please cite the following paper when using BayMBT:
%
% Dearing Crampton-Flood, E., Tierney, J. E., Peterse, F., Kirkels, F. M.,
% & Sinninghe Damsté, J. S. (2020). BayMBT: A Bayesian calibration model
% for branched glycerol dialkyl glycerol tetraethers in soils and peats.
% Geochimica et Cosmochimica Acta, 268, 142-159.

% Ensure column vector
    t=t(:);
    % load appropriate model
    if  strcmp(model,"T")
        load('baymbt_params.mat');
    elseif strcmp(model,"T0")
        load('baymbt0_params.mat');
    else
        error('Model not recognized - choices are "T" and "T0"');
    end
    %
    %define dimensions
    Nobs=length(t);
    % parse parameters
    alpha = b_draws_final(:,2);
    betaT = b_draws_final(:,1);
    sigma = sqrt(tau2_draws_final);

    %vectorized calculation of ensemble.
    mbt5me = normrnd(alpha' + t * betaT', repmat(sigma',Nobs,1));
    %any mbt values outside 0 to 1 are forced to be in that range.
    mbt5me(mbt5me>1)=1;
    mbt5me(mbt5me<0)=0;
end
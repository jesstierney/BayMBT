function output = baymbt_predict(mbt5me,prior_mean,prior_std,model)
% function output = baymbt_predict(mbt5me,prior_mean,prior_std,model)
%
% BAYMBT prediction model for MBT5Me measured in soils and peats.
% Predicts Mean annual air temperature or mean temperatures above zero.
% ----- Inputs -----
% mbt5me: A scalar or vector of mbt5me values (1 x N) or (N x 1)
%
% prior_mean: A scalar prior mean value of T in degrees C.
%
% prior_std: A scalar prior standard deviation value of T in degrees C.
%
% model: A string corresponding the model you want to use. Options are:
% "T" = Calculate mean annual air temperature (BayMBT)
% "T0" = Calculate mean annual temperatures above zero (BayMBT0)
%
% ----- Outputs -----
%
% output is a structure containing:
% output.prior_mean: User choice of prior mean
% output.prior_std: User choice of prior standard deviation
% output.T: 2.5%, 50%, and 97.5% confidence intervals of posterior SST
% output.ens: full ensemble of posterior SST (N x 2000);
%
% ----- Citation -----
% Please cite the following paper when using the BayMBT calibrations:
%
% Dearing Crampton-Flood, E., Tierney, J. E., Peterse, F., Kirkels, F. M.,
% & Sinninghe Damsté, J. S. (2020). BayMBT: A Bayesian calibration model
% for branched glycerol dialkyl glycerol tetraethers in soils and peats.
% Geochimica et Cosmochimica Acta, 268, 142-159.

    % Ensure column vector.
    mbt5me=mbt5me(:);
    % load appropriate model
    if  strcmp(model,"T")
        load('baymbt_params.mat');
    elseif strcmp(model,"T0")
        load('baymbt0_params.mat');
    else
        error('Model not recognized - choices are "T" and "T0"');
    end
    % get dimensions of time series and draws
    nd = length(mbt5me);
    n_draws = length(tau2_draws_final);
    
    % parse parameters
    alpha = b_draws_final(:,2);
    betaT = b_draws_final(:,1);
    sigma = sqrt(tau2_draws_final);
    
    % Prior mean and inverse covariance matrix
    pmu = repmat(ones(nd, 1) * prior_mean,1,n_draws);
    pinv_cov = repmat(prior_std,nd,n_draws).^-2;
    
    % Posterior calculations
    post_mean_num = pinv_cov .* pmu + repmat(sigma',nd,1).^-2 .* repmat(betaT',nd,1) .* (mbt5me - repmat(alpha',nd,1));
    post_mean_den = pinv_cov + repmat(betaT',nd,1).^2 .* repmat(sigma',nd,1).^-2;
    post_mean = post_mean_num ./ post_mean_den;
    post_sig = sqrt(post_mean_den.^-1);
    output.ens = post_mean + randn(nd,n_draws).*post_sig;
    output.prior_mean = prior_mean;
    output.prior_std = prior_std;
    if strcmp(model,"T0")
    % if using BayMBT0, truncate at T < 0
    output.ens(output.ens < 0) = NaN;
    else
    end
    output.T = prctile(sort(output.ens,2),[2.5 50 97.5],2);
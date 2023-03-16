function output = baymbt_predict(mbt5me,prior_mean,prior_std,Type,Tmodel)
% function output = baymbt_predict(mbt5me,prior_mean,prior_std,Type,Tmodel)
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
% Type: A string corresponding to the data type. Options are:
% "soil" = use the soil calibration
% "lake" = use the lake calibration
%
% Tmodel: A string corresponding the model you want to use. Options are:
% "T0" = (default) Calculate mean annual temperatures above zero (BayMBT0,
% only option for lakes) 
% "T" = Calculate mean annual air temperature (BayMBT, for soils only)
%
% ----- Outputs -----
%
% output is a structure containing:
% output.prior_mean: User choice of prior mean
% output.prior_std: User choice of prior standard deviation
% output.T: 2.5%, 50%, and 97.5% confidence intervals of posterior SST
% output.ens: full ensemble of posterior SST (N x 1000);
%
% ----- Citation -----
% Please cite the following papers when using the BayMBT calibrations:
%
% For soils:
% Dearing Crampton-Flood, E., Tierney, J. E., Peterse, F., Kirkels, F. M.,
% & Sinninghe Damsté, J. S. (2020). BayMBT: A Bayesian calibration model
% for branched glycerol dialkyl glycerol tetraethers in soils and peats.
% Geochimica et Cosmochimica Acta, 268, 142-159. https://doi.org/10.1016/j.gca.2019.09.043
%
% For lakes:
% Martínez-Sosa, P., Tierney, J. E., Stefanescu, I. C., Dearing
% Crampton-Flood, E., Shuman, B. N., Routson, C. (2021) A global Bayesian
% temperature calibration for lacustrine brGDGTs. Geochimica et
% Cosmochimica Acta, 305, 87-10. https://doi.org/10.1016/j.gca.2021.04.038
%

    % Ensure column vector.
    mbt5me=mbt5me(:);
    % assign Tmodel if not given
    if ~exist('Tmodel','var')
        Tmodel = "T0";
    else
    end
    % check if lake is used w/ T
    if strcmp(Type,"lake") && strcmp(Tmodel,"T")
        error('Model T is not supported for lakes. Use T0 instead')
    else
    end
    % load appropriate model
    if strcmp(Type,"lake")
        load('baymbt0_params_lake.mat','b_draws_final','tau2_draws_final');
    elseif strcmp(Type,"soil")
        if  strcmp(Tmodel,"T")
            load('baymbt_params_soil.mat','b_draws_final','tau2_draws_final');
        elseif strcmp(Tmodel,"T0")
            load('baymbt0_params_soil.mat','b_draws_final','tau2_draws_final');
        else
        error('TModel not recognized - choices are "T" and "T0"');
        end
    else
        error('Type not recognized - choices are "soil" and "lake"');
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
    if strcmp(Tmodel,"T0")
    % if using BayMBT0, truncate at T < 0
    output.ens(output.ens < 0) = NaN;
    else
    end
    output.T = prctile(sort(output.ens,2),[2.5 50 97.5],2);
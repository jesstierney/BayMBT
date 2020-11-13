function mbt5me = baymbt_forward(t,Tmodel,Type)
% function output = baymbt_forward(t,model)
%
% BAYMBT forward model for MBT5Me measured in soils and peats.
% Predicts MBT5Me based on mean annual air temperature or mean temperatures
% above zero.
% ----- Inputs -----
% t: A scalar or vector of temperature values (1 x N) or (N x 1)
%
% Tmodel: A string corresponding the model you want to use. Options are:
% "T" = assumes mean annual air temperature (BayMBT)
% "T0" = assumes mean annual temperatures above zero (BayMBT0)
%
% Type: A string corresponding to the data type. Options are:
% "soil" = use the soil calibration
% "lake" = use the lake calibration
%
% ----- Outputs -----
%
% mbt5me: A 1000-member ensemble of mbt5me values (N x 1000)
%
% ----- Citation -----
% Please cite the following papers when using the BayMBT calibrations:
%
% For soils:
% Dearing Crampton-Flood, E., Tierney, J. E., Peterse, F., Kirkels, F. M.,
% & Sinninghe Damsté, J. S. (2020). BayMBT: A Bayesian calibration model
% for branched glycerol dialkyl glycerol tetraethers in soils and peats.
% Geochimica et Cosmochimica Acta, 268, 142-159.
%
% For lakes:
% Martínez-Sosa, P., Tierney, J. E., Stefanescu, I. C., Dearing
% Crampton-Flood, E., Shuman, B. N., Routson, C. A global Bayesian
% temperature calibration for lacustrine brGDGTs. EarthArXiV, 
% doi: 10.31223/X5PS3P
%

% Ensure column vector
    t=t(:);
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
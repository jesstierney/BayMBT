function baymbt_model(Tmodel,Type)
% function baymbt_model(model)
%
% BAYMBT calibration model for MBT5Me measured in soils and peats.
% Computes slope, intercept, error variance terms using Bayesian linear
% regression.
%
% ----- Inputs -----
% Tmodel: A string corresponding the temperature model you want to calculate. Options are:
% "T" = Calculate mean annual air temperature (BayMBT)
% "T0" = Calculate mean annual temperatures above zero (BayMBT0)
%
% Type: A string corresponding to the data type. Options are:
% "soil" = use the soil calibration data
% "lake" = use the lake calibration data
%
% ----- Citation -----
% This BayMBT calibration method is published here:
%
% Dearing Crampton-Flood, E., Tierney, J. E., Peterse, F., Kirkels, F. M.,
% & Sinninghe Damsté, J. S. (2020). BayMBT: A Bayesian calibration model
% for branched glycerol dialkyl glycerol tetraethers in soils and peats.
% Geochimica et Cosmochimica Acta, 268, 142-159.

if strcmp(Type,"soil")
    calibs = readtable('soil_calibration_data.csv');
elseif strcmp(Type,"lake")
    calibs = readtable('lake_calibration_data.csv');
else
    error('Type not recognized - choices are "soil" and "lake"');
end


%% set up for Bayesian regression of form, Y=X*B + e
if  strcmp(Tmodel,"T")
    X = [calibs.MAT ones(size(calibs,1),1)];
elseif strcmp(Tmodel,"T0")
    X = [calibs.MAT_0 ones(size(calibs,1),1)];
else
    error('TModel not recognized - choices are "T" and "T0"');
end

Y = calibs.MBT5Me;

if strcmp(Type,"lake")
    indUse = calibs.Outliers==0;
    X = X(indUse,:);
    Y = Y(indUse,:);
else
end


iObs=size(X,1);
iR=size(X,2);

%least squares approximation of coefficients B (Bhat)
Bhat = (X'*X)\(X'*Y);
tau2hat=1/(iObs-iR)*(Y-X*Bhat)'*(Y-X*Bhat);

%% priors:
%we use a multivariate normal prior for B
m_p=Bhat; %prior on the mean, let's use least squares values

%prior covariance matrix
var_p=(2.*Bhat).^2; %variance, choose a large number, scale by Bhat
cov_p=diag(var_p); %variance times identity matrix

%we use an inverse gamma prior for tau^2. For your own applications,
%you might have to do some experimentation to find good starting values.
a=.05;
b=.006;

%% GIBBS SAMPLER

% set number of chains to do, and draws for each chain:
Nc=5; %number of chains
Nd=1000; %number of draws
b_draws=NaN(Nc,Nd,iR);
tau2_draws=NaN(Nc,Nd,1);
N=length(X);

for i=1:Nc
%start each chain by setting the initial tau^2 value.
tau2_val=tau2hat*10.*rand(1,1);

%now start the draws
for kk=1:Nd
   
   %calculate conditional posterior distribution of b:
   b_covar=(((X'*X)/tau2_val + inv(cov_p)))^(-1); %covariance
   b_mn=(cov_p\m_p+X'*Y/tau2_val)'*b_covar; %mean
   %make a draw from the distribution
   b_val=mvnrnd(b_mn,b_covar)';
   
   %now calculate marginal posterior distribution of tau^2:
   p_a = a+N/2; %shape parameter
   p_b=b+(1/2)*(Y-X*b_val)'*(Y-X*b_val); %scale parameter
   %make a draw from the distribution
   tau2_val= 1/gamrnd(p_a,1/p_b);
   %put in a ceiling in case this blows up (it can happen sometimes)
   tau2_val= min(20,tau2_val);

   %fill in the vectors:
    b_draws(i,kk,:)=b_val;
    tau2_draws(i,kk)=tau2_val;  
end

end

%% discard first 200 samples of each chain as burn-in:
b_draws=b_draws(:,201:end,:);
tau2_draws=tau2_draws(:,201:end);

% calculate the Rhat statistic to assess convergence.
m=size(b_draws,1);
n=size(b_draws,2);
Bs_var=n/(m-1)*sum((mean(b_draws(:,:,1),2)-mean(mean(b_draws(:,:,1),2))).^2);
Ws_var=mean(var(b_draws(:,:,1),0,2));
Rhat=sqrt(((n-1)/n*Ws_var+1/n*Bs_var)/Ws_var)
%Rhat should be close to 1 if convergence is achieved.

%% plot posteriors with priors.
figure(5); clf; hold on;

set(gcf,'pos',[100 500 900 230]);
subplot(1,3,1)
slope=b_draws(:,:,1);
xt=[min(slope(:))-.01:.0005:max(slope(:))+.01];
post=ksdensity(slope(:),xt); hold on;
prior=normpdf(xt,m_p(1),sqrt(mean(cov_p(1,:))));
p1=plot(xt,prior,'k','linewidth',1.5); hold on;
p2=plot(xt,post,'r','linewidth',1.5);
set(gca,'box','on','linewidth',.75);
legend([p1 p2],'Prior','Posterior');
legend('boxoff');
title('Slope');
ylabel('Prob. density');

subplot(1,3,2)
intercept=b_draws(:,:,2);
xt=[min(intercept(:)-.1):.001:max(intercept(:)+.1)];
post=ksdensity(intercept(:),xt); hold on;
prior=normpdf(xt,m_p(2),sqrt(mean(cov_p(2,:))));
p1=plot(xt,prior,'k','linewidth',1.5); hold on;
p2=plot(xt,post,'r','linewidth',1.5);
set(gca,'box','on','linewidth',.75);
legend([p1 p2],'Prior','Posterior');
legend('boxoff');
title('Intercept');

subplot(1,3,3)
xt=[0:.0001:.02];
post=ksdensity(tau2_draws(:),xt);
prior = b^a/gamma(a)*xt.^(-a-1).*exp(-b./xt);
p1=plot(xt,prior,'k','linewidth',1.5); hold on;
p2=plot(xt,post,'r','linewidth',1.5);
set(gca,'box','on','linewidth',.75);
legend([p1 p2],'Prior','Posterior');
legend('boxoff');
title('Error variance (\tau^2)');

%% reshape and thin draws, save parameters
b_draws_r=reshape(b_draws,size(b_draws,1)*size(b_draws,2),size(b_draws,3));
tau2_draws_r=reshape(tau2_draws,size(tau2_draws,1)*size(tau2_draws,2),1);
b_draws_final=b_draws_r(1:4:end,:);
tau2_draws_final=tau2_draws_r(1:4:end);

if  strcmp(Tmodel,"T") && strcmp(Type,"soil")
    filename = 'baymbt_params_soil.mat';
elseif strcmp(Tmodel,"T0") && strcmp(Type,"soil")
    filename = 'baymbt0_params_soil.mat';
elseif strcmp(Tmodel,"T") && strcmp(Type,"lake")
    filename = 'baymbt_params_lake.mat';
elseif strcmp(Tmodel,"T0") && strcmp(Type,"lake")
    filename = 'baymbt0_params_lake.mat';
else
end

save(filename,'b_draws_final','tau2_draws_final');

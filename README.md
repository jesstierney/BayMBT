# BayMBT

BayMBT is a Bayesian calibration for the branched GDGT MBT5Me proxy in soils, peats, and lakes. In this package you'll find the calibration data, the code to construct the calibration, code to forward model MBT5Me from temperatures, and code to predict temperatures from MBT5Me. Please cite the following publication when using BayMBT for soils and lakes:

Dearing Crampton-Flood, E., Tierney, J. E., Peterse, F., Kirkels, F. M.,
& Sinninghe Damsté, J. S. (2020). BayMBT: A Bayesian calibration model
for branched glycerol dialkyl glycerol tetraethers in soils and peats.
Geochimica et Cosmochimica Acta, 268, 142-159. <https://doi.org/10.1016/j.gca.2019.09.043>

And the following publication when using BayMBT for lakes:

Martínez-Sosa, P., Tierney, J. E., Stefanescu, I. C., Dearing Crampton-Flood, E., Shuman, B. N., Routson, C. A global Bayesian temperature calibration for lacustrine brGDGTs. EarthArXiV <https://doi.org/10.31223/X5PS3P>

## Description of contents

**baymbt_model.m** This contains the code to construct the calibration model. It loads either the soil or lake core top dataset and calculates model parameters using Bayesian linear regression. For everyday use you don't need to run this as the parameters are already calculated and provided here as .mat files (_params_.mat files). This is provided here for transparency.

**baymbt_forward.m** This calculates MBT5Me values from user-provided temperatures, using the parameters in the params.mat files. There are no dependent functions except for the .mat files.

**baymbt_predict.m** This calculates temperatures from user-provided MBT5Me values, using the parameters in the params.mat files, through Bayesian inference. It requires a prior mean and standard deviation on temperature. For soils, one can use either the "BayMBT" model (which is calibrated to mean annual temperatures) or the "BayMBT0" model (which is calibrated to mean annual temperatures above 0C). For lakes the default is the "BayMBT0" model. For both soils and lakes, we recommend using BayMBT0 for most applications, as it is a better fit to the core top data. The implicit assumption in BayMBT0 is that brGDGTs are not synthesized at temperatures below freezing.

**baymbt_example_soils.m** This is an example of how to use **baymbt_predict.m** to predict temperatures from soil brGDGT data. It loads fractional GDGT abundance data from the Dearing Crampton-Flood et al. Hank Core study (Hank_core.txt) <https://doi.org/10.1016/j.epsl.2018.03.030>, calculates MBT5Me, and calibrates it to temperature. These data are shown in Figure 9 of Dearing Crampton-Flood et al., (2020).

**baymbt_example_lakes.m** This is an example of how to use **baymbt_predict.m** to predict temperatures from lake brGDGT data. It loads data from the Miller et al. Basin Pond study (BasinPondMiller.csv) <https://doi.org/10.5194/cp-14-1653-2018> and calibrates MBT5Me to temperature. These data are shown in Figure 9 of Martínez-Sosa et al., (2021).

To cite the code repository directly use:

*Tierney, Jessica E., 2021. BayMBT. <https://github.com/jesstierney/BayMBT>.*

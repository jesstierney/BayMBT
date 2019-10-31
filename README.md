# BayMBT

BayMBT is a Bayesian calibration for the branched GDGT MBT5Me proxy in soils and peats. In this package you'll find the calibration data, the code to construct the calibration, code to forward model MBT5Me from temperatures, and code to predict temperatures from MBT5Me. Please cite the following publication when using BayMBT:

Dearing Crampton-Flood, E., Tierney, J. E., Peterse, F., Kirkels, F. M.,
& Sinninghe Damsté, J. S. (2020). BayMBT: A Bayesian calibration model
for branched glycerol dialkyl glycerol tetraethers in soils and peats.
Geochimica et Cosmochimica Acta, 268, 142-159.

## Description of contents

**baymbt_model.m** This contains the code to construct the calibration model. It loads the core top dataset in calibration_data.csv and calculates model parameters using Bayesian linear regression. For everyday use you don't need to run this as the parameters are already calculated and provided here as .mat files (baymbt_params.mat and baymbt0_params.mat). This is provided here for transparency.

**baymbt_forward.m** This calculates MBT5Me values from user-provided temperatures, using the parameters in the params.mat files. There are no dependent functions except for the .mat files.

**baymbt_predict.m** This calculates temperatures from user-provided MBT5Me values, using the parameters in the params.mat files, through Bayesian inference. It requires a prior mean and standard deviation on temperature. One can use either the "BayMBT" model (which is calibrated to mean annual temperatures) or the "BayMBT0" model (which is calibrated to mean annual temperatures above 0C). We recommend using the latter for most applications, as it is a better fit to the core top data. The implicit assumption in BayMBT0 is that brGDGTs are not synthesized at temperatures below freezing.

**baymbt_example.m** This is an example of how to use **baymbt_predict.m**. It loads fractional GDGT abundance data from Hank_core.txt, calculates MBT5Me, and calibrates it to temperature. These data are shown in Figure 9 of Dearing Crampton-Flood et al., (2020).

To cite the code repository directly use:

*Tierney, Jessica E., 2019. BayMBT. \<https://github.com/jesstierney/BayMBT \>.*
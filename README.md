# LossCalculator

A python script to calculate the loss in amorphous oxides.  The code currently supports the use of bootstrapping, importance sampling, 
and the calculations of both mechanical and dielectric loss in amorphous oxides.  It can be run using

python CalcLoss.py [suffix]

If no suffix is specified, it reads Data/tls.tot.  If a suffix is specified, it reads Data/tls.tot_[suffix].

The functions shouldn't need to be changed, but there are several options within CalcLoss.py.

plotLoss = True/False, whether the python code plots the calculated loss

bSamples = [integer], the number of bootstrap samples to use

percentiles = True/False, whether or not to use the percentile bootstrap method

saveLoss = [string], filename in which the loss and, if calculated, the 95% confidence intervals are saved

saveDistros = True/False, whether or not to save the loss value calculated at every temperature

BoltzmannCorrection = True/False, whether or not to use importance sampling to correct for the Boltzmann factor bias in classical MD

constG = True/False, whether to use the average coupling constant in calculations

constT = True/False, whether to use the average relaxation time in calculations

constY = True/False, whether to use the average modulus in calculations

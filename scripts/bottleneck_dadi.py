#!/usr/bin/env python3

"""
Matt Gibson
Aug. 2019
Indiana University Departement of Biology
Moyle Lab

Dadi functions and code for performing demographic model inference and hypothesis testing.
This code takes in a SFS file, fits both the bottleneck and neutral models, and then performs
a nonparametric bootstrap on the bottleneck model. 


Arg 1: input sfs file (dadi format)
Arg 2: output file name (for bootstraps)
Arg 3: whether or not to perform bootstrapping ('boot' if yes)
"""


import dadi
from dadi import Godambe
from dadi import Inference
import matplotlib
import sys
import numpy, pylab

####Models#################
###########################
def bottleneckF(params, ns, pts):
    """
    Simple instantaneous bottlneck model with inbreeding
    """
	nuB, nuF, TB, TF, F = params
	xx = dadi.Numerics.default_grid(pts)

	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.Integration.one_pop(phi, xx, TB, nuB)
	phi = dadi.Integration.one_pop(phi, xx, TF, nuF)

	fs = dadi.Spectrum.from_phi_inbreeding(phi, ns, (xx,), (F,), (2,))
	return(fs)

def constant(notused, ns, pts):
    """
    Standard neutral model.

    ns = (n1,)

    n1: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs
###########################
###########################


#Read in SFS
pop = dadi.Spectrum.from_file(sys.argv[1])


##Tajima's D##
print("TAJIMA:")
print(pop.Tajima_D())

ns = pop.sample_sizes




####Bottleneck Fit#########
###########################

print("BOTTLENECK MODEL")

##Define starting parameters and bounds
params = [1,1,1,1,0.5]
lowerbounds = [1e-2, 1e-2, 0, 0,0]
upperbounds = [100, 100, 10, 20,1]
pts_l = [30]

func_ex = dadi.Numerics.make_extrap_func(bottleneckF)

model = func_ex(params, ns, pts_l)

ll_model = dadi.Inference.ll_multinom(model, pop)

print('Model log-likelihood:', str(ll_model))
theta0 = dadi.Inference.optimal_sfs_scaling(model, pop)
print('Theta0_1:', str(theta0))

p0 = dadi.Misc.perturb_params(params,
	fold=1,
	upper_bound=upperbounds,
	lower_bound=lowerbounds)

popt = dadi.Inference.optimize_log_fmin(p0, pop, func_ex, pts_l,
	lower_bound=lowerbounds,
	upper_bound=upperbounds,
	verbose=len(params),
	maxiter=50)

print('Optimized parameters', str(repr(popt)))

model = func_ex(popt[0], ns, pts_l)
ll_opt = dadi.Inference.ll_multinom(model, pop)

print('Optimized log-likelihood:', str(ll_opt))
print('Theta0_2:', str(theta0))



##Print out the scaled SFS and the Anscombe Poisson residuals
print("==============================================")

print(Inference.optimally_scaled_sfs(model, pop))
rescaled = Inference.optimally_scaled_sfs(model, pop)
print(Inference.Anscombe_Poisson_residual(rescaled,pop))
dadi.Plotting.plot_1d_comp_multinom(model,pop)

print("==============================================")




####Neutral Fit#########
###########################

print("NEUTRAL MODEL")

params = [1]
lowerbounds = [100]
upperbounds = [0.001]
pts_l = [30]

func_ex = dadi.Numerics.make_extrap_func(constant)

model = func_ex(params, ns, pts_l)

ll_model = dadi.Inference.ll_multinom(model, pop)

print('Model log-likelihood:', str(ll_model))
theta0 = dadi.Inference.optimal_sfs_scaling(model, pop)
print('Theta0_1:', str(theta0))


p0 = dadi.Misc.perturb_params(params,
	fold=1,
	upper_bound=upperbounds,
	lower_bound=lowerbounds)

popt = dadi.Inference.optimize_log_fmin(p0, pop, func_ex, pts_l,
	lower_bound=  lowerbounds,
	upper_bound=upperbounds,
	verbose=len(params),
	maxiter=50)


print('Optimized parameters', str(repr(popt)))

model = func_ex(popt[0], ns, pts_l)
ll_opt = dadi.Inference.ll_multinom(model, pop)

print('Optimized log-likelihood:', str(ll_opt))

print('Theta0_2:', str(theta0))

print("==============================================")

print(pop)
print(Inference.optimally_scaled_sfs(model, pop))
rescaled = Inference.optimally_scaled_sfs(model, pop)
print(Inference.Anscombe_Poisson_residual(rescaled,pop))
dadi.Plotting.plot_1d_comp_multinom(model,pop)

print("==============================================")




####Bootstrapping##########
###########################

out1 = open(sys.argv[2], "w", buffering = 1)

#Redefine bottleneck model

params = [1,1,1,1,0.5]
lowerbounds = [1e-2, 1e-2, 0, 0,0]
upperbounds = [100, 100, 10, 20,1]
pts_l = [30]

func_ex = dadi.Numerics.make_extrap_func(bottleneckF)

if (sys.argv[3] == 'boot'):
    nboot = 2000
    bootstrapped_data = []
    bootstrap_results = []
    for i in range(0,nboot):
            s = pop.sample()

            #Perturb parameters a bit each time
            p0 = dadi.Misc.perturb_params(params,
            fold=1,
            upper_bound=upperbounds,
            lower_bound=lowerbounds)

            popt = dadi.Inference.optimize_log_fmin(p0, s, func_ex, pts_l,
            lower_bound=lowerbounds,
            upper_bound=upperbounds,
            verbose=len(params),
            maxiter=50)

            new = []
            for p in popt[0]:
                    new.append(str(p))

            print(new)
            print(list(popt))

            out1.write('\t'.join(new) + '\t'+ str(popt[1]) + '\n')

    out1.close()



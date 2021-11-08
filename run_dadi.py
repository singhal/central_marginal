import dadi
import pylab
import Optimize_Functions
import re
import numpy
import os
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='run dadi')
parser.add_argument('--file', help="file to run this on")
parser.add_argument('--folded', action='store_true', help="should it be folded?")
args = parser.parse_args()

filename = args.file
folded = args.folded

# decide projection based on mean of data
d = pd.read_csv(filename, sep = "\t")
d['total'] = d['species'] + d['species.1']
projection = int(round(d.total.mean() / 2, 0) * 2)

dd = dadi.Misc.make_data_dict(filename)
# if it is folded, then it is NOT polarized
fs = dadi.Spectrum.from_data_dict(dd, pop_ids = [ 'species' ], 
    projections = [projection], polarized = not folded)
ns = [projection]

# set sampling grid dynamically
rounded = round(projection, -1) + 10
pts = [rounded, rounded + 10, rounded + 20]

sp = re.sub("^.*/", "", filename)
sp = re.sub(".txt", "", sp)
if not os.path.isdir(sp):
    os.mkdir(sp)

#############
# can change this
# but good for now
#############
rounds = 3
reps = [10, 20, 50]
maxiters = [5, 10, 20]
folds = [3, 2, 1]


# for each model
# upper, lower, and initial guess
param_names = "na"
upper_bound = [0]
lower_bound = [1]
p_init = [1]
Optimize_Functions.Optimize_Routine(fs, pts, "%s/neutral" % sp, 
    "neutral", dadi.Demographics1D.snm,
    rounds, len(upper_bound), fs_folded = folded,
    param_labels = param_names, in_upper = upper_bound, 
    in_lower = lower_bound,
    in_params = p_init, reps = reps, 
    maxiters = maxiters, folds = folds)


param_names = "nu, T"
lower_bound = [1e-2 , 0]
upper_bound = [10000 , 3]
p_init = [2, 0.2]
Optimize_Functions.Optimize_Routine(fs, pts, "%s/expgrowth" % sp, 
    "expgrowth", dadi.Demographics1D.growth,
    rounds, len(upper_bound), fs_folded = folded,
    param_labels = param_names, in_upper = upper_bound, 
    in_lower = lower_bound,
    in_params = p_init, reps = reps,
    maxiters = maxiters, folds = folds)


param_names = "nu, T"
lower_bound = [1e-2 , 0]
upper_bound = [10000 , 3]
p_init = [2, 0.2]
Optimize_Functions.Optimize_Routine(fs, pts, "%s/two_epoch" % sp, 
    "two_epoch", dadi.Demographics1D.two_epoch,
    rounds, len(upper_bound), fs_folded = folded,
    param_labels = param_names, in_upper = upper_bound, 
    in_lower = lower_bound,
    in_params = p_init, reps = reps,
    maxiters = maxiters,  folds = folds)


param_names = "nuB, nuF, T"
lower_bound = [1e-2 , 1e-2, 0]
upper_bound = [10000 , 10000, 3]
p_init = [2, 2, 0.2]
Optimize_Functions.Optimize_Routine(fs, pts, "%s/bottlegrowth" % sp, 
    "bottlegrowth", dadi.Demographics1D.bottlegrowth,
    rounds, len(upper_bound), fs_folded = folded,
    param_labels = param_names, in_upper = upper_bound, 
    in_lower = lower_bound,
    in_params = p_init, reps = reps,
    maxiters = maxiters,  folds = folds)
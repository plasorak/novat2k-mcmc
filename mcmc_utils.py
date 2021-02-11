import uproot as up
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='MCMC Utils')
parser.add_argument('input', nargs="+", help="the input file")
parser.add_argument('--burnin', nargs='?', help='number of burnin steps', default=20000, type=int)
argument = parser.parse_args()
standard_burnin = argument.burnin
input_file = argument.input[0]

class Variable:
    def __init__(self, tree_var_name, nice_name, the_range, nbins, change_of_var, unit):
        self.tree_var_name = tree_var_name
        self.nice_name = nice_name
        self.range = the_range
        self.nbins = nbins
        self.change_of_var = change_of_var
        self.unit = unit

def divide_by_2pi(val):
    return val/np.pi

def sinsq2(val):
    return np.power(np.sin(val), 2.)

def noop(val):
    return val
syst_name_nova = [r"CCQE z-exp EV shift \#1",                                    # 0
                  r"CCQE z-exp EV shift \#2",                                    # 1 
                  r"CCQE z-exp EV shift \#3",                                    # 2 
                  r"CCQE z-exp EV shift \#4",                                    # 3 
                  r"MaCCRES",                                                    # 4 
                  r"MvCCRES",                                                    # 5 
                  r"MaNCRES",                                                    # 6 
                  r"MvNCRES",                                                    # 7 
                  r"ZNormCCQE",                                                  # 8 
                  r"RPA shape: higher-$Q^{2}$ enhancement (2020)",               # 9 
                  r"RPA shape: low-$Q^{2}$ suppression (2020)",                  #10 
                  r"RES low-$Q^{2}$ suppression",                                #11 
                  r"DIS $\nu nCC1\pi$",                                          #12 
                  r"hN FSI mean free path",                                      #13 
                  r"hN FSI fate fraction eigenvector",                           #14 
                  r"MEC E$_{\nu}$ shape",                                        #15 
                  r"MEC E$_{\bar{\nu}}$ shape",                                  #16 
                  r"MEC 2020 ($q_{0}$, $|\vec{q}|$) response, neutrinos",        #17 
                  r"MEC 2020 ($q_{0}$, $|\vec{q}|$) response, antineutrinos",    #18 
                  r"MEC initial state np fraction, neutrinos",                   #19 
                  r"MEC initial state np fraction, antineutrinos",               #20 
                  r"Radiative corrections for $\nu_{e}$",                        #21 
                  r"Radiative corrections for $\bar{\nu}_{e}$",                  #22 
                  r"Second class currents",                                      #23 
                  r"Genie PC 0",                                                 #24 
                  r"Genie PC 1",                                                 #25 
                  r"Genie PC 2",                                                 #26 
                  r"Genie PC 3",                                                 #27 
                  r"Genie PC 4",                                                 #28 
                  r"Genie PC 5",                                                 #29 
                  r"Genie PC 6",                                                 #30 
                  r"Genie PC 7",                                                 #31 
                  r"Genie PC 8",                                                 #32 
                  r"Genie PC 9",                                                 #33 
                  r"Genie PC 10",                                                #34 
                  r"Genie PC 11",                                                #35 
                  r"$\nu_{\tau}$ Scale",                                         #36 
                  r"PPFX Flux Component 00",                                     #37 
                  r"PPFX Flux Component 01",                                     #38 
                  r"PPFX Flux Component 02",                                     #39 
                  r"PPFX Flux Component 03",                                     #40 
                  r"PPFX Flux Component 04",                                     #41 
                  r"Absolute Calibration",                                       #42 
                  r"Relative Calibration",                                       #43 
                  r"Calibration Shape",                                          #44 
                  r"Calibration Drift",                                          #45 
                  r"Light Level FD",                                             #46 
                  r"Light Level ND",                                             #47 
                  r"Cherenkov",                                                  #48 
                  r"Uncorr ND Mu Energy Scale 2020",                             #49 
                  r"Uncorr MuCat Mu Energy 2020",                                #50 
                  r"Neutron Pile-up 2020",                                       #51 
                  r"Corr Mu Energy Scale 2020",                                  #52 
                  r"Uncorr FD Mu Energy Scale 2020",                             #53 
                  r"Neutron visible energy systematic 2018",                     #54 
                  r"Acceptance ND to FD Kinematics Signal FHC 2020",             #55 
                  r"Acceptance ND to FD Kinematics Signal RHC 2020",             #56 
                  r"Michel Electrons Tagging Uncertainty"]                       #57 

osc_variables = [Variable("val_dcp" , "$\delta_{CP}$"          , [ 0     , 2.    ],  50, divide_by_2pi, "/ $\pi$"), # 0
                 Variable("val_mh"  , "Mass hierarchy"         , [-1     , 1.    ],   2, noop         , ""),        # 1
                 Variable("val_th13", r"$\sin^{2} \theta_{13}$", [ 0.018 , 0.025 ], 100, sinsq2       , ""),        # 2
                 Variable("val_th13", r"$\theta_{13}$"         , [ 0.12  , 0.17  ], 100, noop         , ""),        # 3
                 Variable("val_th23", r"$\sin^{2} \theta_{23}$", [ 0.4   , 0.65  ], 100, sinsq2       , ""),        # 4
                 Variable("val_th23", r"$\theta_{23}$"         , [ 0.65  , 0.95  ],  50, noop         , ""),        # 5
                 Variable("val_dm32", r"$\left|\Delta m^2_{32}\right|$"      , [ 0.0023, 0.0027],  50, noop         , "eV$^2$")]  # 6

nova_syst_variables = [Variable("val_nova_syst_"+str(i), syst_name_nova[i], [-5.,5.], 100, noop, "") for  i in range(0,58)]
grouped_nova_syst_variables = {"largest cross section": nova_syst_variables[:24],
                               "other cross section"  : nova_syst_variables[24:37],
                               "flux"                 : nova_syst_variables[37:41],
                               "calibration"          : nova_syst_variables[42:46],
                               "other"                : nova_syst_variables[46:]}


def acf(x, lag):
    # Slice the relevant subseries based on the lag
    if len(x)-lag<=0:
        print("x (length: "+str(len(x))+") is smaller than lag ("+str(lag)+").")
        raise
    y1 = x[:(len(x)-lag)]
    y2 = x[lag:]
    # print(lag,len(y1), len(y2))
    # Subtract the mean of the whole series x to calculate Cov
    mean=np.mean(x)
    sum_product = np.sum((y1-mean)*(y2-mean))
    # Normalize with var of whole series
    return sum_product / ((len(x) - lag) * np.var(x))

def autocorrelation(x):
    max_lag=20000

    x = (x-np.mean(x))/np.sqrt(np.var(x))
    n = len(x)
    
    this_max_lag = min(len(x)-1, max_lag)
    autocorr = np.zeros(this_max_lag)

    for i in range(0, n):
        for j in range(0, min(this_max_lag, n-i)):
            autocorr[j]=autocorr[j] + x[i] * x[i+j] / (n-j)
    
    return autocorr



def get_var_array(tree, variable, burnin=None):
    
    if variable.tree_var_name.find("nova_syst") == -1:
        v = tree[variable.tree_var_name].array()
        #print("v is:"+ str(v))
    elif variable.tree_var_name.find("nova_syst") > -1:
        index = int(variable.tree_var_name.split("_")[-1])
        v = tree["val_nova_syst"].array()
        v = np.array(v)
        v = np.transpose(v)
        v = v[index]
        
    if burnin is None:
        return np.array(v[standard_burnin:])
    else:
        return np.array(v[burnin:])

    
def hdi(bin_val, prob):
    initial_shape = bin_val.shape
    bin_val = bin_val.transpose().flatten()
    
    s = np.sum(bin_val)
    if s>1.0001 or s<0.9999:
        print ("histogram hasnt been normalised!!", s)
        raise
    
    max_value = np.max(bin_val)
    min_value = np.min(bin_val)
    n_prob_step = 10000

    prob_values = np.linspace(min_value, max_value, n_prob_step)
    best_step = None
    previous_ci=0
    previous_mask=np.zeros(len(bin_val))
    
    for i_prob_step in range(n_prob_step):
        current_prob = prob_values[i_prob_step]
        this_ci=0
        this_mask=np.zeros(len(bin_val))
        
        for i_val in range(len(bin_val)):
            bin_content = bin_val[i_val]
            
            if bin_content>current_prob:
                this_ci += bin_content
                this_mask[i_val] = 1.
                
        if previous_ci>prob and prob>this_ci:
            return previous_mask.reshape(initial_shape)
        
        previous_ci = this_ci
        previous_mask = this_mask

    return
        
        
def get_2dhisto(tree, varx, vary, apply_change_of_var=True):
    stepx = (varx.range[1] - varx.range[0])/ varx.nbins
    stepy = (vary.range[1] - vary.range[0])/ vary.nbins
    
    bx = np.arange(varx.range[0], varx.range[1]+stepx, stepx)
    by = np.arange(vary.range[0], vary.range[1]+stepy, stepy)
    bx = bx[:varx.nbins+1]
    by = by[:vary.nbins+1]
    
    valuesx = np.array(get_var_array(tree, varx))
    valuesy = np.array(get_var_array(tree, vary))
    
    if apply_change_of_var:
        valuesx = varx.change_of_var(valuesx)
        valuesy = vary.change_of_var(valuesy)

    weights = [1./len(valuesx)] * len(valuesx)
    
    histo, binsx, binsy = np.histogram2d(valuesx, valuesy, bins=[bx,by], weights=weights)
    return histo, binsx, binsy

def get_histo(tree, var, apply_change_of_var=True):
    step = (var.range[1] - var.range[0])/ var.nbins
    b = np.arange(var.range[0], var.range[1]+step, step)

    values = np.array(get_var_array(tree, var))
    if apply_change_of_var:
        values = var.change_of_var(values)

    weights = [1./len(values)] * len(values)
    b = [np.array([-np.inf]), b, np.array([np.inf])]
    b = np.concatenate(b)
    histo, bins = np.histogram(values, bins=b, weights=weights)
    histo[1]  += histo[0]
    histo[-2] += histo[-1]
    histo = histo[1:-1]
    bins = bins[1:-1]
    return histo, bins

import uproot as up
import numpy as np
import argparse

class mcmc_parameters:
    def __init__(self):
        self.burnin=0
        self.nsteps_max=0
        self.nsteps_noburnin=0
        self.nsteps=0
        self.input_trees=[]
        self.do_compare=False
        
params=mcmc_parameters()

class arguments:
    
    def __init__(self, name):
        self.name = name
        self.parser = argparse.ArgumentParser(description='MCMC Utils')
        self.parser.add_argument('--cafana'    ,                      help="Whether input file is cafana/stan format or not", action="store_true")
        self.parser.add_argument('--input'     , nargs="+",           help="Input file, if you are comparing different chain, add --comparison")
        self.parser.add_argument('--burnin'    , nargs='?', type=int, help='Number of burnin steps', default=20000)
        self.parser.add_argument('--nsubfile'  , nargs='?', type=int, help='Number of subfiles before hadding', default=1)
        # self.parser.add_argument("--stride"  , nargs=1,   type=int, help="Stride parameter over the MCMC")ll
        self.parser.add_argument("--nstepmax"  , nargs=1,   type=int, help="Maximum number of steps to consider", default=[-1])
        self.parser.add_argument('--output'    , nargs=1,   type=str, help='Output name (default is some clever combination of the name of the plotter you are running and the name of the file)')
        # self.parser.add_argument('--list'      ,                      help='Whether the input is a list', action='store_true')
        self.parser.add_argument('--comparison',                      help='Whether the input files are to be compared', action='store_true')
        self.parser.add_argument('--label'     , nargs="+",           help='Labels')
        
    def parse(self):
        self.argument    = self.parser.parse_args()
        params.burnin    = self.argument.burnin
        input_file       = self.argument.input
        self.do_compare  = self.argument.comparison
        params.do_compare=self.do_compare

        params.nsteps          = 0
        params.nsteps_noburnin = 0
        params.nsteps_max      = self.argument.nstepmax[0]
        if self.argument.label is None:
            self.label=input_file
        else:
            self.label=self.argument.label
        
        self.output = self.name+"_"
        if type(input_file)==str:
            self.output += input_file.split(r".root")[0]+"_plot.pdf"
            if not self.argument.cafana:
                t=up.open(input_file)["result"]
            else:
                t=up.open(input_file)["run/samples"]
            params.input_trees.append(t)
            if not self.do_compare:
                params.nsteps+=t.num_entries
                params.nsteps_noburnin+=t.num_entries-params.burnin
                
        if type(input_file)==list:
            for i in input_file:
                self.output+=i.split(r".root")[0]+"_"
                if not self.argument.cafana:
                    t=up.open(i)["result"]
                else:
                    t=up.open(i)["run/samples"]
                params.input_trees.append(t)
                if not self.do_compare:
                    params.nsteps+=t.num_entries
                    params.nsteps_noburnin+=t.num_entries-params.burnin
            self.output+="plot.pdf"
            
        if self.argument.output is not None:
            self.output = self.argument.output[0]

        if not self.do_compare:
            print(f"Parsed {len(params.input_trees)} input file(s), with a total of {params.nsteps/1000000.}M entries in them (and {params.nsteps_noburnin/1000000.}M excluding {params.burnin/1000}k burnin steps for each of them).")
        else:
            print(f"Parsed {len(params.input_trees)} input file(s) for comparison")
            count=0
            for t in params.input_trees:
                print(f"file {count+1} ({input_file[count]}) has {t.num_entries/1000000.}M entries (and {(t.num_entries-params.nsteps_noburnin)/1000000.}M excluding {params.burnin/1000}k burnin steps).")
                count+=1
            
        print(f"Output is {self.output}")
        
        
class Variable:
    def __init__(self, tree_var_name, nice_name, the_range, nbins, change_of_var, unit):
        self.tree_var_name = tree_var_name
        self.nice_name = nice_name
        self.range = the_range
        self.nbins = nbins
        self.change_of_var = change_of_var
        self.unit = unit

def divide_by_pi(val):
    return val/np.pi

def sinsq2(val):
    return np.power(np.sin(val), 2.)

def noop(val):
    return val

def mabs(val):
    return abs(val)

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

llh    = Variable("loglikelihood", "log likelihood"            , [ -10 , 0    ], 100, noop, "")
dcp    = Variable("val_dcp" , "$\delta_{CP}$"                  , [ 0    , 2.    ], 50, divide_by_pi, "/ $\pi$")
mh     = Variable("val_mh"  , "Mass hierarchy"                 , [-1    , 1.    ],   2, noop        , ""       )
s2th13 = Variable("val_th13", r"$\sin^{2} \theta_{13}$"        , [ 0.1, 0.2  ], 100, sinsq2      , ""       )
th13   = Variable("val_th13", r"$\theta_{13}$"                 , [ 0.14, 0.16   ], 100, noop, "/ $\pi$")
s2th23 = Variable("val_th23", r"$\sin^{2} \theta_{23}$"        , [ 0.3  , 0.7   ], 100, sinsq2      , ""       )
th23   = Variable("val_th23", r"$\theta_{23}$"                 , [ 0.65 , 0.95  ], 50, noop, "/ $\pi$")
dm32   = Variable("val_dm32", r"$\left|\Delta m^2_{32}\right|$", [ 0.002, 0.0028], 100, noop        , "eV$^2$" )

cafana_llh    = Variable("logprob"  , "log likelihood"                 , [-10  , 0    ], 100, noop, "")
cafana_dcp    = Variable("delta(pi)", "$\delta_{CP}$"                  , [ 0    , 2.    ], 50, noop, "/ $\pi$")
cafana_mh     = Variable("MH"       , "Mass hierarchy"                 , [-1    , 1.    ],   2, noop, ""       )
cafana_th13   = Variable("th13"     , r"$\theta_{13}$"                 , [ 0.14  , 0.16   ], 100, noop, "/ $\pi$")
cafana_th23   = Variable("th23"     , r"$\theta_{23}$"                 , [ 0.65 , 0.95  ], 50, noop, "/ $\pi$")
cafana_dm32   = Variable("dmsq32"   , r"$\left|\Delta m^2_{32}\right|$", [ 0.002, 0.0028], 100, mabs, "eV$^2$" )

nova_syst_variables = [Variable("val_nova_syst_"+str(i), syst_name_nova[i], [-5.,5.], 100, noop, "") for i in range(0,58)]
nova_syst_name = ["ZExpAxialFFSyst2020_EV1",
                  "ZExpAxialFFSyst2020_EV2",
                  "ZExpAxialFFSyst2020_EV3",
                  "ZExpAxialFFSyst2020_EV4",
                  "MaCCRES",
                  "MvCCRES",
                  "MaNCRES",
                  "MvNCRES",
                  "ZNormCCQE",
                  "RPAShapeenh2020",
                  "RPAShapesupp2020",
                  "LowQ2RESSupp2020",
                  "DISvnCC1pi_2020",
                  "hNFSI_MFP_2020",
                  "hNFSI_FateFracEV1_2020",
                  "MECEnuShape2020Nu",
                  "MECEnuShape2020AntiNu",
                  "MECShape2020Nu",
                  "MECShape2020AntiNu",
                  "MECInitStateNPFrac2020Nu",
                  "MECInitStateNPFrac2020AntiNu",
                  "radcorrnue",
                  "radcorrnuebar",
                  "2ndclasscurr",
                  "genie_small_pc00",
                  "genie_small_pc01",
                  "genie_small_pc02",
                  "genie_small_pc03",
                  "genie_small_pc04",
                  "genie_small_pc05",
                  "genie_small_pc06",
                  "genie_small_pc07",
                  "genie_small_pc08",
                  "genie_small_pc09",
                  "genie_small_pc10",
                  "genie_small_pc11",
                  "NuTauScale",
                  "ppfx_hadp_beam_pc00",
                  "ppfx_hadp_beam_pc01",
                  "ppfx_hadp_beam_pc02",
                  "ppfx_hadp_beam_pc03",
                  "ppfx_hadp_beam_pc04",
                  "Calibration",
                  "RelativeCalib",
                  "CalibShape",
                  "CalibDrift",
                  "Light_Level_FD",
                  "Light_Level_ND",
                  "Cherenkov",
                  "UnCorrNDMuEScaleSyst2020",
                  "UnCorrMuCatMuESyst2020",
                  "PileupMuESyst2020",
                  "CorrMuEScaleSyst2020",
                  "UnCorrFDMuEScaleSyst2020",
                  "LeptonAngleSystNDXZ2020",
                  "LeptonAngleSystNDYZ2020",
                  "LeptonAngleSystFDXZ2020",
                  "LeptonAngleSystFDYZ2020",
                  "NeutronEvisPrimariesSyst2018",
                  "NormHornCorr",
                  "NormFHC2020",
                  "NormRHC2020",
                  "accept_signalkin_pTextrap_FHC_2020",
                  "accept_signalkin_pTextrap_RHC_2020",
                  "michel_tagging2020",
                  "RockScale",
                  "cosmicScale"]
cafana_nova_syst_variables = [Variable(nova_syst_name[i], syst_name_nova[i], [-5.,5.], 100, noop, "") for i in range(58)]

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



def get_var_array(variable, ntree=None):
    ret=[]
    trees=[]
    if ntree is None:
        trees=params.input_trees
    else:
        trees.append(params.input_trees[ntree])
        
    if variable.tree_var_name.find("nova_syst") == -1:
        for tree in params.input_trees:
            if params.nsteps_max==-1:
                ret.append(tree[variable.tree_var_name].array()[params.burnin:])
            else:
                r=tree[variable.tree_var_name].array()[params.burnin:params.nsteps_max+params.burnin]
                ret.append(r)
                    
            
    elif variable.tree_var_name.find("nova_syst") > -1:
        index = int(variable.tree_var_name.split("_")[-1])
        for tree in trees:
            v = tree["val_nova_syst"].array()
            v = np.array(v)
            v = np.transpose(v)
            ret.append([index][params.burnin:])
            
    if params.do_compare:
        ret = np.array(ret)
        for i in ret:
            print (len(i)/1000, "k entries")
        return ret
    else:
        ret = np.array(ret).flatten()
        print (len(ret)/1000, "k entries")
        return ret
        

    
def hdi(bin_val, prob):
    initial_shape = bin_val.shape
    bin_val = bin_val.flatten()
    
    s = np.sum(bin_val)
    if s>1.01 or s<0.99:
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
        
        
def get_2dhisto(varx, vary, apply_change_of_var=True):
    stepx = (varx.range[1] - varx.range[0])/ varx.nbins
    stepy = (vary.range[1] - vary.range[0])/ vary.nbins
    
    bx = np.arange(varx.range[0], varx.range[1]+stepx, stepx)
    by = np.arange(vary.range[0], vary.range[1]+stepy, stepy)
    bx = bx[:varx.nbins+1]
    by = by[:vary.nbins+1]
    
    valuesx = np.array(get_var_array(varx))
    valuesy = np.array(get_var_array(vary))
    
    if apply_change_of_var:
        valuesx = varx.change_of_var(valuesx)
        valuesy = vary.change_of_var(valuesy)

    weights = [1./len(valuesx)] * len(valuesx)
    
    histo, binsx, binsy = np.histogram2d(valuesx, valuesy, bins=[bx,by], weights=weights)
    return histo, binsx, binsy

def get_histo(var, apply_change_of_var=True):
    step = (var.range[1] - var.range[0])/ var.nbins
    b = np.arange(var.range[0], var.range[1]+step, step)

    values = np.array(get_var_array(var))
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



# Print iterations progress
def printProgressBar (iteration, total, prefix='', suffix='', decimals=1, length=100, fill='█', printEnd="\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()

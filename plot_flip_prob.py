import mcmc_utils as mcmcu
import uproot as up
import matplotlib
import numpy as np
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

tree = up.open("main.root")["result"]

with PdfPages('mh_variation.pdf') as pdf:
    var = mcmcu.osc_variables[1]
    values = var.change_of_var(mcmcu.get_var_array(tree, var))
    plt.plot(values, label="MH time serie")
    values = (values+1)/2
    values = np.abs(np.diff(values))
    N=10000
    factor=50.
    prob = np.convolve(values, np.ones(N)/N, mode='valid')
    plt.plot(prob*factor, label=r"Average flipping probability $\times "+str(factor)+"$") # I think this is correct...
    plt.xlabel('Step')
    plt.ylabel(var.nice_name+" "+var.unit)
    plt.legend()
    plt.title(var.nice_name+' variations NOvA-only (all steps after burnin)')
    pdf.savefig()
    plt.close()
    

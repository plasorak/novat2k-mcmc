import uproot as up
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import mcmc_utils as mcmcu


tree = up.open(mcmcu.input_file)["result"]
likelihood = tree["likelihood"     ].array()
acceptance = tree["acceptance_rate"].array()
step       = tree["step"           ].array()

def window_rms(a, window_size):
    a2 = np.power(a,2)
    window = np.ones(window_size)/float(window_size)
    return np.sqrt(np.convolve(a2, window, 'valid'))

with PdfPages('burnin.pdf') as pdf:
    
    for var in mcmcu.osc_variables:
        plt.plot(step, tree[var.tree_var_name].array())
        plt.xlabel('Step')
        plt.ylabel(var.nice_name+" "+var.unit)
        pdf.savefig()
        plt.close()

    
    likelihood = np.array(likelihood)
    plt.plot(step, likelihood)
    average_llh = np.convolve(likelihood, np.ones(100)/100, mode='valid')
    plt.plot(step[:len(average_llh)], average_llh)
    plt.xlabel('Step')
    plt.ylabel('Likelihood')
    pdf.savefig()
    plt.close()

    plt.plot(step, acceptance)
    plt.xlabel('Step')
    plt.ylabel('acceptance rate')
    pdf.savefig()
    plt.close()

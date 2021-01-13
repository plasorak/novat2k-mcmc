import mcmc_utils as mcmcu
import numpy as np
import uproot as up
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

tree = up.open("result.root")["result"]
color="black"
color68="black"
color90="dimgrey"
color95="lightgrey"

with PdfPages('nova_syst_pars.pdf') as pdf:
    for var in mcmcu.nova_syst_variables:
        print(var.tree_var_name)
        (hist, bins) = mcmcu.get_histo(tree, var)
        ci68 = mcmcu.hdi(hist, prob=0.682689)
        ci90 = mcmcu.hdi(hist, prob=0.90    )
        ci95 = mcmcu.hdi(hist, prob=0.9545  )
        hist68 = np.zeros(len(hist))
        hist90 = np.zeros(len(hist))
        hist95 = np.zeros(len(hist))
            
        for i in range(len(ci68)):
            hist68[i] = ci68[i]*hist[i]
            hist90[i] = ci90[i]*hist[i]
            hist95[i] = ci95[i]*hist[i]
                
        plt.hist(bins[:-1], bins, weights=hist  , histtype='step', color="black")
        plt.hist(bins[:-1], bins, weights=hist95, color=color95, label="CI=95%")
        plt.hist(bins[:-1], bins, weights=hist90, color=color90, label="CI=90%")
        plt.hist(bins[:-1], bins, weights=hist68, color=color68, label="CI=68%")
        plt.legend()
        
        plt.xlabel(var.nice_name+" "+var.unit)
        plt.ylabel('Posterior Probability')
        plt.title(var.nice_name+' NOvA-only')
        pdf.savefig()
        plt.close()


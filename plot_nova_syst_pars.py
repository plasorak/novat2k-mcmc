import mcmc_utils as mcmcu
import numpy as np
import uproot as up
import matplotlib
import math
import scipy.stats as stats
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

tree = up.open(mcmcu.output)["result"]
color="black"
colorPrior="lightslategrey"
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
        plt.hist(bins[:-1], bins, weights=hist95, color=color95, label="Cred Interv=95\%")
        plt.hist(bins[:-1], bins, weights=hist90, color=color90, label="Cred Interv=90\%")
        plt.hist(bins[:-1], bins, weights=hist68, color=color68, label="Cred Interv=68\%")
        mu = 0
        variance = 1
        sigma = math.sqrt(variance)
        x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
        y = stats.norm.pdf(x, mu, sigma)
        y_max = y.max()
        y = y/y_max*hist.max()
        plt.plot(x, y, label="Prior: Gauss(0,1)", color=colorPrior)
        plt.legend()
        
        plt.xlabel(var.nice_name)
        plt.ylabel('Posterior probability')
        plt.title(var.nice_name)
        pdf.savefig()
        plt.close()

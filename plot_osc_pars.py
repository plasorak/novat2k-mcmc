import mcmc_utils as mcmcu
import uproot as up
import matplotlib
import numpy as np
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

tree = up.open("result.root")["result"]
color="black"
color68="black"
color90="dimgrey"
color95="lightgrey"

with PdfPages('osc_pars.pdf') as pdf:
    for var in mcmcu.osc_variables:

        print (var.tree_var_name)
        dcp = var.tree_var_name.find("cp")>0
        (hist, bins) = mcmcu.get_histo(tree, var)
        
        if var.nice_name.find("hierarchy")>0:
            _, ax = plt.subplots()
            ax.bar([0,1], hist, 1)
            bin_label = ['Inverted','Normal']
            ax.set_xticks([0,1])
            ax.set_xticklabels(bin_label,rotation=45, rotation_mode="anchor", ha="right")
            print("Inverted ordering probability",hist[0])
            print("Normal ordering probability",hist[1])
            
        else:
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
            
        print(len(hist),len(bins),np.sum(hist))
        plt.xlabel(var.nice_name+" "+var.unit)
        plt.ylabel('Posterior Probability')
        plt.title(var.nice_name+' NOvA-only')
        pdf.savefig()
        plt.close()

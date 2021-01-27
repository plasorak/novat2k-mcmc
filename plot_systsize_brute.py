
import mcmc_utils as mcmcu
import uproot as up
import hist
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


tree_no_syst = up.open("no_syst.root")["result"]
tree_calib   = up.open("only_syst42.root")["result"]
#color95="grey"

with PdfPages('systsize_brute.pdf') as pdf:

    for var_osc in mcmcu.osc_variables:
        if (var_osc.tree_var_name != "val_dm32"):
            continue
        var_osc.nbins = 100
        
        osc_0_syst, bins = mcmcu.get_histo(tree_no_syst, var_osc)
        osc_1_syst, bins = mcmcu.get_histo(tree_calib  , var_osc)
        
        ci68_0_syst = mcmcu.hdi(osc_0_syst, prob=0.682689)
        ci68_1_syst = mcmcu.hdi(osc_1_syst, prob=0.682689)

        hist_0_syst = np.zeros(len(osc_0_syst))
        hist_1_syst = np.zeros(len(osc_1_syst))

        min_0_syst  = 100000
        max_0_syst  = 0
        best_0_syst = 0
        min_1_syst  = 100000
        max_1_syst  = 0
        best_1_syst = 0
        
        for i in range(len(ci68_0_syst)):
            if ci68_0_syst[i]>0:
                if bins[i]>max_0_syst: max_0_syst = bins[i]
                if bins[i]<min_0_syst: min_0_syst = bins[i]
                
            if ci68_1_syst[i]>0:
                if bins[i]>max_1_syst: max_1_syst = bins[i]
                if bins[i]<min_1_syst: min_1_syst = bins[i]
            
            if osc_0_syst[i]>osc_0_syst[best_0_syst]: best_0_syst = i
            if osc_1_syst[i]>osc_1_syst[best_1_syst]: best_1_syst = i
                
            hist_0_syst[i] = ci68_0_syst[i]*osc_0_syst[i]
            hist_1_syst[i] = ci68_1_syst[i]*osc_1_syst[i]

        best_0_syst = bins[best_0_syst]
        best_1_syst = bins[best_1_syst]

        error_pos_0_syst = max_0_syst-best_0_syst
        error_neg_0_syst = best_0_syst-min_0_syst
        error_pos_1_syst = max_1_syst-best_1_syst
        error_neg_1_syst = best_1_syst-min_1_syst

        error_on_error_pos_0_syst = 2.*(bins[1]-bins[0])
        error_on_error_neg_0_syst = 2.*(bins[1]-bins[0])
        error_on_error_pos_1_syst = 2.*(bins[1]-bins[0])
        error_on_error_neg_1_syst = 2.*(bins[1]-bins[0])

        print("Error:")
        print(" - 0 syst neg: ", error_neg_0_syst, "+/-", error_on_error_pos_0_syst)
        print(" - 0 syst pos: ", error_pos_0_syst, "+/-", error_on_error_neg_0_syst)
        print(" - 1 syst neg: ", error_neg_1_syst, "+/-", error_on_error_pos_1_syst)
        print(" - 1 syst pos: ", error_pos_1_syst, "+/-", error_on_error_neg_1_syst)
        
        plt.hist(bins[:-1], bins, weights=osc_0_syst, histtype='step', label="0 syst")
        plt.hist(bins[:-1], bins, weights=osc_1_syst, histtype='step', label="1 syst")

        
        plt.xlabel(var_osc.nice_name+" "+var_osc.unit)
        plt.ylabel('Posterior probability')
        plt.legend()
        
        pdf.savefig()
        plt.close()

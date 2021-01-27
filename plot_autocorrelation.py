import mcmc_utils as mcmcu
import uproot as up
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

tree = up.open("result.root")["result"]


with PdfPages('pars_autocorr.pdf') as pdf:
    def for_each_var(var):
        print(var.tree_var_name)
        x = mcmcu.get_var_array(tree, var)
        acf_data=[mcmcu.acf(x, lag) for lag in range(5000)]
        plt.plot(acf_data)
        plt.title(var.nice_name+' ACF')
        plt.xlabel('Step lag')
        plt.ylabel(var.nice_name+' auto-correlation')
        pdf.savefig()
        plt.close()

    for var in mcmcu.osc_variables:
        for_each_var(var)

    for var in mcmcu.nova_syst_variables:
        for_each_var(var)

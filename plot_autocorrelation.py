import mcmc_utils as mcmcu
import uproot as up
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

tree = up.open("result.root")["result"]
    
def acf(x, lag):
    # Slice the relevant subseries based on the lag
    if len(x)-lag<=0:
        pass
    y1 = x[:(len(x)-lag)]
    y2 = x[lag:]
    # print(lag,len(y1), len(y2))
    # Subtract the mean of the whole series x to calculate Cov
    mean=np.mean(x)
    sum_product = np.sum((y1-mean)*(y2-mean))
    # Normalize with var of whole series
    return sum_product / ((len(x) - lag) * np.var(x))

with PdfPages('pars_autocorr.pdf') as pdf:
    def for_each_var(var):
        x = mcmcu.get_var_array(tree, var)
        acf_data=[acf(x, lag) for lag in range(5000)]
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

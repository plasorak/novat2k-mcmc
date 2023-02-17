import mcmc_utils as mcmcu
import uproot as up
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

arg=mcmcu.arguments("autocorrelation")
arg.parser.add_argument("--max-lag", nargs=1, type=int, help="maximum lag for autocorrelation plot", default=1000)
arg.parse()
maxlag=arg.argument.max_lag[0]

output=arg.output

with PdfPages(output) as pdf:

    def get_acf(var):
        print(var.tree_var_name)
        data = mcmcu.get_var_array(var)

        if arg.do_compare:
            count=0

            for x in data:
                acf_data=[]
                mcmcu.printProgressBar(0, maxlag, prefix='Progress:', suffix='Complete', length=50)
                for lag in range(maxlag):
                    acf_data.append(mcmcu.acf(x, lag))
                    mcmcu.printProgressBar(lag, maxlag, prefix='Progress:', suffix='Complete', length=50)
                print("\n")
                plt.plot(acf_data, label=arg.label[count])
                count+=1
                plt.legend()
        else:
            acf_data=[]

            mcmcu.printProgressBar(0, maxlag, prefix='Progress:', suffix='Complete', length=50)
            for lag in range(maxlag):
                acf_data.append(mcmcu.acf(data, lag))
                mcmcu.printProgressBar(lag, maxlag, prefix='Progress:', suffix='Complete', length=50)
            print("\n")
            plt.plot(acf_data)
            
        plt.title(var.nice_name+' ACF')
        plt.xlabel('Step lag')
        plt.ylabel(var.nice_name+' auto-correlation')
        pdf.savefig()
        plt.close()

    interesting_vars = [mcmcu.llh, mcmcu.dcp, mcmcu.s2th23, mcmcu.s2th13, mcmcu.dm32]
    
    for var in interesting_vars:
        get_acf(var)

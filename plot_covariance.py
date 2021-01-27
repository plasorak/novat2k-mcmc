import mcmc_utils as mcmcu
import uproot as up
import hist
import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


tree = up.open("main.root")["result"]

with PdfPages('correlation.pdf') as pdf:
    
    for var_osc in mcmcu.osc_variables:
        if (var_osc.tree_var_name == "val_mh" or
            var_osc.nice_name == r"$\theta_{13}$" or
            var_osc.nice_name == r"$\theta_{23}$"):
            continue
        
        osc = mcmcu.get_var_array(tree, var_osc)
        print(var_osc.nice_name)
        for type_syst, systs in mcmcu.grouped_nova_syst_variables.items():
            covs = {}
            
            print("   "+type_syst)

            for i, var_syst in enumerate(systs):
                syst = mcmcu.get_var_array(tree, var_syst)
                cov = np.cov(osc, syst)
                covs[var_syst.nice_name] = cov[0,1]/np.sqrt(cov[0,0]*cov[1,1])
                #print(var_osc.tree_var_name,var_syst.nice_name,cov, covs[var_syst.nice_name])
                                
            _, ax = plt.subplots()
            keys   = list(reversed(covs.keys()))
            values = list(reversed(covs.values()))
            plt.barh(range(len(covs)), values, 0.7)
            ax.xaxis.grid(True, linestyle='--', which='major',
                           color='grey', alpha=.25)
            ax.yaxis.grid(True, linestyle='--', which='major',
                           color='grey', alpha=.25)
            ax.axvline(0, color='black')
            ax.set_yticks(range(len(covs)))
            ax.set_yticklabels(keys)
            plt.xlabel('Correlation')
            plt.title('Correlations of '+var_osc.nice_name+" to "+type_syst+" systematics" )
            plt.tight_layout()
            
            pdf.savefig()
            plt.close()
        

    

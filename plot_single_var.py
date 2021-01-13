import mcmc_utils as mcmcu
import uproot as up
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# tree = up.open("result.root")["result"]
tree = up.open("result.root")["result"]

with PdfPages('osc_pars_var.pdf') as pdf:
    for var in mcmcu.osc_variables:
        print (var.tree_var_name)
        values = var.change_of_var(mcmcu.get_var_array(tree, var))[:1000]
        plt.plot(values)
        plt.xlabel('Step')
        plt.ylabel(var.nice_name+" "+var.unit)
        plt.title(var.nice_name+' variations NOvA-only')
        pdf.savefig()
        plt.close()


    for var in mcmcu.nova_syst_variables:
        print (var.tree_var_name)
        values = var.change_of_var(mcmcu.get_var_array(tree, var))[:1000]
        plt.plot(values)
        plt.xlabel('Step')
        plt.ylabel(var.nice_name+" "+var.unit)
        plt.title(var.nice_name+' variations NOvA-only')
        pdf.savefig()
        plt.close()


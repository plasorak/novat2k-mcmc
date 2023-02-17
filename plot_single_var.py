import mcmc_utils as mcmcu
import uproot as up
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# tree = up.open("result.root")["result"]
arg=mcmcu.arguments("traces")
arg.parse()
output=arg.output

with PdfPages(output) as pdf:
    interesting_vars=[]
    if arg.argument.cafana:
        interesting_vars = [mcmcu.cafana_llh ,
                            mcmcu.cafana_dcp ,
                            mcmcu.cafana_mh  ,
                            mcmcu.cafana_th13,
                            mcmcu.cafana_th23,
                            mcmcu.cafana_dm32]
    else:
        interesting_vars = [mcmcu.llh, mcmcu.dcp, mcmcu.s2th23, mcmcu.s2th13, mcmcu.dm32]
    # interesting_vars = [mcmcu.llh]
    for var in interesting_vars:
        print (var.tree_var_name)
        print(mcmcu.params.nsteps_max)
        values = var.change_of_var(mcmcu.get_var_array(var))
        plt.plot(values)
        plt.xlabel('Step')
        plt.ylabel(var.nice_name+" "+var.unit)
        plt.title(var.nice_name+f' variations NOvA-only burn-in: {mcmcu.params.burnin} n steps: {mcmcu.params.nsteps_max}')
        pdf.savefig()
        plt.close()
        
    # for var in mcmcu.nova_syst_variables:
    #     print (var.tree_var_name)
    #     values = var.change_of_var(mcmcu.get_var_array(var))[mcmcu.params.burnin:mcmcu.params.nsteps_max]
    #     plt.plot(values)
    #     plt.xlabel('Step')
    #     plt.ylabel(var.nice_name+" "+var.unit)
    #     plt.title(var.nice_name+f' variations NOvA-only burn-in: {mcmcu.params.burnin} n step max: {mcmcu.params.nsteps_max}')
    #     pdf.savefig()
    #     plt.close()

import mcmc_utils as mcmcu
import uproot as up
import hist
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


ss = [0.01,
      0.02,
      0.05,
      0.07,
      0.1,
      0.15,
      0.20,
      0.30]
sss = ["0.01",
       "0.02",
       "0.05",
       "0.07",
       "0.1",
       "0.15",
       "0.20",
       "0.30"]

trees = []
legname = []
acc= []

for i, s in enumerate(sss):
    fn = "result_500k_ss_"+str(s)+".root"
    file = up.open(fn)
    tr = file["result"]
    ac = file["acceptance_rate"].to_boost()
    trees.append(tr)
    print(ac[0])# unlike root histogram, there isn't any overflow and underflow, so [BOOST-HIST bin 0] -> [ROOT bin 1]
    acc.append(ac[0])
    legname.append("Step-size = "+str(ss[i]))

n_steps = 50000

with PdfPages('stepsize.pdf') as pdf:
    for i, tree in enumerate(trees):
        likelihood = tree["likelihood"].array()[:n_steps]
        plt.plot(likelihood, label=legname[i])
        plt.xlabel('Step')
        plt.ylabel('Likelihood')
    plt.legend()
    pdf.savefig()
    plt.close()
    
    plt.errorbar(ss, [d.value for d in acc], yerr=[d.variance for d in acc])
    plt.xlabel('NOvA systematics step size')
    plt.ylabel('Acceptance rate')
    pdf.savefig()
    plt.close()

    for i, var in enumerate(mcmcu.nova_syst_variables):
        print(i)
        acf_data_1k = []
        for i, tree in enumerate(trees):
            x = mcmcu.get_var_array(tree, var)
            acf_data=[mcmcu.acf(x, lag) for lag in range(10000)]
            acf_data_1k.append(mcmcu.acf(x, 1000))
            plt.plot(acf_data, label=legname[i])
            plt.title(var.nice_name+' ACF')
            plt.xlabel('Step lag')
            plt.ylabel(var.nice_name+' auto-correlation')
        plt.legend()
        pdf.savefig()
        plt.close()


        plt.errorbar(ss, acf_data_1k)
        plt.xlabel('NOvA systematics step size')
        plt.ylabel('ACF on '+var.nice_name+' after 1000 steps')
        pdf.savefig()
        plt.close()

    
        

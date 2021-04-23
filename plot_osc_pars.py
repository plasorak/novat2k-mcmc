import mcmc_utils as mcmcu
import uproot as up
import matplotlib
import numpy as np
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import ROOT as root

arg=mcmcu.arguments("osc_pars")
#print(arg.calculate_yourself(60))

color="black"
colorPrior="lightslategrey"
color68="black"
color90="dimgrey"
color95="lightgrey"

file_out = root.TFile(arg.output.replace("pdf","root"), "RECREATE")

output=arg.output

with PdfPages(output) as pdf:
    interesting_vars = [mcmcu.dcp, mcmcu.s2th23, mcmcu.s2th13, mcmcu.dm32, mcmcu.mh]
    # interesting_vars = [mcmcu.mh]
    
    for var in interesting_vars:
        
        print (var.tree_var_name)
        dcp = var.tree_var_name.find("cp")>0
        (hist, bins) = mcmcu.get_histo(var, arg.input_trees, arg.burnin)
        
        histo = root.TH1D(var.tree_var_name, var.nice_name, len(bins), bins.min(), bins.max())
        for i,content in enumerate(hist):
            histo.SetBinContent(i+1, content)
        histo.Write()
        
        if var.nice_name.find("hierarchy")>0:
            _, ax = plt.subplots()
            ax.bar([0,1], hist, 1, color=color90, label="MH posterior")
            bin_label = ['Inverted','Normal']
            ax.set_xticks([0,1])
            ax.set_xticklabels(bin_label,rotation=45, rotation_mode="anchor", ha="right")
            print("Inverted ordering probability",hist[0])
            print("Normal ordering probability",hist[1])
            ax.plot([-0.5, 1.5], [0.5, 0.5], label="Prior: 50\%-50\% Normal/Inverted", color=colorPrior)
            plt.legend(loc="upper left")
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
            plt.hist(bins[:-1], bins, weights=hist95, color=color95, label="Cred Interv=95\%")
            plt.hist(bins[:-1], bins, weights=hist90, color=color90, label="Cred Interv=90\%")
            plt.hist(bins[:-1], bins, weights=hist68, color=color68, label="Cred Interv=68\%")
            if var.nice_name.find("\sin^{2}")>0:
                prior = 1. / (np.sqrt(bins * (1-bins)))
                skimmed_bins = []
                skimmed_prior = []
                for i in range(len(prior)):
                    if prior[i] <10000 and prior[i] > -10000:
                        skimmed_bins .append(bins[i])
                        skimmed_prior.append(prior[i])
                prior = np.array(skimmed_prior)
                binss = np.array(skimmed_bins )
                middle_bin = int((len(prior)-1)/2)
                # prior = prior / pir(0.5*(prior[middle_bin]+prior[middle_bin+1])) * hist[middle_bin] * 0.1
                prior = prior / prior.max() * hist.max() * 0.1
                
                plt.plot(binss, prior, label="Prior: Flat in angle", color=colorPrior)
            else:
                plt.plot([bins[0], bins[-1]], [hist.max()/10, hist.max()/10], label="Prior: Flat", color=colorPrior)
            plt.legend()
            
        plt.xlabel(var.nice_name+" "+var.unit)
        plt.ylabel('Posterior Probability')
        plt.title(var.nice_name+' NOvA-only')
        pdf.savefig()
        plt.close()
        chi2 = []
        for d in hist:
            if d>0:
                chi2.append(-2*np.log(d))
            else:
                chi2.append(np.inf)
                
        chi2 -= np.min(chi2)
        maxchi2=np.max(chi2[chi2!=np.inf])
        for i, d in enumerate(chi2):
            if d==np.inf:
                chi2[i]=maxchi2
            else:
                chi2[i]=d
        
        plt.hist(bins[:-1], bins, weights=chi2, histtype='step', color="black")
        plt.xlabel(var.nice_name+" "+var.unit)
        plt.ylabel('\"$\chi^2$\"')
        plt.title(var.nice_name+' NOvA-only')
        pdf.savefig()
        plt.close()

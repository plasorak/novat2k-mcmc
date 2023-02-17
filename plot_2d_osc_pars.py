import mcmc_utils as mcmcu
import uproot as up
import matplotlib
import numpy as np
import scipy.stats as stats
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

color68="black"
color90="dimgrey"
color95="lightgrey"

import ROOT as root
file_out = root.TFile("for_artur_2d.root", "RECREATE")

with PdfPages('osc_pars_2d.pdf') as pdf:
    interesting_pairs = [[mcmcu.s2th23, mcmcu.dm32  ],
                         [mcmcu.dcp   , mcmcu.s2th23],
                         [mcmcu.s2th13, mcmcu.dcp   ]]
    
    for pair in interesting_pairs:
        hist, binsx, binsy = mcmcu.get_2dhisto(pair[0], pair[1])

        ci68 = mcmcu.hdi(hist, prob=0.682689)
        ci90 = mcmcu.hdi(hist, prob=0.90    )
        ci95 = mcmcu.hdi(hist, prob=0.9545  )
        
        n_binx = len(binsx)-1
        n_biny = len(binsy)-1
        
        mock_data_x = [binsx[:-1]]*n_biny
        mock_data_y = [binsy[:-1]]*n_binx
        mock_data_x = np.array(mock_data_x).transpose().flatten()
        mock_data_y = np.array(mock_data_y).flatten()
        
        cnts68,_,_,_ = plt.hist2d(x=mock_data_x, y=mock_data_y, weights=ci68.flatten(), bins=[binsx, binsy])
        cnts90,_,_,_ = plt.hist2d(x=mock_data_x, y=mock_data_y, weights=ci90.flatten(), bins=[binsx, binsy])
        cnts95,_,_,_ = plt.hist2d(x=mock_data_x, y=mock_data_y, weights=ci95.flatten(), bins=[binsx, binsy])
        
        plt.close()
        
        cnts68 = cnts68.transpose()
        cnts90 = cnts90.transpose()
        cnts95 = cnts95.transpose()

        cs68 = plt.contour(cnts68,levels=[0.5],extent=[binsx.min(),binsx.max(),binsy.min(),binsy.max()], linewidths=1, colors=color68)
        cs90 = plt.contour(cnts90,levels=[0.5],extent=[binsx.min(),binsx.max(),binsy.min(),binsy.max()], linewidths=1, colors=color90)
        cs95 = plt.contour(cnts95,levels=[0.5],extent=[binsx.min(),binsx.max(),binsy.min(),binsy.max()], linewidths=1, colors=color95)
        h68,_ = cs68.legend_elements()
        h90,_ = cs90.legend_elements()
        h95,_ = cs95.legend_elements()
        plt.legend([h68[0],h90[0],h95[0]], ['CI=68%','CI=90%','CI=95%'])

        histo68 = root.TH2D("68CL_"+pair[0].tree_var_name+"_"+pair[1].tree_var_name,
                            ";"+pair[0].nice_name+";"+pair[1].nice_name,
                            len(binsx)-1, binsx.min(), binsx.max(),
                            len(binsy)-1, binsy.min(), binsy.max())
        histo90 = root.TH2D("90CL_"+pair[0].tree_var_name+"_"+pair[1].tree_var_name,
                            ";"+pair[0].nice_name+";"+pair[1].nice_name,
                            len(binsx)-1, binsx.min(), binsx.max(),
                            len(binsy)-1, binsy.min(), binsy.max())
        histo95 = root.TH2D("95CL_"+pair[0].tree_var_name+"_"+pair[1].tree_var_name,
                            ";"+pair[0].nice_name+";"+pair[1].nice_name,
                            len(binsx)-1, binsx.min(), binsx.max(),
                            len(binsy)-1, binsy.min(), binsy.max())
        
        for x in range(ci68.shape[0]):
            for y in range(ci68.shape[1]):
                bx = histo68.GetXaxis().FindBin(binsx[x]+0.0000001)
                by = histo68.GetYaxis().FindBin(binsy[y]+0.0000001)
                histo68.SetBinContent(bx,by, ci68[x, y])
                histo90.SetBinContent(bx,by, ci90[x, y])
                histo95.SetBinContent(bx,by, ci95[x, y])
                
        histo68.SetOption("COLZ")
        histo90.SetOption("COLZ")
        histo95.SetOption("COLZ")
        
        histo68.Write()
        histo90.Write()
        histo95.Write()

        plt.xlabel(pair[0].nice_name+" "+pair[0].unit)
        plt.ylabel(pair[1].nice_name+" "+pair[1].unit) 
        plt.title(pair[0].nice_name+"-"+pair[1].nice_name+' NOvA-only')
        pdf.savefig()
        plt.close()
        

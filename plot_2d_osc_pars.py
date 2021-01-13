import mcmc_utils as mcmcu
import uproot as up
import matplotlib
import numpy as np
import scipy.stats as stats
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

tree = up.open("result.root")["result"]
color68="black"
color90="dimgrey"
color95="lightgrey"


with PdfPages('osc_pars_2d.pdf') as pdf:
    interesting_pairs = [[mcmcu.osc_variables[4],mcmcu.osc_variables[6]],
                         [mcmcu.osc_variables[0],mcmcu.osc_variables[4]]]
    
    for pair in interesting_pairs:
        hist, binsx, binsy = mcmcu.get_2dhisto(tree, pair[0], pair[1])
        ci68 = mcmcu.hdi(hist, prob=0.682689).flatten()
        ci90 = mcmcu.hdi(hist, prob=0.90    ).flatten()
        ci95 = mcmcu.hdi(hist, prob=0.9545  ).flatten()
        
        n_binx = len(binsx)-1
        n_biny = len(binsy)-1
        
        mock_data_x = [binsx[:-1]]*n_biny
        mock_data_y = [binsy[:-1]]*n_binx
        mock_data_x = np.array(np.array(mock_data_x).flatten())
        mock_data_y = np.array(np.array(mock_data_y).transpose().flatten())
        
        cnts68,_,_,_ = plt.hist2d(x=mock_data_x, y=mock_data_y, weights=ci68, bins=[binsx, binsy])
        cnts90,_,_,_ = plt.hist2d(x=mock_data_x, y=mock_data_y, weights=ci90, bins=[binsx, binsy])
        cnts95,_,_,_ = plt.hist2d(x=mock_data_x, y=mock_data_y, weights=ci95, bins=[binsx, binsy])
        plt.close()
        
        # God only knows why I need to transpose cnts...
        cs68 = plt.contour(cnts68.transpose(),levels=[0.5],extent=[binsx.min(),binsx.max(),binsy.min(),binsy.max()], linewidths=1, colors=color68)
        cs90 = plt.contour(cnts90.transpose(),levels=[0.5],extent=[binsx.min(),binsx.max(),binsy.min(),binsy.max()], linewidths=1, colors=color90)
        cs95 = plt.contour(cnts95.transpose(),levels=[0.5],extent=[binsx.min(),binsx.max(),binsy.min(),binsy.max()], linewidths=1, colors=color95)
        h68,_ = cs68.legend_elements()
        h90,_ = cs90.legend_elements()
        h95,_ = cs95.legend_elements()
        plt.legend([h68[0],h90[0],h95[0]], ['CI=68%','CI=90%','CI=95%'])
        
        plt.xlabel(pair[0].nice_name+" "+pair[0].unit)
        plt.ylabel(pair[1].nice_name+" "+pair[1].unit) 
        plt.title(pair[0].nice_name+"-"+pair[1].nice_name+' NOvA-only')
        pdf.savefig()
        plt.close()


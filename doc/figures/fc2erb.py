import os
import numpy as np
from matplotlib import pyplot

scriptDirectory = os.path.dirname(os.path.realpath(__file__)) # must run entire file

numfc = 100 # number of center frequencies
fc = np.logspace(2.0, 4.0, num=numfc) / 1000 # log-spaced center frequencies in kHz

w1 = 24.7 * (4.37 * fc + 1)
w2 = 6.23 * np.power(fc, 2) + 93.39 * fc + 28.52

#%% Plot curves and save figure
fc2erbPlot = pyplot.figure()
pyplot.rcParams["font.family"] = "Times New Roman"
pyplot.rcParams["axes.formatter.useoffset"] = False
pyplot.plot(fc, w1, 'k')
pyplot.plot(fc, w2, 'k--')
pyplot.xscale('log')
pyplot.yscale('log')
pyplot.xlim(fc[1]/1.5,fc[-1]*1.5)
pyplot.ylim(10,5000)
pyplot.xticks((0.1,0.2,0.5,1,2,5,10), ('0.1', '0.2', '0.5', '1', '2', '5', '10'), fontsize=12)
pyplot.yticks((10,50,100,500,1000,5000), ('10', '50', '100', '500', '1000', '5000'), fontsize=12)
pyplot.xlabel('Center Frequency (kHz)', fontsize=16)
pyplot.ylabel('ERB (Hz)', fontsize=16)
pyplot.text(6, 400, 'Linear', fontsize=14)
pyplot.text(3, 1300, 'Quadratic', fontsize=14)
pyplot.show()
fc2erbPlot.savefig(os.path.join(scriptDirectory,'fc2erb.eps'), format='eps', dpi=1200)
import uproot as up
import numpy as np
import matplotlib.pyplot as plt

tree = up.open("result.root")["result"]

likelihood = tree["likelihood"].array()[:1000]
step       = tree["step"      ].array()[:1000]

plt.plot(step, likelihood)
plt.xlabel('Step')
plt.ylabel('Likelihood')

plt.savefig("likelihood.pdf")


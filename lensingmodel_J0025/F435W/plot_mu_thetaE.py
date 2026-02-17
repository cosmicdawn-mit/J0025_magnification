import os, sys
import numpy as np
import matplotlib.pyplot as plt

sigma, mu, thetaE = np.loadtxt('muthetaE.txt')

plt.plot(thetaE, mu)
plt.show()




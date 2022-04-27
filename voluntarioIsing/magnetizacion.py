# -*- coding: utf-8 -*-



# Importamos el modulo matplotlib y numpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

data=np.loadtxt('magnetizacion.txt')
y0 = data[:,0]
y1 = data[:,1]
y2 = data[:,2]
y3 = data[:,3]
x = data[:,4]  #T

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x,y0,color='red',label='N=16',linewidth=.25)
ax.plot(x,y1,color='blue',label='N=32',linewidth=.25)
ax.plot(x,y2,color='purple',label='N=64',linewidth=.25)
ax.plot(x,y3,color='yellow',label='N=128',linewidth=.25)
plt.legend()
#ax.title('Magnetizacion en funcion de N')
#ax.xlabel('T')
#ax.ylabel('m')
plt.show()
#fig.tight_layout()
fig.savefig('magnetizacion.pdf')


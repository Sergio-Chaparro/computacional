# Importamos el módulo matplotlib y numpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# plt.style.use('./style.mplstyle.patch')
# plt.style.use('dark_background')

# Leemos datos de un fichero
data = np.loadtxt('Variables.txt') # 6000 puntos del atractor de Lorenz (6000x3)
x1 = data[:,0]
y1 = data[:,1]
x2 = data[:,3]
y2 = data[:,4]
x3 = data[:,5]
y3 = data[:,6]
x4 = data[:,7]
y4 = data[:,8]
x5 = data[:,9]
y5 = data[:,10]

# Creamos una figura y hacemos nuestro primer plot
plt.figure()
plt.plot(x,y)
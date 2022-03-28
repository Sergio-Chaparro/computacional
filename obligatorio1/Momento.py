# Importamos el módulo matplotlib y numpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# plt.style.use('./style.mplstyle.patch')
# plt.style.use('dark_background')

# Leemos datos de un fichero
data = np.loadtxt('momento.txt') # 6000 puntos del atractor de Lorenz (6000x3)
x = data[:,0]
y = data[:,1]

# Creamos una figura y hacemos nuestro primer plot
plt.figure()
plt.plot(x,y)
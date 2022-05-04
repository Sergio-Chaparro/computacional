# Importamos el módulo matplotlib y numpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# Leemos datos de un fichero
data = np.loadtxt('norma.txt') 
x = data[:,0]
y = data[:,1]

# Creamos una figura y hacemos nuestro primer plot
plt.figure()
plt.plot(x,y)

plt.savefig('norma.png', dpi=300)
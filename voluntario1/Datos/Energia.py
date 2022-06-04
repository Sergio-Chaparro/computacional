# Importamos el módulo matplotlib y numpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# Leemos datos de un fichero
data = np.loadtxt('Energia.txt') 
x = data[:,0]
y = data[:,1]
z = data[:,2]
t = data[:,3]

# Creamos una figura y hacemos nuestro primer plot
plt.figure()
plt.plot(x,y,label = "Energia usada en el lanzamiento")
plt.plot(x,z,label = "Energia usada durante el viaje")
plt.plot(x,t,label = "Energia total")
plt.xlabel("Lugar de la explosion en distancia entre la Tierra y la Luna")
plt.ylabel("Energia en Kilotones")
plt.legend()

plt.savefig('Energia.png', dpi=300)
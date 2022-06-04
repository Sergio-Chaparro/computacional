# Importamos el módulo matplotlib y numpy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# Leemos datos de un fichero
data = np.loadtxt('Bomba.txt') 
x = data[:,0]
y = data[:,1]

# Creamos una figura y hacemos nuestro primer plot
plt.figure()
plt.plot(x,y,"r")
plt.xlabel("Lugar de la explosion en distancia entre la Tierra y la Luna")
plt.ylabel("Energia en megatones")
plt.title("Energia de las explosiones")

plt.savefig('Bomba.png', dpi=300)
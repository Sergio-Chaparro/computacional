//Programa para la simulacion del sistema solar, con 8 planetas


#include <math.h>
#include <stdio.h>

//defino las constantes que me harán falta para reescalar las variables 
#define G 6.67E-11
#define c 1.496E11
#define Ms 1.99E30


void DatosIniciales(double r[2][],double v[2][], double a[2][],double m[],double t,int N)
void reescalar(double r[2][],double m[],double t,int N)
void desescalar(double r[2][],double m[],double t,int N)
void Algoritmo(double r[2][],double v[2][], double a[2][],double m[],int N)
void Simulacion(double r[2][],double v[2][], double a[2][],double m[],double t,int N, double h, double tmax)



int main(void)
{
    //Declaro el tiempo y el paso que vamos a utilizar y las inicializo
    double h,tmax;
    int N;
    N=8;
    h=0.01;
    tmax=5;

    //Declaramos la posicion, la velocidad y la aceleracion
    //que vamos a utilizar para los N planetas
    double a[2][N],v[2][N],r[2][N], m[N],t;

    //Inicializo todos los vectores
    DatosIniciales(r,v,a,m,t,N);

    //Reescalamos los datod para poder tratarlos con mas facilidad
    //Usamos masas solares y la distancia tierra-sol
    reescalar(r,m,t,N);//A partir de aqui estaran en esas unidades

    //Podemos empezar a simular, con parametros h y tmax
    
    Simulacion(r,v,a,m,t,N,h,tmax);

    // Podemos desescalar los datos tras la simulacion
    desescalar(r,m,t,N);
    
    return 0;
}


//Funcion que aporta los valores iniciales a las variables
void DatosIniciales(double r[2][],double v[2][], double a[2][],double m[],double t)
{

}


//Funcion que reescala las variables para ser usadas con valores mas sencillos
void reescalar(double r[2][],double m[],double t,int N)
{
    int i,j;
    double memoria;
    t=t*(G*Ms)/(c*c*c);
    for(i=1;i<=N;i++)
    {
        m[i]=m[i]/Ms;
        for(j=1;j<=2;j++)
        {
            r[j][i]=r[j][i]/c;
        }
    }
    return;
}
//Funcion que devuelve las variables a la normalidad
void desescalar(double r[2][],double m[],double t,int N)
{
    int i,j;
    t=t/(G*Ms)*(c*c*c);
    for(i=1;i<=N;i++)
    {
        m[i]=m[i]/Ms;
        for(j=1;j<=2;j++)
        {
            r[j][i]=r[j][i]*c;
        }
    }
    return;
}


void Algoritmo(double r[2][],double v[2][], double a[2][],double m[],int N)








void Simulacion(double r[2][],double v[2][], double a[2][],double m[],double t,double h, double tmax)
{
    while(t<tmax)
    {
        Algoritmo(r,v,a,m,N);
        t=t+h;
    }
    return;
}

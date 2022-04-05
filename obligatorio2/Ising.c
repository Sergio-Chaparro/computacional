// Programa para la simulación de un aterial mediante
//el metodo de ising, usando la generacion de numeros aleatorios

#include "gsl_rng.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>


#define N 128
//Puntero para generar números aleatorios
gsl_rng *tau;



void InicializaRed(int red[][N]);
void Algoritmo(int red[][N],double T);
void EscribeResultados(int red[][N],FILE *f);
double EvaluaP(int red[][N],int i,int j,double T);

int main(void)
{
//Las siguientes lineas de codigo preparan
//la generacion de numeros aleatorios mediante
// el uso de funciones gls_rng_uniform por ejemplo
extern gsl_rng *tau;
int semilla=564534843;
tau=gsl_rng_alloc(gsl_rng_taus);
gsl_rng_set(tau,semilla);

//Declaro las variables que me harán falta 
//para el modelo
double T;
int red[N][N];
int iteracion;
FILE *resultados;

//Obtengo unos valores iniciales para la red
InicializaRed(red);


//Inicializo T a un valor en kelvin y hago la simulacion
T=0.5;

resultados=fopen("resultados.txt","w");
for(iteracion=0;iteracion<N*N;iteracion++)
{
    Algoritmo(red,T);
    EscribeResultados(red,resultados);
}
fclose(resultados);

}



//Función que inicializa red con valores aleatorios
//Usa el puntero tau, declarado de forma externa para generar los numeros
void InicializaRed(int red[][N])
{
    extern gsl_rng *tau;
    double x;
    int i,j;

    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            x=gsl_rng_uniform(tau);
            if(x<0.5) red[i][j]=1;
            else red[i][j]=-1;
        }
    }
    return;
}


void Algoritmo(int red[][N],double T)
{
    extern gsl_rng *tau;
    int i,j;
    double e;
    i=gsl_rng_uniform_int(tau,N-1);
    j=gsl_rng_uniform_int(tau,N-1);
    e=gsl_rng_uniform(tau);
    if(e<EvaluaP(red,i,j,T))    red[i][j]=-red[i][j];
    return;
}


void EscribeResultados(int red[][N],FILE *f)
{
    int i,j;
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            if(j!=0) 
            {
                fprintf(f,",");
            }
            fprintf(f, "%i",red[i][j]);            
        }
        fprintf(f,"\n");
    }
    fprintf(f,"\n");
    return;
}


//Funcion que calcula el valor de p para i,j
//con P=min(1,e**(-E/T)))
//con E=2s(i,j)(s(i+1,j)+s(i-1,j)+s(i,j+1)+s(i,j-1))
double EvaluaP(int red[][N],int i,int j,double T)
{
    int E;
    E=2*red[i][j]*(red[(i+1)%N][j]+red[(i-1)%N][j]+red[i][(j+1)%N]+red[i][(j-1)%N]);
    if(exp(-E/T)<1) return E;
    else return 1;
}
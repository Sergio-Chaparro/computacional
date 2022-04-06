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
int condp(int i);

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
int iteracion, pasosMC,paso,fpiter;
FILE *resultados;

//Obtengo unos valores iniciales para la red
InicializaRed(red);


//Variable que uso para representar un fotograma cada x iteraciones
fpiter=2;

//Inicializo T a un valor en kelvin y hago la simulacion
T=0.5;
pasosMC=5000;

resultados=fopen("resultados.txt","w");
for(paso=0;paso<pasosMC;paso++)
{
    for(iteracion=0;iteracion<N*N;iteracion++)
    {
        Algoritmo(red,T);        
    }    
    if((paso%fpiter)==0)
    {
        EscribeResultados(red,resultados);
    }        
}
fclose(resultados);
   
return 0;
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
    i=gsl_rng_uniform_int(tau,N);
    j=gsl_rng_uniform_int(tau,N);
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
    E=2*red[i][j]*(red[condp(i+1)][j]+red[condp(i-1)][j]+red[i][condp(j+1)]+red[i][condp(j-1)]);
    if(exp(-E/T)<1) return exp(-E/T);
    else return 1;
}


//funcion para aplicar las condiciones de contorno periodicas
//se ha declarado el valor maximo como variable global
int condp(int i)
{
    
    if(i==N) 
    {
        return 0;
    }
    else if (i==-1) 
    {
        return (N-1);
    }
    else return i;
}
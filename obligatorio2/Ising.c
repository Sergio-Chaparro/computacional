// Programa para la simulación de un aterial mediante
//el metodo de ising, usando la generacion de numeros aleatorios

#include "gsl_rng.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define N 16
//Puntero para generar números aleatorios
gsl_rng *tau;



void InicializaRed(double red[][N]);
void Algoritmo(double red[][],double T);
void EscribeResultados(double red[][N],File *f);

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
void InicializaRed(double red[][N])
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


void Algoritmo(double red[][],double T)
{
    extern gsl_rng *tau;
    int i,j;
    double p,E;


}


void EscribeResultados(double red[][N],File *f)
{
    int i,j;
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            fprintf(f, "%i,\t",red[i][j]);
            fprintf(f,"\n");
        }
        fprintf(f,"\n");
    }
}
// Programa para la simulación de un material mediante
//el metodo de ising, usando la generacion de numeros aleatorios
//Mediremos el valor promedio de ciertas magnitudes

#include "gsl_rng.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>


#define N 64

//Puntero para generar números aleatorios
gsl_rng *tau;


//Funciones para la simulacion del modelo
void InicializaRed(int red[][N],bool aleatoriedad);
void Algoritmo(int red[][N],double T);
void EscribeResultados(int red[][N],FILE *f);
double EvaluaP(int red[][N],int i,int j,double T);
int condp(int i);

//FUnciones para calcular magnitudes y escribirlas por pantalla
void EscribeMagnitudes(double Magnetizacion, double Energia,double CalorEsp, double FuncCorrelacion);
double CalcMagnetizacion(int red[][N]);  
double CalcEnergia(int red[][N]); //Falta dividir al escribir 2N²
double CalcEnergiaCua(int red[][N]);
double CalcCalorEsp(double Energia,double EnergiaCua, double T); 
void ProductoSpin(int red[][N],double PromedioSpin[][N]);
double CalcFuncCorrelacion(double PromedioSpin[][N],int medidas);   
double PromediarMedidas(double Magnitud,int Medidas);
void Inicializa(double r[][N], double valor);
double EnergiaMedia(double Energia);


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
//para el modelo T[1.5 y 3.5]
double T;
int red[N][N];
double PromedioSpin[N][N];
int iteracion, pasosMC,paso,niter1,niter2,medidas;
FILE *resultados;
bool aleatoriedad=false; //True para aleatorio, false para uniforme
//Declaro las magnitudes que promediare
double Magnetizacion=0, Energia=0, EnergiaCua=0, CalorEsp=0, FuncCorrelacion=0;

//Obtengo unos valores iniciales para la red
InicializaRed(red,aleatoriedad);
Inicializa(PromedioSpin,0.);
//Variable que uso para representar un fotograma cada x iteraciones, 1 para medidas y otro para 
//resultados
niter1=100;
niter2=50;
//Inicializo T a un valor en kelvin y hago la simulacion
T=2.5;
pasosMC=10E6;
medidas=pasosMC/niter1;



resultados=fopen("resultados.txt","w");
for(paso=0;paso<pasosMC;paso++)
{
    for(iteracion=0;iteracion<N*N;iteracion++)
    {
        Algoritmo(red,T);        
    }    
    if((paso%niter2)==0)
    {
        EscribeResultados(red,resultados);
    }
    if((paso%niter1)==0)
    {        
        Magnetizacion=Magnetizacion+CalcMagnetizacion(red);
        Energia=Energia+CalcEnergia(red);
        EnergiaCua=EnergiaCua+CalcEnergiaCua(red);
        ProductoSpin(red,PromedioSpin);
    }        
}
fclose(resultados);

//Finalizo promediando y escribiendo los valores obtenidos

//Promedio las magnitudes totales, anteriormente obtenidas
Magnetizacion=PromediarMedidas(Magnetizacion,medidas);//Magnitud final
Energia=PromediarMedidas(Energia,medidas);
EnergiaCua=PromediarMedidas(EnergiaCua,medidas);
CalorEsp=CalcCalorEsp(Energia,EnergiaCua,T);//Magnitud final
FuncCorrelacion=CalcFuncCorrelacion(PromedioSpin,medidas);//Magnitud final
Energia=EnergiaMedia(Energia);

EscribeMagnitudes(Magnetizacion,Energia,CalorEsp,FuncCorrelacion);

return 0;
}



//Función que inicializa red con valores aleatorios
//Usa el puntero tau, declarado de forma externa para generar los numeros
//Si la variable aleatoriedad es true, entonces usa un 
//patron aleatorio , si no uno uniforme que inicializa a 1
void InicializaRed(int red[][N], bool aleatoriedad)
{
    extern gsl_rng *tau;
    double x;
    int i,j;

    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            if(aleatoriedad)
            {
                x=gsl_rng_uniform(tau);
                if(x<0.5) red[i][j]=1;
                else red[i][j]=-1;
            }
            else red[i][j]=1;
            
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


//Funcion que escribe las magnitudes finales por pantalla
//Ya deben ser los resultados finales
void EscribeMagnitudes(double Magnetizacion, double Energia,double CalorEsp, double FuncCorrelacion)
{
    printf("La magnetización promedio es: %lf\n",Magnetizacion);
    printf("La energia media es: %lf\n",Energia);
    printf("El calor especifico es: %lf\n",CalorEsp);
    printf("La función de correlación: %lf\n",FuncCorrelacion);
    return;
}

//Funcion que calcula la magnetizacion total, que debe ser
//Promediada 
double CalcMagnetizacion(int red[][N])  
{
    int i,j;
    double resultado;
    resultado=0;
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            resultado=resultado+red[i][j];
        }
    }
    return resultado/(N*N);
}

//Funcion que calcula la energia total, que debe ser promediada
//Y debe dividirse por 2N²
double CalcEnergia(int red[][N])
{
    int i,j;
    double resultado;
    resultado=0;
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            resultado=resultado+red[i][j]*(red[condp(i+1)][j]+red[condp(i-1)][j]+red[i][condp(j+1)]+red[i][condp(j-1)]);
        }
    }
    return -resultado/2;
}

//Funcion que calcula la suma de energia al cuadrado, sin promediar
double CalcEnergiaCua(int red[][N])
{
    int i,j;
    double resultado;
    resultado=0;
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            resultado=resultado+pow((red[i][j]*(red[condp(i+1)][j]+red[condp(i-1)][j]+red[i][condp(j+1)]+red[i][condp(j-1)])),2);
        }
    }
    return resultado/4;
}

//Funcion que calcula el calor especifico a partir del promedio de 
//la energia y de la energia al cuadrado
double CalcCalorEsp(double Energia,double EnergiaCua, double T) 
{
    return (EnergiaCua-Energia*Energia)/(T*N*N);
}

//Funcion que calcula el promedio de el producto de un spin
// Con el spin de al lado, y lo suma, para promediar despues
void ProductoSpin(int red[][N],double PromedioSpin[][N])
{
    int i,j;
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            PromedioSpin[i][j]=PromedioSpin[i][j]+red[i][j]*red[condp(i+1)][j];
        }
    }
    return;
}

//Calcula la funcion de correlacion a partir de los productos de spin, que promedia
double CalcFuncCorrelacion(double PromedioSpin[][N],int medidas) 
{
    int i,j;
    double resultado=0;
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            resultado=resultado+PromedioSpin[i][j]/medidas;
        }
    }
    return resultado/(N*N);
}

//Promedia magnitudes con un numero de medidas
double PromediarMedidas(double Magnitud,int Medidas)
{
    return Magnitud/Medidas;
}

//Funcion que inicializa una matriz
void Inicializa(double r[][N], double valor)
{
    int i,j;
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            r[i][j]=valor;
        }
    }
    return;
}

//Funcion que devuelve la energia media, una vez promediada
//la energia 
double EnergiaMedia(double Energia)
{
    return Energia/(2*N*N);
}

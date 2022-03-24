//Programa para la simulacion del sistema solar, con 8 planetas y el sol
//Debe ser capaz de representar graficas, obtener periodos, y obtener datos de energia y momento angular

#include <math.h>
#include <stdio.h>

//defino las constantes que me harán falta para reescalar las variables 
#define G 6.67E-11
#define c 1.496E11
#define Ms 1.99E30

//Funciones que usaremos para obtener y transformar los datod
void DatosIniciales(double **r,double **v,double *m,double t,int N);
void reescalar(double **r,double **v,double *m,double t,int N);
void desescalar(double **raux,double **r,double **v,double *m,double t,int N);
void Escribedatos(double **r,double **v,double t,int N,FILE *f2);

//Funciones que usaremos dentro del algoritmod e verlet, junto con este
void Algoritmo(double **r,double **v, double **a,double **w,double *m,double t,int N, double h);
void aceleracion(double **r,double **a,double *m,int N);
void velocidad(double **r,double **v, double **a,double *m,int N,double h);
void posicion(double **r,double **v, double **a,double *m,int N,double h);
void velocidadauxiliar(double **r,double **v, double **a,double **w,double *m,int N,double h);


//Funciones que usaremos para la comprobacion de resultados
double momentoangulartotal(double **r,double **v, double **a,double *m,int N);



int main(void)
{
    //Declaro el tiempo y el paso que vamos a utilizar y las inicializo
    double h,tmax;
    int N;
    N=9;
    h=0.01;
    tmax=5;

    //Declaramos la posicion, la velocidad y la aceleracion
    //que vamos a utilizar para los N planetas
    double a[2][N],v[2][N],r[2][N], raux[2][N], w[2][N],m[N],t;

    //Declaro el fichero que usaremos para guardar los resultados
    FILE *resultados;

    //Inicializo todos los vectores
    DatosIniciales(r,v,m,t,N);    
    //Reescalamos los datod para poder tratarlos con mas facilidad
    //Usamos masas solares y la distancia tierra-sol
    reescalar(r,v,m,t,N);//A partir de aqui estaran en esas unidades



    //Podemos empezar a simular, con parametros h y tmax
    //Abrimos el fichero para escribir los datos
    resultados=fopen("resultados.txt","w");
    //Obtengo las aceleraciones iniciales
    aceleracion(r,a,m,N);

    //repito el ciclo hasta que nos interese
    while(t<tmax)
    {
        Algoritmo(r,v,a,w,m,t,N,h);
        Escribedatos(r,v,t,N,resultados);
        momentoangulartotal(r,v,a,m,N);//falta printear
    }

    fclose(resultados);
    
    
    return 0;
}


//Funcion que aporta los valores iniciales a las variables
// el archivo esta hecho de la forma
// x,y,vx,vy,m   de el elemento i
// finalmente el valor de t
void DatosIniciales(double **r,double **v,double *m,double t,int N)
{
    int i;
    FILE *f1;
    f1=fopen("datos.txt","r");
    for(i=0;i<N;i++)
    {
        fscanf(f1,"lf\t\tlf\t\tlf\t\tlf\t\tlf",&r[0][i],&r[1][i],&v[0][i],&v[1][i],&m[i]);        
    }
    fscanf(f1,"lf",&t);
    fclose(f1);
    return;
}

//Funcion que escribe datos en un fichero ya abierto 
void Escribedatos(double **r,double **v,double t,int N,FILE *f2);
{
    int i;
    for(i=0;i<N;i++)
    {
        fprintf(f2,"lf,lf",r[0][i],r[1][i]);        
    }
    //fprintf(f2,"lf",t);
    
    return;
}

//Funcion que reescala las variables para ser usadas con valores mas sencillos
void reescalar(double **r,double **v,double *m,double t,int N)
{
    int i,j;
    t=t*(G*Ms)/(c*c*c);
    for(i=0;i<N;i++)
    {
        m[i]=m[i]/Ms;
        for(j=0;j<2;j++)
        {
            r[j][i]=r[j][i]/c;
            v[j][i]=v[j][i]*c*c/(G*Ms);
        }
    }
    return;
}
//Funcion que devuelve las variables a la normalidad, en el vector raux
void desescalar(double **raux,double **r,double **v,double *m,double t,int N)
{
    int i,j;
    
    for(i=0;i<N;i++)
    {
        
        for(j=0;j<2;j++)
        {
            raux[j][i]=r[j][i]*c;
        }
    }
    return;
}




//Funcion que realiza el algoritmo de verlet, necesita la aceleracion inicial ya calculada, 
//junto con la posicion y la velocidad ya dadas
void Algoritmo(double **r,double **v, double **a,double **w,double *m,double t,int N, double h)
{
    posicion(r,v,a,m,N,h);
    velocidadauxiliar(r,v,a,w,m,N,h);
    aceleracion(r,a,m,N);
    velocidad(r,v,a,m,N,h);
    t=t+h;
    return;
}



//funcion para calcular la aceleracion de todos los cuerpos en funcion de susu posiciones
//Se usa la segunda ley de newton
void aceleracion(double **r,double **a,double *m,int N)
{
    int i,j,planeta;
    double suma;

    for(planeta=0;planeta<9;planeta++)
    {
        for(i=0;i<2;i++)
        {
        suma=0;
        for(j=0;j<N;j++)
        {
            if(planeta=!j) 
            {

                suma=suma-m[j]*(r[i][planeta]-r[i][j])/pow(sqrt(pow((r[0][planeta]-r[0][j]),2)+pow((r[1][planeta]-r[1][j]),2)),3);   
            }
        }
        a[i][planeta]=suma;    
        }
    }
    
    return;

}




void velocidad(double **r,double **v, double **a,double *m,int N,double h)
{
    int i,j;
    for(i=0;i<2;i++)
    {
        for(j=0;j<N;j++)
        {
            v[i][j]=v[i][j]+0.5*h*a[i][j];
        }
    }
}


void posicion(double **r,double **v, double **a,double *m,int N,double h)
{
    int i,j;
    for(i=0;i<2;i++)
    {
        for(j=0;j<N;j++)
        {
            r[i][j]=r[i][j]+h*v[i][j]+0.5*h*h*a[i][j];
        }
    }
}




//velocidad auxiliar que usaremos para el algoritmo de verlet 
void velocidadauxiliar(double **r,double **v, double **a,double **w,double *m,int N,double h)
{
    int i,j;
    for(i=0;i<2;i++)
    {
        for(j=0;j<N;j++)
        {
            w[i][j]=v[i][j]+0.5*h*a[i][j];
        }
    }
}




double momentoangulartotal(double **r,double **v, double **a,double *m,int N)
{
    int i,j;
    double suma=0;
    for(j=0;j<N;j++)
    {
        suma=suma+m[j]*(r[0][j]*v[1][j]-r[1][j]*v[0][j]);
    }
    return suma;
}
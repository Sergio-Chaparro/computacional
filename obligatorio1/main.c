//Programa para la simulacion del sistema solar, con 8 planetas y el sol
//Debe ser capaz de representar graficas, obtener periodos, y obtener datos de energia y momento angular


#include <math.h>
#include <stdio.h>

//defino las constantes que me harán falta para reescalar las variables 
#define G 6.67E-11
#define c 1.496E11
#define Ms 1.99E30


void DatosIniciales(double **r,double **v,double *m,double t,int N);
void reescalar(double **r,double **v,double *m,double t,int N);
void desescalar(double **r,double **v,double *m,double t,int N);
void Simulacion(double **r,double **v, double **a,double *m,double t,int N, double h, double tmax);
void aceleracion(double **r,double **a,double *m,int N);
void velocidad(double **r,double **v, double **a,double *m,int N);
void posicion(double **r,double **v, double **a,double *m,int N);
void velocidadauxiliar(double **r,double **v, double **a,double *m,int N);
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
    double a[2][N],v[2][N],r[2][N], w[2][N],m[N],t;

    //Inicializo todos los vectores
    DatosIniciales(r,v,m,t,N);
    

    //Reescalamos los datod para poder tratarlos con mas facilidad
    //Usamos masas solares y la distancia tierra-sol
    reescalar(r,v,m,t,N);//A partir de aqui estaran en esas unidades

    //Podemos empezar a simular, con parametros h y tmax
    
    Simulacion(r,v,a,m,t,N,h,tmax);

    // Podemos desescalar los datos tras la simulacion
    desescalar(r,v,m,t,N);
    
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


//Funcion que reescala las variables para ser usadas con valores mas sencillos
void reescalar(double **r,double **v,double *m,double t,int N)
{
    int i,j;
    t=t*(G*Ms)/(c*c*c);
    for(i=0;i<N;i++)
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
void desescalar(double **r,double **v,double *m,double t,int N)
{
    int i,j;
    t=t/(G*Ms)*(c*c*c);
    for(i=0;i<N;i++)
    {
        m[i]=m[i]/Ms;
        for(j=1;j<=2;j++)
        {
            r[j][i]=r[j][i]*c;
        }
    }
    return;
}





void Simulacion(double **r,double **v, double **a,double *m,double t,int N, double h, double tmax)
{
    int i,j;
    for(i=0;i<2;i++)
    {
        
    }    
    while(t<tmax)
    {
        
        t=t+h;
    }
    return;
}


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




void velocidad(double **r,double **v, double **a,double *m,int N)






void posicion(double **r,double **v, double **a,double *m,int N)





void velocidadauxiliar(double **r,double **v, double **a,double *m,int N)








double momentoangulartotal(double **r,double **v, double **a,double *m,int N)
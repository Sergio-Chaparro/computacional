//Programa para la simulacion del sistema solar, con 8 planetas y el sol
//Debe ser capaz de representar graficas, obtener periodos, y obtener datos de energia y momento angular
//En este caso coloca la tierra en el centro, mediante una funcion
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


//defino las constantes que me harán falta para reescalar las variables 
#define G 6.67E-11
#define c 1.496E11
#define Ms 1.99E30
#define N 9

//Funciones que usaremos para obtener y transformar los datod
void DatosIniciales(double r[][N],double v[][N],double m[],double t);
void reescalar(double r[][N],double v[][N],double m[],double t);
void Escribedatos(double r[][N],double v[][N],double a[][N],double m[],double t,FILE *f1,FILE *f2,FILE *f3,int reduccion,int iteraciones, int reduccion2, double tmax2);
void Tierracentro(double r[][N],double v[][N]);

//Funciones que usare para obtener periodos
void Inicializa(double r[], double valor);
void Iguala(double r1[][N],double r2[][N]);
void Escribeperiodo(FILE *f,int Vueltas[],double Periodo[]);
void Revisaperiodo(double r[][N],double raux[][N],int Vueltas[],double t,double taux[]);
void CalculaPeriodo(int Vueltas[],double taux[],double Periodo[]);

//Funciones que usaremos dentro del algoritmod e verlet, junto con este
void Algoritmo(double r[][N],double v[][N], double a[][N],double w[][N],double m[], double h);
void aceleracion(double r[][N],double a[][N],double m[]);
void velocidad(double r[][N],double v[][N],double w[][N], double a[][N],double m[],double h);
void posicion(double r[][N],double v[][N], double a[][N],double m[],double h);
void velocidadauxiliar(double r[][N],double v[][N], double a[][N],double w[][N],double m[],double h);


//Funciones que usaremos para la comprobacion de resultados
double momentoangulartotal(double r[][N],double v[][N], double a[][N],double m[]);



int main(void)
{
    //Declaro el tiempo y el paso que vamos a utilizar y las inicializo
    double h,tmax,tmax2;
    h=0.001;
    tmax=70;
    tmax2=tmax/10;

    //Defino una variable reduccion, cantidad entre
    //la cual dividire el numero de resultados obtenidos
    int reduccion=70;
    int reduccion2=7;
    int iteraciones=0;

    //Declaro variables para calcular los periodos de rotacion;
    int Planeta,Vueltas[N];
    double taux[N],Periodo[N];

    //Declaramos la posicion, la velocidad y la aceleracion
    //que vamos a utilizar para los N planetas
    double a[2][N],v[2][N],r[2][N], raux[2][N], w[2][N],m[N],t;

    //Declaro el fichero que usaremos para guardar los resultados
    FILE *resultados,*momento,*resultados2,*Periodos;

    //Inicializo todos los vectores
    DatosIniciales(r,v,m,t); 
    //Pongo la tierra en el centro
    Tierracentro(r,v);   
    //Reescalamos los datod para poder tratarlos con mas facilidad
    //Usamos masas solares y la distancia tierra-sol
    reescalar(r,v,m,t);//A partir de aqui estaran en esas unidades



    //Podemos empezar a simular, con parametros h y tmax
    //Abrimos el fichero para escribir los datos
    resultados=fopen("resultados_tierra.txt","w");
    momento=fopen("momento_tierra.txt","w");
    resultados2=fopen("resultados2_tierra.txt","w");
    Periodos=fopen("Periodos_tierra.txt","w");
    //Obtengo las aceleraciones iniciales
    aceleracion(r,a,m);
    //Inicio taux como t
    Inicializa(taux,t);
    Inicializa(Periodo,0.);
    for(Planeta=0;Planeta<N;Planeta++)
    {
        Vueltas[Planeta]=0;
    }
    //repito el ciclo hasta que nos interese
    while(t<tmax)
    {
        //raux guarda las posiciones inmediatamente anteriores
        Iguala(raux,r);
        Algoritmo(r,v,a,w,m,h);
        Escribedatos(r,v,a,m,t,resultados,momento,resultados2,reduccion,iteraciones,reduccion2,tmax2);
        //Compruebo si es el principio de un avuelta de un planeta
        Revisaperiodo(r,raux,Vueltas,t,taux);              
        t=t+h;
        iteraciones++;
    }

    //Escribo los resultado sdel periodo
    CalculaPeriodo(Vueltas,taux,Periodo);
    Escribeperiodo(Periodos,Vueltas,Periodo);

    fclose(resultados);
    fclose(momento);
    fclose(resultados2);
    fclose(Periodos);
    
    
    return 0;
}


//Funcion que aporta los valores iniciales a las variables
// el archivo esta hecho de la forma
// x,y,vx,vy,m   de el elemento i
// finalmente el valor de t
void DatosIniciales(double r[][N],double v[][N],double m[],double t)
{
    int i;
    FILE *f1;
    f1=fopen("datos.txt","r");
    for(i=0;i<N;i++)
    {
        fscanf(f1,"%lf\t\t%lf\t\t%lf\t\t%lf\t\t%lf",&r[0][i],&r[1][i],&v[0][i],&v[1][i],&m[i]);        
    }
    fscanf(f1,"%lf",&t);
    fclose(f1);
    return;
}

//Funcion que escribe datos en ficheros ya abiertos
//f1 es para las posiciones de todos los cuerpos,
//f2 es para el momento angular total del sistema
// y f3 es para representar el resto de variables 
void Escribedatos(double r[][N],double v[][N],double a[][N],double m[],double t,FILE *f1,FILE *f2,FILE *f3,int reduccion,int iteraciones, int reduccion2, double tmax2)
{
    int i;
    
    if(iteraciones%reduccion==0)
    {
        
        for(i=0;i<N;i++)
        {   
            fprintf(f1,"%lf,\t%lf\t\n",r[0][i],r[1][i]);
            
        }
        fprintf(f1,"\n");
        fprintf(f2,"%lf %E\n",t,momentoangulartotal(r,v,a,m));
       
    }
    if((iteraciones%reduccion2==0)&&(t<tmax2))
    {
        
        for(i=0;i<N;i++)
        {   
            if(i<5)
            {
                fprintf(f3,"%lf,\t%lf\t\n",r[0][i],r[1][i]); 
            }
        }
        fprintf(f3,"\n");        
    }
    return;
}

//Funcion que reescala las variables para ser usadas con valores mas sencillos
void reescalar(double r[][N],double v[][N],double m[],double t)
{
    int i,j;
    t=t*sqrt((G*Ms)/(c*c*c));
    for(i=0;i<N;i++)
    {
        m[i]=m[i]/Ms;
        for(j=0;j<2;j++)
        {
            r[j][i]=r[j][i]/c;
            v[j][i]=v[j][i]*sqrt(c/(G*Ms));
        }
    }
    return;
}


//Funcion que coloca la tierra en el centro
void Tierracentro(double r[][N],double v[][N])
{
    int i,j;
    double pos[2],vel[2];
    pos[0]=r[0][3];
    pos[1]=r[1][3];
    vel[0]=v[0][3];
    vel[1]=v[1][3];
    for(i=0;i<2;i++)
    {
        for(j=0;j<N;j++)
        {
            r[i][j]=r[i][j]-pos[i];
            //v[i][j]=v[i][j]-vel[i];
        }
    }
    return;
}


//Funcion que inicializa un vector a un valor concreto
//debe haberse definido el tamaño del vector como variable N
void Inicializa(double r[], double valor)
{
    int i;
    for(i=0;i<N;i++)
    {
        r[i]=valor;
    }
    return;
}


//Funcion que copia el contenido de r2 en r1

void Iguala(double r1[][N],double r2[][N])
{
    int i,j;
    for(i=0;i<2;i++)
    {
        for(j=0;j<N;j++)
        {
            r1[i][j]=r2[i][j];
        }
    }

}



//Funcion que comprueba si el planeta se encuentra en el principio de 
//una vuelta y en ese caso guarda  el tiempo y las vueltas que lleva
void Revisaperiodo(double r[][N],double raux[][N],int Vueltas[],double t,double taux[])
{
    
    int i;
    for(i=1;i<N;i++)
    {
        if((r[1][i]>0.)&&(raux[1][i]<0.))
        {
            taux[i]=t;
            Vueltas[i]++;
        }
    }
}

//Funcion que calcula el periodo de los planetas, apartir de las 
//vueltas que han dado y el tiempo que han tardado
void CalculaPeriodo(int Vueltas[],double taux[],double Periodo[])
{
    int i;
    for(i=1;i<N;i++)
    {
        if(Vueltas[i]!=0)
        {
            Periodo[i]=taux[i]/Vueltas[i]*sqrt((c*c*c/(G*Ms)))/24/60/60;
        }
        else Periodo[i]=0;
    }
    return;
}


//Funcion que escribe las vueltas y los periodos de los 
//Planetas, sin contar el Sol

void Escribeperiodo(FILE *f,int Vueltas[],double Periodo[])
{
    int i;
    for(i=1;i<N;i++)
    {
        fprintf(f,"Planeta %i ha dado %i vueltas, con un periodo de :%lf dias.\n",i,Vueltas[i],Periodo[i]);
    
    }
    return;
}

//Funcion que realiza el algoritmo de verlet, necesita la aceleracion inicial ya calculada, 
//junto con la posicion y la velocidad ya dadas
void Algoritmo(double r[][N],double v[][N], double a[][N],double w[][N],double m[], double h)
{
    posicion(r,v,a,m,h);
    velocidadauxiliar(r,v,a,w,m,h);
    aceleracion(r,a,m);
    velocidad(r,v,w,a,m,h);
    return;
}



//funcion para calcular la aceleracion de todos los cuerpos en funcion de susu posiciones
//Se usa la segunda ley de newton
void aceleracion(double r[][N],double a[][N],double m[])
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
                if(planeta!=j) 
                {
                    suma=suma-m[j]*(r[i][planeta]-r[i][j])/pow(sqrt(pow((r[0][planeta]-r[0][j]),2)+pow((r[1][planeta]-r[1][j]),2)),3);   
                }           
            }
        a[i][planeta]=suma;    
        }
    }
    
    return;

}



//Calcula la velocidad 
void velocidad(double r[][N],double v[][N],double w[][N], double a[][N],double m[],double h)
{
    int i,j;
    for(i=0;i<2;i++)
    {
        for(j=0;j<N;j++)
        {
            v[i][j]=w[i][j]+0.5*h*a[i][j];
        }
    }
}


void posicion(double r[][N],double v[][N], double a[][N],double m[],double h)
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
void velocidadauxiliar(double r[][N],double v[][N], double a[][N],double w[][N],double m[],double h)
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




double momentoangulartotal(double r[][N],double v[][N], double a[][N],double m[])
{
    int i,j;
    double suma=0;
    for(j=0;j<N;j++)
    {
        suma=suma+m[j]*(r[0][j]*v[1][j]-r[1][j]*v[0][j]);
    }
    return suma;
}

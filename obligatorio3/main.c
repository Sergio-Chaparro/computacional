//Programa para la resolución de la ecuación de schrödinger 
//en un pozo cuadrado

#include "complex.h"
#include <math.h>
#include <stdio.h>

#define pi 3.141592
#define N 1000
#define nciclos 250


void CalcPotencial(double lambda,fcomplex V[]);
void FuncInicial(fcomplex FuncOnda[]);
void CondContorno(fcomplex FuncOnda[]);
void CalcA(fcomplex V[], double st, fcomplex A[][3]);
void Calcb(fcomplex FuncOnda[],fcomplex b[], double st, int contador);
void CalcB(fcomplex B[], fcomplex A[][3], fcomplex b[],int contador,fcomplex alfa[]);
void CalcXi(fcomplex A[][3], fcomplex B[], fcomplex Xi[],int contador, fcomplex[]);
void Algoritmo(fcomplex B[], fcomplex b[], fcomplex A[][3], fcomplex Xi[], double st, fcomplex V[],int contador, fcomplex alfa[N],fcomplex FuncOnda[]);
void Escribe(fcomplex FuncOnda[], int contador,FILE *f,double h);


int main(void)
{
    //Defino las variables que necesito A- es la primera colunma , A0 la segunda y A+ la tercera
    int contador=0;
    fcomplex FuncOnda[N], V[N];
    fcomplex A[N][3], b[N], B[N], Xi[N],alfa[N];
    double lambda=0.3, kt,st, h=1.;
    FILE *f;

    //Calculo K y S para usarlos posteriormente
    kt=2*pi*nciclos/N;
    st=1/(4*kt*kt);

    //genero la función inicial y uso condiciones de contorno
    FuncInicial(FuncOnda);
    CondContorno(FuncOnda);
    CalcPotencial(lambda,V);
    f=fopen("resultados.txt","w");
    Escribe(FuncOnda,contador,f,h);

    //ejetuco el algoritmo    
    while (contador<nciclos-1)
    {
       Algoritmo(B,b,A,Xi,st,V,contador,alfa,FuncOnda);
       Escribe(FuncOnda,contador+1,f,h);
       contador++;
    }
    fclose(f);

    return 0;
}

//Función que calcula el potencial especificadpo, si esta en una zona es lambda*kt**2
void CalcPotencial(double lambda,fcomplex V[])
{
    int j;
    for(j=0;j<N;j++)
    {
        if((j<=3.*N/5)&&(j>=2.*N/5))
        {
            V[j]=Cmul(Complex(lambda,0.),Complex((2*pi*nciclos*1./N)*(2*pi*nciclos*1./N),0.));
        }
        else 
        {
            V[j]=Complex(0.,0.);
        }
            
    }
}

//Función que inicializa la función de onda
void FuncInicial(fcomplex FuncOnda[])
{
    int j;
    for(j=1;j<N-2;j++)
    {
        FuncOnda[j]=Cgauss(j*2*pi*nciclos*1./N,exp(-8.*(4.*j-N)*(4.*j-N)/(N*N)));        
    }
    return;
}

//Función que hace cero los extremos
void CondContorno(fcomplex FuncOnda[])
{
    FuncOnda[0]=Complex(0.,0.);
    FuncOnda[N-1]=Complex(0.,0.);
    return;
}

//Funcion que calcula los vectores A
void CalcA(fcomplex V[], double st, fcomplex A[][3])
{
    int j;
    for(j=0;j<N;j++)
    {
        A[j][0]=Complex(1.,0.);
        A[j][1]=Csub(Complex(-2,2/st),V[j]);
        A[j][2]=Complex(1.,0.);
    }
    return;
}

//FUncion que calcula el vector b a partir de la función de onda
void Calcb(fcomplex FuncOnda[],fcomplex b[], double st,int contador)
{
    int j;
    for(j=0;j<N;j++)
    {        
        b[j]=Cmul(FuncOnda[j],Complex(0,4/st));          
    }
    return;
}

//funcion que calcula el vector B
void CalcB(fcomplex B[], fcomplex A[][3], fcomplex b[],int contador,fcomplex alfa[N])
{

    int j;
    fcomplex  gamma[N];

    B[N-1]=Complex(0.,0.);
    alfa[N-1]=Complex(0,0);
    
    for(j=N-1;j>0;j--)
    {
        gamma[j]=Cdiv(Complex(1.,0),Cadd(A[j][1],Cmul(A[j][2],alfa[j])));
        B[j-1]=Cmul(gamma[j],Csub(b[j],Cmul(A[j][2],B[j])));
        alfa[j-1]=Cmul(gamma[j],Csub(Complex(0.,0.),A[j][0]));
    }
    return;
}

//funcion que calcula Xi apartir de A y B
void CalcXi(fcomplex A[][3], fcomplex B[], fcomplex Xi[],int contador,fcomplex alfa[N])
{
    int j;
    Xi[0]=Complex(0,0);
    for(j=0;j<N-2;j++)
    {
        Xi[j+1]=Cadd(Cmul(alfa[j],Xi[j]),B[j]);
    }
    return;
}

//FUncion que realiza el algoritmo especificado en el guion
void Algoritmo(fcomplex B[], fcomplex b[], fcomplex A[][3], fcomplex Xi[], double st, fcomplex V[],int contador,fcomplex alfa[N], fcomplex FuncOnda[])
{
    int j;
    CalcA(V,st,A);
    Calcb(FuncOnda,b,st,contador);
    CalcB(B,A,b,contador,alfa);
    CalcXi(A,B,Xi,contador,alfa);
    for(j=0;j<N;j++)
    {
        FuncOnda[j]=Csub(Xi[j],FuncOnda[j]);
    }
    return;    
}

//Funcion que escribe los resultados, escribe posicion, parte real y parte imaginaria
//Lo escribe en  un fichero f
void Escribe(fcomplex FuncOnda[], int contador, FILE *f, double h)
{
    int i,j;
    for(i=0;i<N;i++)
    {
        fprintf(f,"%.16lf,\t%.16lf\n",h*i,Cabs(FuncOnda[i]));        
    }
    fprintf(f,"\n");
    return;
}
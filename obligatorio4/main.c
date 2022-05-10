//Programa que simula el envio de una nave espacial a la luna

//Declaro una gran cantidad de variables que usaré
//en todo el programa y que no necesitaré cambiar
#define Mt  5.9736E24
#define Ml  0.07349E24
#define w   2.6617E-6
#define Dtl 3.844E8
#define G   6.67E-11
#define Rt  6.37816E6
#define Rl  1.7374E6
//librerias
#include <math.h>
#include <stdio.h>
#include <stdbool.h>

void RK4(double *r,double *phi,double *pr,double *pphi,double h,double *t);
double  DistLuna(double r,double phi,double t);
double DistTierra(double r,double phi,double t);
void CondIniciales(double *r,double *phi,double *pr,double *pphi,double t,double v, double theta);
void EscribeAnimacion(double r,double phi,double pr,double pphi,double t,FILE *f);

//defino las cuatro funciones que luego pasaré al algoritmo RK4
double rpunto(double r,double phi,double pr,double pphi,double t);
double phipunto(double r,double phi,double pr,double pphi,double t);
double prpunto(double r,double phi,double pr,double pphi,double t);
double pphipunto(double r,double phi,double pr,double pphi,double t);

void ReescalaVariables(double *r,double *phi,double *pr,double *pphi,double t);
//void ComprobarH(double h,double t,);


int main(void)
{
    double r,phi,pr,pphi,t=0;
    double h, m;
    double tmax,ve,v,theta;
    FILE *resultados;
    bool Comprobacion;
    int contador=0;
    //calculo la velocidad de escape
    ve=sqrt(2*G*Mt/Rt)/Dtl;
    
    //variables principales a cambiar
    v=ve*1000000;
    theta=-0.2;
    h=0.01;
    phi=-1.5;
    tmax=10E3;
    Comprobacion=false;

    CondIniciales(&r,&phi,&pr,&pphi,t,v,theta);
    ReescalaVariables(&r,&phi,&pr,&pphi,t);

    //abro los ficheros
    resultados=fopen("resultados.txt","w");
    while (t<tmax)
    {        
        if((contador%1000)==0)
        {
            EscribeAnimacion(r,phi,pr,pphi,t,resultados);            
        }        
        RK4(&r,&phi,&pr,&pphi,h,&t);
        contador++;
    }
    //cierro los ficheros
    fclose(resultados);

    return 0;
}

//funcion que realiza el algoritmo de RK4 para cuatro variables 
void RK4(double *r,double *phi,double *pr,double *pphi,double h,double *t)
{
    double K1[4],K2[4],K3[4],K4[4];
    double R,PHI,PR,PPHI;
    
    K1[0]=h*rpunto(*r,*phi,*pr,*pphi,*t);
    K1[1]=h*phipunto(*r,*phi,*pr,*pphi,*t);
    K1[2]=h*prpunto(*r,*phi,*pr,*pphi,*t);
    K1[3]=h*pphipunto(*r,*phi,*pr,*pphi,*t);

    K2[0]=h*rpunto(*r+K1[0]/2,*phi+K1[1]/2,*pr+K1[2]/2,*pphi+K1[3]/2,*t+h/2);
    K2[1]=h*phipunto(*r+K1[0]/2,*phi+K1[1]/2,*pr+K1[2]/2,*pphi+K1[3]/2,*t+h/2);
    K2[2]=h*prpunto(*r+K1[0]/2,*phi+K1[1]/2,*pr+K1[2]/2,*pphi+K1[3]/2,*t+h/2);
    K2[3]=h*pphipunto(*r+K1[0]/2,*phi+K1[1]/2,*pr+K1[2]/2,*pphi+K1[3]/2,*t+h/2);

    K3[0]=h*rpunto(*r+K2[0]/2,*phi+K2[1]/2,*pr+K2[2]/2,*pphi+K2[3]/2,*t+h/2);
    K3[1]=h*phipunto(*r+K2[0]/2,*phi+K2[1]/2,*pr+K2[2]/2,*pphi+K2[3]/2,*t+h/2);
    K3[2]=h*prpunto(*r+K2[0]/2,*phi+K2[1]/2,*pr+K2[2]/2,*pphi+K2[3]/2,*t+h/2);
    K3[3]=h*pphipunto(*r+K2[0]/2,*phi+K2[1]/2,*pr+K2[2]/2,*pphi+K2[3]/2,*t+h/2);

    
    K4[0]=h*rpunto(*r+K3[0],*phi+K3[1],*pr+K3[2],*pphi+K3[3],*t+h);
    K4[1]=h*phipunto(*r+K3[0],*phi+K3[1],*pr+K3[2],*pphi+K3[3],*t+h);
    K4[2]=h*prpunto(*r+K3[0],*phi+K3[1],*pr+K3[2],*pphi+K3[3],*t+h);
    K4[3]=h*pphipunto(*r+K3[0],*phi+K3[1],*pr+K3[2],*pphi+K3[3],*t+h);

    *r=*r+1./6*(K1[0]+2*K2[0]+2*K3[0]+K4[0]);
    *phi=*phi+1./6*(K1[1]+2*K2[1]+2*K3[1]+K4[1]);
    *pr=*pr+1./6*(K1[2]+2*K2[2]+2*K3[2]+K4[2]);
    *pphi=*pphi+1./6*(K1[3]+2*K2[3]+2*K3[3]+K4[3]);
    *t=*t+h;
    
    return;
}




void CondIniciales(double *r,double *phi,double *pr,double *pphi,double t,double v, double theta)
{
    *r=Rt;
    *pr=v*cos(theta-*phi);
    *pphi=(*r)*v*sin(theta-*phi);
    return;
}


void EscribeAnimacion(double r,double phi,double pr,double pphi,double t,FILE *f)
{
    double x,y,xl,yl;
    x=r*cos(phi);
    y=r*sin(phi);
    xl=cos(w*t);
    yl=sin(w*t);
    fprintf(f,"0,\t0\n");
    fprintf(f,"%lf,\t%lf\n",xl,yl);
    fprintf(f,"%lf,\t%lf\n\n",x,y);
    return;
}


double rpunto(double r,double phi,double pr,double pphi,double t)
{
    return pr;
}


double phipunto(double r,double phi,double pr,double pphi,double t)
{
    return (pphi/(r*r));
}


double prpunto(double r,double phi,double pr,double pphi,double t)
{
    double nu, delta, rprima;
    delta=G*Mt/(Dtl*Dtl*Dtl);
    nu=Ml/Mt;
    rprima=sqrt(1+r*r-2*r*cos(phi-w*t));
    return (pphi*pphi/(r*r*r)-delta*(1/(r*r)+nu/(rprima*rprima*rprima)*(r-cos(phi-w*t))));
}


double pphipunto(double r,double phi,double pr,double pphi,double t)
{
    double nu, delta, rprima;
    delta=G*Mt/(Dtl*Dtl*Dtl);
    nu=Ml/Mt;
    rprima=sqrt(1+r*r-2*r*cos(phi-w*t));
    return (-delta*nu*r)/(rprima*rprima*rprima)*sin(phi-w*t);
}



void ReescalaVariables(double *r,double *phi,double *pr,double *pphi,double t)
{
    *r=*r/Dtl;
    *phi=*phi;
    *pr=*pr/(Dtl);
    *pphi=*pphi/(Dtl*Dtl);
    return;
}

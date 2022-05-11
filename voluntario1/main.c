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
#define Masteroid   8.71E15
#define Vasteroid   5000
#define Rasteroid 3.46E-5
//librerias
#include <math.h>
#include <stdio.h>
#include <stdbool.h>

void RK4(double *r,double *phi,double *pr,double *pphi,double h,double *t, double angulolunainicial);

void CondIniciales(double *r,double *phi,double *pr,double *pphi,double t,double v, double theta, double m);
void EscribeAnimacion(double r,double phi,double pr,double pphi,double t,FILE *f, double angulolunainicial, double xasteroid);

//defino las cuatro funciones que luego pasaré al algoritmo RK4
double rpunto(double r,double phi,double pr,double pphi,double t, double angulolunainicial);
double phipunto(double r,double phi,double pr,double pphi,double t, double angulolunainicial);
double prpunto(double r,double phi,double pr,double pphi,double t, double angulolunainicial);
double pphipunto(double r,double phi,double pr,double pphi,double t, double angulolunainicial);

void ReescalaVariables(double *r,double *phi,double *pr,double *pphi,double t, double m);
void impulso(double *pr, double velocidad, double *energia, int numeroimpulsos, double m);

int main(void)
{
    double r,phi,pr,pphi,t=0, angulolunainicial, xasteroid;
    double h, m=235000;
    double tmax,ve,v,theta;
    FILE *resultados;
    bool Impacto=false;
    int contador=0;
    double energiausada=0, vimpulsos=660;
    int impulsos=0;
    //calculo la velocidad de escape
    ve=sqrt(2*G*Mt/Rt);
    
    //variables principales a cambiar
    v=ve*0.993;
    theta=0.46;
    h=0.01;
    phi=1.2;
    tmax=8*10E4;
    angulolunainicial=-0.85;
    
    

    CondIniciales(&r,&phi,&pr,&pphi,t,v,theta,m);
    ReescalaVariables(&r,&phi,&pr,&pphi,t,m);
    xasteroid=10.;
    //abro los ficheros
    resultados=fopen("resultados.txt","w");
   
    while ((t<tmax)&&!Impacto)
    {        
        if((contador%100000)==0)
        {
            EscribeAnimacion(r,phi,pr,pphi,t,resultados,angulolunainicial, xasteroid);                       
        }        
        RK4(&r,&phi,&pr,&pphi,h,&t,angulolunainicial);
        //el movimiento del asteroide será simplemente un movimiento acelerado en una dimensión
        xasteroid=10-Vasteroid/Dtl*t-G*Mt*Masteroid/(xasteroid*xasteroid)*t*t/2;
        contador++;
        if(((xasteroid-r*cos(phi))*(xasteroid-r*cos(phi))+r*r*sin(phi)*sin(phi))<(Rasteroid*Rasteroid)) 
        {
            Impacto=true;
            printf("velocidad radial relativa=%lf\t velocidad angular relativa=%lf\n",pr*Dtl*m,pphi*Dtl*m*Dtl);
        }
        //impulsos
        if(fabs(t-0.7*tmax)<h)
        {
            impulso(&pr,-vimpulsos,&energiausada,2,m);
            impulsos=impulsos+2;            
        } 

    }
    if(Impacto) printf("Se ha impactado con %i impulsos.\n Energia usada: %lf\n",impulsos,energiausada);

    //cierro los ficheros
    fclose(resultados);


    return 0;
}

//funcion que realiza el algoritmo de RK4 para cuatro variables 
void RK4(double *r,double *phi,double *pr,double *pphi,double h,double *t, double angulolunainicial)
{
    double K1[4],K2[4],K3[4],K4[4];
    double R,PHI,PR,PPHI;
    
    K1[0]=h*rpunto(*r,*phi,*pr,*pphi,*t,angulolunainicial);
    K1[1]=h*phipunto(*r,*phi,*pr,*pphi,*t,angulolunainicial);
    K1[2]=h*prpunto(*r,*phi,*pr,*pphi,*t,angulolunainicial);
    K1[3]=h*pphipunto(*r,*phi,*pr,*pphi,*t,angulolunainicial);

    K2[0]=h*rpunto(*r+K1[0]/2,*phi+K1[1]/2,*pr+K1[2]/2,*pphi+K1[3]/2,*t+h/2,angulolunainicial);
    K2[1]=h*phipunto(*r+K1[0]/2,*phi+K1[1]/2,*pr+K1[2]/2,*pphi+K1[3]/2,*t+h/2,angulolunainicial);
    K2[2]=h*prpunto(*r+K1[0]/2,*phi+K1[1]/2,*pr+K1[2]/2,*pphi+K1[3]/2,*t+h/2,angulolunainicial);
    K2[3]=h*pphipunto(*r+K1[0]/2,*phi+K1[1]/2,*pr+K1[2]/2,*pphi+K1[3]/2,*t+h/2,angulolunainicial);

    K3[0]=h*rpunto(*r+K2[0]/2,*phi+K2[1]/2,*pr+K2[2]/2,*pphi+K2[3]/2,*t+h/2,angulolunainicial);
    K3[1]=h*phipunto(*r+K2[0]/2,*phi+K2[1]/2,*pr+K2[2]/2,*pphi+K2[3]/2,*t+h/2,angulolunainicial);
    K3[2]=h*prpunto(*r+K2[0]/2,*phi+K2[1]/2,*pr+K2[2]/2,*pphi+K2[3]/2,*t+h/2,angulolunainicial);
    K3[3]=h*pphipunto(*r+K2[0]/2,*phi+K2[1]/2,*pr+K2[2]/2,*pphi+K2[3]/2,*t+h/2,angulolunainicial);

    
    K4[0]=h*rpunto(*r+K3[0],*phi+K3[1],*pr+K3[2],*pphi+K3[3],*t+h,angulolunainicial);
    K4[1]=h*phipunto(*r+K3[0],*phi+K3[1],*pr+K3[2],*pphi+K3[3],*t+h,angulolunainicial);
    K4[2]=h*prpunto(*r+K3[0],*phi+K3[1],*pr+K3[2],*pphi+K3[3],*t+h,angulolunainicial);
    K4[3]=h*pphipunto(*r+K3[0],*phi+K3[1],*pr+K3[2],*pphi+K3[3],*t+h,angulolunainicial);

    *r=*r+1./6*(K1[0]+2*K2[0]+2*K3[0]+K4[0]);
    *phi=*phi+1./6*(K1[1]+2*K2[1]+2*K3[1]+K4[1]);
    *pr=*pr+1./6*(K1[2]+2*K2[2]+2*K3[2]+K4[2]);
    *pphi=*pphi+1./6*(K1[3]+2*K2[3]+2*K3[3]+K4[3]);
    *t=*t+h;
    
    return;
}




void CondIniciales(double *r,double *phi,double *pr,double *pphi,double t,double v, double theta, double m)
{
    *r=Rt;
    *pr=m*v*cos(theta-*phi);
    *pphi=m*(Rt)*v*sin(theta-*phi);
    return;
}


void EscribeAnimacion(double r,double phi,double pr,double pphi,double t,FILE *f, double angulolunainicial,double xasteroid)
{
    double x,y,xl,yl;
    x=r*cos(phi);
    y=r*sin(phi);
    xl=cos(w*t+angulolunainicial);
    yl=sin(w*t+angulolunainicial);    
    
    fprintf(f,"%lf,\t%lf\n",0.,0.);
    fprintf(f,"%lf,\t%lf\n",xl,yl);
    fprintf(f,"%lf,\t%lf\n",x,y);
    fprintf(f,"%lf,\t%lf\n\n",xasteroid,0.);
    return;
}


double rpunto(double r,double phi,double pr,double pphi,double t, double angulolunainicial)
{
    return pr;
}


double phipunto(double r,double phi,double pr,double pphi,double t, double angulolunainicial)
{
    return (pphi/(r*r));
}


double prpunto(double r,double phi,double pr,double pphi,double t, double angulolunainicial)
{
    double nu, delta, rprima;
    delta=G*Mt/(Dtl*Dtl*Dtl);
    nu=Ml/Mt;
    rprima=sqrt(1+r*r-2*r*cos(phi-w*t-angulolunainicial));
    return (pphi*pphi/(r*r*r)-delta*(1/(r*r)+nu/(rprima*rprima*rprima)*(r-cos(phi-w*t-angulolunainicial))));
}


double pphipunto(double r,double phi,double pr,double pphi,double t, double angulolunainicial)
{
    double nu, delta, rprima;
    delta=G*Mt/(Dtl*Dtl*Dtl);
    nu=Ml/Mt;
    rprima=sqrt(1+r*r-2*r*cos(phi-w*t-angulolunainicial));
    return (-delta*nu*r)/(rprima*rprima*rprima)*sin(phi-w*t-angulolunainicial);
}



void ReescalaVariables(double *r,double *phi,double *pr,double *pphi,double t,double m)
{
    *r=*r/Dtl;
    *phi=*phi;
    *pr=*pr/(Dtl*m);
    *pphi=*pphi/(Dtl*Dtl*m);
    return;
}

//funcion que genera un impulso en la direccion radial, se considera
//que se pueden hacer varios impulsos a la vez
void impulso(double *pr, double velocidad, double *energia, int numeroimpulsos, double m)
{
    int i;
    for(i=0;i<numeroimpulsos;i++)
    {
        *pr=*pr+velocidad/Dtl;
        *energia=*energia+0.5*m*velocidad*velocidad;
    }
}
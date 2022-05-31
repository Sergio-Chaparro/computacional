//Programa que simula el envio de una nave espacial a destruir un asteroide

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
#define Rasteroid 10000
//librerias
#include <math.h>
#include <stdio.h>
#include <stdbool.h>

void RK4(double *r,double *phi,double *pr,double *pphi,double h,double *t, double angulolunainicial);
void RK4asteroid(double *r,double *phi,double *pr,double *pphi,double h,double *t, double angulolunainicial,double m);

void CondIniciales(double *r,double *phi,double *pr,double *pphi,double t,double v, double theta, double m);
void EscribeAnimacion(double r,double phi,double t,FILE *f, double angulolunainicial, double ra,double phia);

//defino las cuatro funciones que luego pasaré al algoritmo RK4
double rpunto(double r,double phi,double pr,double pphi,double t, double angulolunainicial);
double phipunto(double r,double phi,double pr,double pphi,double t, double angulolunainicial);
double prpunto(double r,double phi,double pr,double pphi,double t, double angulolunainicial);
double pphipunto(double r,double phi,double pr,double pphi,double t, double angulolunainicial);

void ReescalaVariables(double *r,double *phi,double *pr,double *pphi,double t, double m);
void impulsor(double *pr, double velocidad, double *energia, int numeroimpulsos, double m);
void impulsophi(double *pphi,double r,double velocidad, double *energia, int numeroimpulsos, double m);
void VariablesCometa(double *ra,double *phia,double *pra,double *pphia,double m);
double DistanciaColision(double r,double phi,double ra,double phia);

void variablesfragmentos(double *rf1,double *phif1,double *prf1,double *pphif1,double *rf2,double *phif2,double *prf2,double *pphif2,double energiabomba, double rasteroid, double phiasteroid, double prasteroid,double pphiasteroid, double m);
void escribeanimacion2(double rf1,double phif1,double rf2,double phif2,double t,FILE *f2, double angulolunainicial);



int main(void)
{
    //variables globales
    double  h,tmax,ve,t=0;
    bool Impacto=false;
    int contador=0;

    //variables de el cohete
    double r,phi,pr,pphi;
    double  m=2350000;
    double v,theta;
    double energiausada=0, vimpulsos=17;
    int MegatonesBomba=10000;
    double energiabomba=MegatonesBomba*4.18E15; //primer numero megatones
    //variables de la luna
    double  angulolunainicial;
    //variables del cometa
    double  ra,pra,phia,pphia;
    double DistanciaMinima=10*Dtl;

    //variables de fragnmentos
    double rf1,rf2,phif1,phif2,prf1,prf2,pphif1,pphif2;

    //otras variables
    FILE *resultados,*resultados2;
    double DistanciExtincion=10*Dtl;
    
    int impulsos=0;
    //calculo la velocidad de escape
    ve=sqrt(2*G*Mt/Rt);
    
    //variables principales a cambiar
    v=ve*0.994;
    theta=0.083;
    h=0.01;
    phi=0.613;
    tmax=8*10E4;
    angulolunainicial=-0.793;
    printf("Velocidad de despegue(Ve) : %lf \nAngulo de despegue : %lf\nLugar de lanzamiento : %lf\nAngulo inicial de la Luna : %lf \n",v/ve,theta,phi,angulolunainicial);
    printf("Energia del lanzamiento en megatones: %lf\n",((0.5*m*v*v)/(4.18E15)));
    
    //condiciones iniciales 
    CondIniciales(&r,&phi,&pr,&pphi,t,v,theta,m);
    ReescalaVariables(&r,&phi,&pr,&pphi,t,m);
    VariablesCometa(&ra,&phia,&pra,&pphia,m);
    //abro los ficheros
    resultados=fopen("resultados.txt","w");
   
    double Tiempoimpulso=0.578*tmax;
    int anchuraimpuslo=25;
    int ImpulsosR=0, ImpulsosPhi1=0, ImpulsosPhi2=0;;
    printf("Impulsos antes de aproximación: %i\n Tiempo del impulso central: %lf \n",anchuraimpuslo*2+1,Tiempoimpulso/tmax);
        
   //ciclo principal donde se desarrolla la simulacion
    while ((t<tmax)&&!Impacto)
    {        
        //escribo resultados
        if(((contador%100000)==0))
        {
            EscribeAnimacion(r,phi,t,resultados,angulolunainicial,ra,phia);                       
        }        
        //simulo los movimientos
        RK4(&r,&phi,&pr,&pphi,h,&t,angulolunainicial);
        RK4asteroid(&ra,&phia,&pra,&pphia,h,&t,angulolunainicial,m);
        t=t+h;
        contador++;


        //Calculo distancia minima
        if(DistanciaMinima>DistanciaColision(r,phi,ra,phia))    DistanciaMinima=DistanciaColision(r,phi,ra,phia);

        //compruebo si la nave ha impactado
        if(DistanciaColision(r,phi,ra,phia)<(Rasteroid/Dtl)) 
        {
            DistanciaMinima=DistanciaColision(r,phi,ra,phia);
            Impacto=true;
            printf("Tiempo de choque: %lf\nVelocidad radial relativa=%lf\t Velocidad angular relativa=%lf\n",t/tmax,((pr-pra*m/Masteroid)*Dtl),((pphi/r-pphia/ra*m/Masteroid)*Dtl));
            printf("Distancia explosion: %lf\n",ra);
            printf("Energia de la bomba en megatones: %i\n",MegatonesBomba);
        }

        //Momentos de impulsar la nave
        if((t>0.15*tmax)&&(r>(Rt/Dtl)))
        {
            if(fabs(t-Tiempoimpulso)<anchuraimpuslo*h)
            {
                impulsor(&pr,-vimpulsos,&energiausada,1,m);
                impulsos=impulsos+1;            
            } 
            
            if(fabs(ra-r)<(1.15*Rasteroid/Dtl))
            {
                if((fabs(pr-pra*m/Masteroid)*Dtl)>(1.1*vimpulsos))
                {
                    impulsor(&pr,-vimpulsos,&energiausada,1,m);
                    impulsos=impulsos+1;
                    ImpulsosR++;            
                }                
                           
            }

            if(fabs(phia-phi)>(atan(Rasteroid/Dtl/ra)))
            {
                if(((phia-phi)<0)&&((fabs(pphi/r-pphia/ra*m/Masteroid)*Dtl)<(5*vimpulsos)))
                {                    
                    impulsophi(&pphi,r,-vimpulsos,&energiausada,1,m);
                    impulsos=impulsos+1; 
                    ImpulsosPhi2++;   
                }
                if(((phia-phi)>0)&&((fabs(pphi/r-pphia/ra*m/Masteroid)*Dtl)<(5*vimpulsos)))
                {                    
                    impulsophi(&pphi,r,vimpulsos,&energiausada,1,m);
                    impulsos=impulsos+1;   
                    ImpulsosPhi2++; 
                }
            }

            if(fabs(phia-phi)<(atan(Rasteroid/Dtl/ra)))
            {
                if(((pphi/r-pphia/ra*m/Masteroid)*Dtl)>vimpulsos)
                {
                    impulsophi(&pphi,r,-vimpulsos,&energiausada,1,m);
                    impulsos=impulsos+1;
                    ImpulsosPhi1++;            
                }
            }
        
            
        } 

        

    }
    printf("Impulsos radiales:%i\nImpulsos Angulares: %i\nImpulsos de ajuste: %i\n",ImpulsosR,ImpulsosPhi1,ImpulsosPhi2);


    //Distancia a colision
    printf("Distancia a colision en radios del asteroide: %lf\n",DistanciaMinima*Dtl/Rasteroid);         


    //resultados finales
    printf("%i impulsos.\nEnergia usada en megatones: %lf\n",impulsos,energiausada/(4.18E15));


    //abro los ficheros
    resultados2=fopen("resultados2.txt","w");
    if(Impacto)     variablesfragmentos(&rf1,&phif1,&prf1,&pphif1,&rf2,&phif2,&prf2,&pphif2,energiabomba,ra,phia,pra,pphia,m);
    
    while((t<tmax)&&(Impacto))
    {
        //escribo resultados
        if(((contador%20000)==0))
        {
            escribeanimacion2(rf1,phif1,rf2,phif2,t,resultados2,angulolunainicial);                       
        }        
        //simulo los movimientos
        
        RK4asteroid(&rf1,&phif1,&prf1,&pphif1,h,&t,angulolunainicial,m*2);
        RK4asteroid(&rf2,&phif2,&prf2,&pphif2,h,&t,angulolunainicial,m*2);
        t=t+h;
        contador++;
        if((rf1*Dtl)<DistanciExtincion)
        {
             DistanciExtincion=rf1*Dtl;
        }
        if((rf2*Dtl)<DistanciExtincion)
        {
             DistanciExtincion=rf2*Dtl;
        }
        
    }

    printf("Distancia extincion(Km):%lf\n",DistanciExtincion/1000);
    //cierro los ficheros
    fclose(resultados);
    fclose(resultados2);
    printf("Final de ejecucion, tiempo: %lf\n",t/tmax);

    return 0;
}

//funcion que realiza el algoritmo de RK4 para cuatro variables 
void RK4(double *r,double *phi,double *pr,double *pphi,double h,double *t, double angulolunainicial)
{
    double K1[4],K2[4],K3[4],K4[4];
    
    
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
    //*t=*t+h;
    
    return;
}

//funciona en otras unidades
void RK4asteroid(double *r,double *phi,double *pr,double *pphi,double h,double *t, double angulolunainicial,double m)
{
    double K1[4],K2[4],K3[4],K4[4];

    //reescala variables para poder usar estas ecuaciones de movimiento con un m distinto    
    *pr=*pr*m/Masteroid;
    *pphi=*pphi*m/Masteroid;
    
    
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
    //*t=*t+h;
    

    //desescalo variables 
    *pr=*pr/m*Masteroid;
    *pphi=*pphi/m*Masteroid;
    return;
}

void CondIniciales(double *r,double *phi,double *pr,double *pphi,double t,double v, double theta, double m)
{
    *r=Rt;
    *pr=m*v*cos(theta-*phi);
    *pphi=m*(Rt)*v*sin(theta-*phi);
    return;
}


void EscribeAnimacion(double r,double phi,double t,FILE *f, double angulolunainicial, double ra,double phia)
{
    double x,y,xl,yl,xa,ya;
    x=r*cos(phi);
    y=r*sin(phi);
    xl=cos(w*t+angulolunainicial);
    yl=sin(w*t+angulolunainicial);  
    xa=ra*cos(phia);
    ya=ra*sin(phia);  
    
    fprintf(f,"%lf,\t%lf\n",0.,0.);
    fprintf(f,"%lf,\t%lf\n",xl,yl);
    fprintf(f,"%lf,\t%lf\n",x,y);
    fprintf(f,"%lf,\t%lf\n\n",xa,ya);
    return;
}

//funciones para resolver la atraccion gravitatoria de la tierra y la luna con el objeto en cuestion
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
void impulsor(double *pr, double velocidad, double *energia, int numeroimpulsos, double m)
{
    int i;
    for(i=0;i<numeroimpulsos;i++)
    {
        if((velocidad*(*pr*Dtl)<0)&&(fabs(velocidad)>fabs(*pr*Dtl)))
        {
            *energia=*energia+0.5*m*(*pr)*Dtl*(*pr)*Dtl;
            *energia=*energia+0.5*m*(velocidad*velocidad);
            *pr=*pr+velocidad/Dtl;            
        }
        else
        {
            *energia=*energia+0.5*m*fabs((velocidad+(*pr)*Dtl)*(velocidad+(*pr)*Dtl)-(*pr)*(*pr)*(Dtl*Dtl));
            *pr=*pr+velocidad/Dtl;            
        }
        
    }
    return;
}

void impulsophi(double *pphi,double r,double velocidad, double *energia, int numeroimpulsos, double m)
{
    int i;
    for(i=0;i<numeroimpulsos;i++)
    {
        *energia=*energia+0.5*m*fabs((velocidad+(*pphi)*Dtl/r)*(velocidad+(*pphi)*Dtl/r)-(*pphi)*(*pphi)*(Dtl*Dtl)/(r*r));
        *pphi=*pphi+velocidad/Dtl*r;
        
    }
    return;
}


//funcion que proporciona las variables del cometa ya reescaladas
void VariablesCometa(double *ra,double *phia,double *pra,double *pphia,double m)
{
    *ra=10;
    *phia=0;
    *pphia=0;
    *pra=1/(m*Dtl)*Masteroid*-Vasteroid;
    return;
}



double DistanciaColision(double r,double phi,double ra,double phia)
{
    double distancia;
    distancia=sqrt((ra*cos(phia)-r*cos(phi))*(ra*cos(phia)-r*cos(phi))+(ra*sin(phia)-r*sin(phi))*(ra*sin(phia)-r*sin(phi)));
    return distancia;
}


void variablesfragmentos(double *rf1,double *phif1,double *prf1,double *pphif1,double *rf2,double *phif2,double *prf2,double *pphif2,double energiabomba, double rasteroid, double phiasteroid, double prasteroid,double pphiasteroid, double m)
{
    double v;    
    v=sqrt(2*energiabomba/Masteroid);    
    *rf1=rasteroid;
    *rf2=rasteroid;
    *phif1=phiasteroid;
    *phif2=phiasteroid;
    *prf1=prasteroid;
    *prf2=prasteroid;
    *pphif1=v*Masteroid/(2*m/rasteroid*Dtl);
    *pphif2=-v*Masteroid/(2*m/rasteroid*Dtl);
    
    printf("Velocidad fragmentos: %lf\n",v);
    return;
}



void escribeanimacion2(double rf1,double phif1,double rf2,double phif2,double t,FILE *f2, double angulolunainicial)
{
    double x1,y1,x2,y2,xl,yl;
    x1=rf1*cos(phif1);
    y1=rf1*sin(phif1);
    xl=cos(w*t+angulolunainicial);
    yl=sin(w*t+angulolunainicial);  
    x2=rf2*cos(phif2);
    y2=rf2*sin(phif2);  
    
    fprintf(f2,"%lf,\t%lf\n",0.,0.);
    fprintf(f2,"%lf,\t%lf\n",xl,yl);
    fprintf(f2,"%lf,\t%lf\n",x1,y1);
    fprintf(f2,"%lf,\t%lf\n\n",x2,y2);
    return;

}


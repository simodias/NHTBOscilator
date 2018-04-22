#include <iostream>
#include <fstream>
#include <math.h>
#include <conio.h>

//------------------------------------------------------------------------------
//FUNKCJE MATEMATYCZNE
//------------------------------------------------------------------------------

//dzielenie---------------------------------------------------------------------
double dziel (double D1, double D2)
{
double iloraz;
iloraz=D1*pow(D2,-1);
return iloraz;       
}
//pierwiastkowanie--------------------------------------------------------------
double pierw (double D1)
{
double pierwiastek;
pierwiastek=pow(D1,0.5);
return pierwiastek;       
}
//kwadratowanie-----------------------------------------------------------------
double kw (double D1)
{
double kwadrat;
kwadrat=pow(D1,2);
return kwadrat;       
}
//zaokraglanie------------------------------------------------------------------
double przybliz (double D1)
{
double zaokraglenie;
zaokraglenie=int(D1+0.5);
return zaokraglenie;
} 

//------------------------------------------------------------------------------
//DEKLARACJE ZMIENNYCH
//------------------------------------------------------------------------------

//FUNKCJE FIZYCZNE DLA KLASY OSCYLATOR


//------------------------------------------------------------------------------
//KLASA OSCYLATOR
//------------------------------------------------------------------------------
class Oscl
{
      
//------------------------------------------------------------------------------
      public:
             
//liczby opisujace stan czastki------------------------------------------------- 
             double x, y, z, x0, y0, z0, xp, yp, zp;
             double px, py, pz, px0, py0, pz0, ppx, ppy, ppz;
             double fx, fy, fz, fxp, fyp, fzp;
             double m, k, CKFp, CKF; //calkowity kwadrat sily dzialajacy na i-ta czastke

//liczby opisujace stan calego ukladu niezbedne do posrednich obliczen----------
             double UC, UCp, D2FC, D2FCp, dzeta, dzeta0, dp, s, s0, sp;

//metody------------------------------------------------------------------------             


//USTAWIENIA POCZATKOWE---------------------------------------------------------
void ustaw(double pozX, double pozY, double pozZ, double prX, double prY, double prZ, double Masa, double stalaK)
{       
       k=stalaK;
       m=Masa;        
       x0=pozX;
       y0=pozY;
       z0=pozZ;
    
       px0=prX*m;
       py0=prY*m;
       pz0=prZ*m;
       
       dzeta0=0;
       s0=1;
    
       x=x0;
       y=y0;
       z=z0;
       px=px0;
       py=py0;
       pz=pz0;
       
       dzeta=dzeta0;
       s=s0;

} 
            
//STAN CZASTKI------------------------------------------------------------------             

void stan()
{
double dx, dy, dz;   
double EU;
double Fx, Fy;
double Fz, F2;
double F2calk, Fcalk, EUsum;
double Fxx, Fyy, Fzz;
double stalaK;

stalaK=k;
   
dx=x;
dy=y;
dz=z; 
    
EU=UP(dx,dy,dz,stalaK);				
					
Fx=FX(dx,stalaK);
Fy=FY(dy,stalaK);
Fz=FZ(dz,stalaK);
                    
F2=D2F(stalaK); 
					 
Fxx=0.0;
Fyy=0.0;
Fzz=0.0;
F2calk=0.0;
Fcalk=0.0; 
EUsum=0.0;

Fxx=Fx;
Fyy=Fy;
Fzz=Fz;
    
F2calk=F2;
EUsum=EU;
    
Fcalk=Fxx*Fxx+Fyy*Fyy+Fzz*Fzz;   
  
fx=Fxx;
fy=Fyy;
fz=Fzz;
  
D2FC=F2calk;
CKF=Fcalk;
UC=EUsum;   
     
}             
//STAN PROBNY CZASTKI-----------------------------------------------------------

void stanP()
{
double dx, dy, dz;   
double EU;
double Fx, Fy;
double Fz, F2;
double F2calk, Fcalk, EUsum;
double Fxx, Fyy, Fzz;
double stalaK;

stalaK=k;

dx=xp;
dy=yp;
dz=zp; 
    
EU=UP(dx,dy,dz,stalaK);				
					
Fx=FX(dx,stalaK);
Fy=FY(dy,stalaK);
Fz=FZ(dz,stalaK);
                    
F2=D2F(stalaK); 
				    
Fxx=0.0;
Fyy=0.0;
Fzz=0.0;
F2calk=0.0;
Fcalk=0.0; 
EUsum=0.0;

Fxx=Fx;
Fyy=Fy;
Fzz=Fz;
    
F2calk=F2;
EUsum=EU;
    
    
Fcalk=Fxx*Fxx+Fyy*Fyy+Fzz*Fzz;   
  
fxp=Fxx;
fyp=Fyy;
fzp=Fzz;
  
D2FCp=F2calk;
CKFp=Fcalk;
UCp=EUsum;   
    
} 
//------------------------------------------------------------------------------
void NewtonRK4(double KrokCzasowy)
{
       
double k1X, k2X, k3X, k4X;
double k1Y, k2Y, k3Y, k4Y;
double k1Z, k2Z, k3Z, k4Z;
       
double w1X, w2X, w3X, w4X;
double w1Y, w2Y, w3Y, w4Y;
double w1Z, w2Z, w3Z, w4Z;       
       
       
//BLOK 1------------------------------------------------------------------------       
      
       
       x0=x;    
       y0=y;
       z0=z; 
      
       px0=px;    
       py0=py;
       pz0=pz;
      
       xp=x0;    
       yp=y0;
       zp=z0; 
      
       ppx=px0;    
       ppy=py0;
       ppz=pz0;
       
      
       stanP();
      
                   
       w1X=dziel(ppx,m);
       w1Y=dziel(ppy,m);
       w1Z=dziel(ppz,m);
       k1X=fxp;
       k1Y=fyp;
       k1Z=fzp;
       
       
//BLOK 2------------------------------------------------------------------------       
       
       
       xp=x0+0.5*KrokCzasowy*w1X;
       yp=y0+0.5*KrokCzasowy*w1Y;
       zp=z0+0.5*KrokCzasowy*w1Z;
       
       ppx=px0+0.5*KrokCzasowy*k1X;
       ppy=py0+0.5*KrokCzasowy*k1Y;
       ppz=pz0+0.5*KrokCzasowy*k1Z;
       

       stanP();
       
                  
       w2X=dziel(ppx,m);
       w2Y=dziel(ppy,m);
       w2Z=dziel(ppz,m);
       k2X=fxp;
       k2Y=fyp;
       k2Z=fzp;
       
       
//BLOK 3------------------------------------------------------------------------

       
       xp=x0+0.5*KrokCzasowy*w2X;
       yp=y0+0.5*KrokCzasowy*w2Y;
       zp=z0+0.5*KrokCzasowy*w2Z;
       
       ppx=px0+0.5*KrokCzasowy*k2X;
       ppy=py0+0.5*KrokCzasowy*k2Y;
       ppz=pz0+0.5*KrokCzasowy*k2Z;
       
       
       stanP();

                 
       w3X=dziel(ppx,m);
       w3Y=dziel(ppy,m);
       w3Z=dziel(ppz,m);
       k3X=fxp;
       k3Y=fyp;
       k3Z=fzp;
       
       
//BLOK 4------------------------------------------------------------------------  
       
       
       xp=x0+KrokCzasowy*w3X;
       yp=y0+KrokCzasowy*w3Y;
       zp=z0+KrokCzasowy*w3Z;
       
       ppx=px0+KrokCzasowy*k3X;
       ppy=py0+KrokCzasowy*k3Y;
       ppz=pz0+KrokCzasowy*k3Z;
       
       stanP();
       
                   
       w4X=dziel(ppx,m);
       w4Y=dziel(ppy,m);
       w4Z=dziel(ppz,m);
       k4X=fxp;
       k4Y=fyp;
       k4Z=fzp;
       
       
//NOWE POZYCJE I PREDKOSCI------------------------------------------------------
       
       
       px=px0+dziel(KrokCzasowy,6)*(k1X+2*k2X+2*k3X+k4X);
       py=py0+dziel(KrokCzasowy,6)*(k1Y+2*k2Y+2*k3Y+k4Y);
       pz=pz0+dziel(KrokCzasowy,6)*(k1Z+2*k2Z+2*k3Z+k4Z);
       
       x=x0+dziel(KrokCzasowy,6)*(w1X+2*w2X+2*w3X+w4X);
       y=y0+dziel(KrokCzasowy,6)*(w1Y+2*w2Y+2*w3Y+w4Y);
       z=z0+dziel(KrokCzasowy,6)*(w1Z+2*w2Z+2*w3Z+w4Z);

       
}
//ROWNANIA NOSEGO-HOOVERA-------------------------------------------------------

void NoseHooverRK4(double Q, double temp, double g, double kB, double KrokCzasowy)
{
double PP=0;       
       
double kd1, ks1, kd2, ks2, kd3, ks3, kd4, ks4;
      
double k1X, k2X, k3X, k4X;
double k1Y, k2Y, k3Y, k4Y;
double k1Z, k2Z, k3Z, k4Z;
       
double w1X, w2X, w3X, w4X;
double w1Y, w2Y, w3Y, w4Y;
double w1Z, w2Z, w3Z, w4Z;       
       
       
//BLOK 1------------------------------------------------------------------------       
      
       
       x0=x;    
       y0=y;
       z0=z; 
      
       px0=px;    
       py0=py;
       pz0=pz;
      
       xp=x0;    
       yp=y0;
       zp=z0; 
      
       ppx=px0;    
       ppy=py0;
       ppz=pz0;
       
       dzeta0=dzeta;
       s0=s;
       
       dp=dzeta0;
       sp=s0;
       
       
       
      
       stanP();
       PP=EnKinP(); 
      
                 
       w1X=dziel(ppx,m);
       w1Y=dziel(ppy,m);
       w1Z=dziel(ppz,m);
       k1X=fxp-dp*ppx;
       k1Y=fyp-dp*ppy;
       k1Z=fzp-dp*ppz;
       
       
   
   
kd1=dziel(1,Q)*(2*PP-g*kB*temp);
ks1=sp*dp;
//BLOK 2------------------------------------------------------------------------       
       
       
       xp=x0+0.5*KrokCzasowy*w1X;
       yp=y0+0.5*KrokCzasowy*w1Y;
       zp=z0+0.5*KrokCzasowy*w1Z;
       
       ppx=px0+0.5*KrokCzasowy*k1X;
       ppy=py0+0.5*KrokCzasowy*k1Y;
       ppz=pz0+0.5*KrokCzasowy*k1Z;
       
       dp=dzeta0+0.5*KrokCzasowy*kd1;
       sp=s0+0.5*KrokCzasowy*ks1;
       
       


       stanP();
       PP=EnKinP(); 
       
                 
       w2X=dziel(ppx,m);
       w2Y=dziel(ppy,m);
       w2Z=dziel(ppz,m);
       k2X=fxp-dp*ppx;
       k2Y=fyp-dp*ppy;
       k2Z=fzp-dp*ppz;
       
kd2=dziel(1,Q)*(2*PP-g*kB*temp);
ks2=sp*dp;
       
//BLOK 3------------------------------------------------------------------------

       xp=x0+0.5*KrokCzasowy*w2X;
       yp=y0+0.5*KrokCzasowy*w2Y;
       zp=z0+0.5*KrokCzasowy*w2Z;
       
       ppx=px0+0.5*KrokCzasowy*k2X;
       ppy=py0+0.5*KrokCzasowy*k2Y;
       ppz=pz0+0.5*KrokCzasowy*k2Z;
       
       dp=dzeta0+0.5*KrokCzasowy*kd2;
       sp=s0+0.5*KrokCzasowy*ks2;
       

       stanP();
       PP=EnKinP();
          
       w3X=dziel(ppx,m);
       w3Y=dziel(ppy,m);
       w3Z=dziel(ppz,m);
       k3X=fxp-dp*ppx;
       k3Y=fyp-dp*ppy;
       k3Z=fzp-dp*ppz;
       
kd3=dziel(1,Q)*(2*PP-g*kB*temp);
ks3=sp*dp;
       
//BLOK 4------------------------------------------------------------------------  
       
       
       xp=x0+KrokCzasowy*w3X;
       yp=y0+KrokCzasowy*w3Y;
       zp=z0+KrokCzasowy*w3Z;
       
       ppx=px0+KrokCzasowy*k3X;
       ppy=py0+KrokCzasowy*k3Y;
       ppz=pz0+KrokCzasowy*k3Z;
       
       dp=dzeta0+KrokCzasowy*kd3;
       sp=s0+KrokCzasowy*ks3;
       
       

       stanP();
       PP=EnKinP();
       
                 
       w4X=dziel(ppx,m);
       w4Y=dziel(ppy,m);
       w4Z=dziel(ppz,m);
       k4X=fxp-dp*ppx;
       k4Y=fyp-dp*ppy;
       k4Z=fzp-dp*ppz;
       
kd4=dziel(1,Q)*(2*PP-g*kB*temp);
ks4=sp*dp;
       
//NOWE POZYCJE I PREDKOSCI------------------------------------------------------
       
       
       px=px0+dziel(KrokCzasowy,6)*(k1X+2*k2X+2*k3X+k4X);
       py=py0+dziel(KrokCzasowy,6)*(k1Y+2*k2Y+2*k3Y+k4Y);
       pz=pz0+dziel(KrokCzasowy,6)*(k1Z+2*k2Z+2*k3Z+k4Z);
       
       x=x0+dziel(KrokCzasowy,6)*(w1X+2*w2X+2*w3X+w4X);
       y=y0+dziel(KrokCzasowy,6)*(w1Y+2*w2Y+2*w3Y+w4Y);
       z=z0+dziel(KrokCzasowy,6)*(w1Z+2*w2Z+2*w3Z+w4Z);
       
       dzeta=dzeta0+dziel(KrokCzasowy,6)*(kd1+2*kd2+2*kd3+kd4);
       s=s0+dziel(KrokCzasowy,6)*(ks1+2*ks2+2*ks3+ks4);
       
       
       
}
//ROWNANIA TRAVISA-BRAGA--------------------------------------------------------

double TravisBragaRK4(double Q, double temp, double kB, double KrokCzasowy)
{      
       
double kd1, ks1, kd2, ks2, kd3, ks3, kd4, ks4;
double SumCKF;
      
double k1X, k2X, k3X, k4X;
double k1Y, k2Y, k3Y, k4Y;
double k1Z, k2Z, k3Z, k4Z;
       
double w1X, w2X, w3X, w4X;
double w1Y, w2Y, w3Y, w4Y;
double w1Z, w2Z, w3Z, w4Z;       
       
       
//BLOK 1------------------------------------------------------------------------       
      
      
       x0=x;    
       y0=y;
       z0=z; 
      
       px0=px;    
       py0=py;
       pz0=pz;
      
       xp=x0;    
       yp=y0;
       zp=z0; 
      
       ppx=px0;    
       ppy=py0;
       ppz=pz0;
       
       dzeta0=dzeta;
       s0=s;
       
       dp=dzeta0;
       sp=s0;
       
       
       
       stanP();
      
      
       SumCKF=0;
                  
       w1X=dziel(ppx,m)+dp*fxp;
       w1Y=dziel(ppy,m)+dp*fyp;
       w1Z=dziel(ppz,m)+dp*fzp;
       k1X=fxp;
       k1Y=fyp;
       k1Z=fzp;
       SumCKF=SumCKF+CKFp;
      
       
     
kd1=dziel(1,Q)*(SumCKF-temp*D2FCp);
ks1=sp*dp*D2FCp;

//BLOK 2------------------------------------------------------------------------       
       
      
       xp=x0+0.5*KrokCzasowy*w1X;
       yp=y0+0.5*KrokCzasowy*w1Y;
       zp=z0+0.5*KrokCzasowy*w1Z;
       
       ppx=px0+0.5*KrokCzasowy*k1X;
       ppy=py0+0.5*KrokCzasowy*k1Y;
       ppz=pz0+0.5*KrokCzasowy*k1Z;
       
       dp=dzeta0+0.5*KrokCzasowy*kd1;
       sp=s0+0.5*KrokCzasowy*ks1;
       
       

       stanP();
       
       
       SumCKF=0;
                   
       w2X=dziel(ppx,m)+dp*fxp;
       w2Y=dziel(ppy,m)+dp*fyp;
       w2Z=dziel(ppz,m)+dp*fzp;
       k2X=fxp;
       k2Y=fyp;
       k2Z=fzp;
       SumCKF=SumCKF+CKFp;
       
       
kd2=dziel(1,Q)*(SumCKF-temp*D2FCp);
ks2=sp*dp*D2FCp;       
     
       
//BLOK 3------------------------------------------------------------------------

       
       xp=x0+0.5*KrokCzasowy*w2X;
       yp=y0+0.5*KrokCzasowy*w2Y;
       zp=z0+0.5*KrokCzasowy*w2Z;
       
       ppx=px0+0.5*KrokCzasowy*k2X;
       ppy=py0+0.5*KrokCzasowy*k2Y;
       ppz=pz0+0.5*KrokCzasowy*k2Z;
       
       dp=dzeta0+0.5*KrokCzasowy*kd2;
       sp=s0+0.5*KrokCzasowy*ks2;
       


       stanP();


       SumCKF=0;
                
       w3X=dziel(ppx,m)+dp*fxp;
       w3Y=dziel(ppy,m)+dp*fyp;
       w3Z=dziel(ppz,m)+dp*fzp;
       k3X=fxp;
       k3Y=fyp;
       k3Z=fzp;
       SumCKF=SumCKF+CKFp;
       
       
kd3=dziel(1,Q)*(SumCKF-temp*D2FCp);
ks3=sp*dp*D2FCp;       
       
//BLOK 4------------------------------------------------------------------------  
       
       
       xp=x0+KrokCzasowy*w3X;
       yp=y0+KrokCzasowy*w3Y;
       zp=z0+KrokCzasowy*w3Z;
       
       ppx=px0+KrokCzasowy*k3X;
       ppy=py0+KrokCzasowy*k3Y;
       ppz=pz0+KrokCzasowy*k3Z;
       
       dp=dzeta0+KrokCzasowy*kd3;
       sp=s0+KrokCzasowy*ks3;
       
       

       stanP();
       
       
       SumCKF=0;
                   
       w4X=dziel(ppx,m)+dp*fxp;
       w4Y=dziel(ppy,m)+dp*fyp;
       w4Z=dziel(ppz,m)+dp*fzp;
       k4X=fxp;
       k4Y=fyp;
       k4Z=fzp;
       SumCKF=SumCKF+CKFp;
       
       
kd4=dziel(1,Q)*(SumCKF-temp*D2FCp);
ks4=sp*dp*D2FCp; 
       
//NOWE POZYCJE I PREDKOSCI------------------------------------------------------
       
      
       px=px0+dziel(KrokCzasowy,6)*(k1X+2*k2X+2*k3X+k4X);
       py=py0+dziel(KrokCzasowy,6)*(k1Y+2*k2Y+2*k3Y+k4Y);
       pz=pz0+dziel(KrokCzasowy,6)*(k1Z+2*k2Z+2*k3Z+k4Z);
       
       x=x0+dziel(KrokCzasowy,6)*(w1X+2*w2X+2*w3X+w4X);
       y=y0+dziel(KrokCzasowy,6)*(w1Y+2*w2Y+2*w3Y+w4Y);
       z=z0+dziel(KrokCzasowy,6)*(w1Z+2*w2Z+2*w3Z+w4Z);
       
       dzeta=dzeta0+dziel(KrokCzasowy,6)*(kd1+2*kd2+2*kd3+kd4);
       s=s0+dziel(KrokCzasowy,6)*(ks1+2*ks2+2*ks3+ks4);
       
}       


//------------------------------------------------------------------------------
//DODATKOWE FUNKCJE NA KLASIE OSCYLATOR   


//ENERGIA POTENCJALNA UKLADU CZASTEK--------------------------------------------
double energiaP()
{
double potencjalna;
potencjalna=0;

potencjalna=UC;

return potencjalna;        
}
//ENERGIA KINETYCZNA UKLADU CZASTEK---------------------------------------------
double energiaK()
{
double kinetyczna; 
kinetyczna=0;

kinetyczna=dziel(kw(px)+kw(py)+kw(pz),2*m);

return kinetyczna;        
}
//PROBNA ENERGIA KINETYCZNA UKLADU CZASTEK--------------------------------------
double EnKinP()
{
double kinP; 
kinP=0;

kinP=dziel(kw(ppx)+kw(ppy)+kw(ppz),2*m);

return kinP;        
}
//ENERGIA CALKOWITA UKLADU CZASTEK---------------------------------------------
double energiaC()
{
double calkowita;
calkowita=0; 
calkowita=dziel(kw(px)+kw(py)+kw(pz),2*m)+UC;
return calkowita;       
} 



private:
//energia potencjalna-----------------------------------------------------------
double UP (double dx, double dy, double dz, double k)
{
double r;
double energia;
r=dx*dx+dy*dy+dz*dz;
              
        energia=0.5*k*r;				              
    
return energia;           
}
//sila X------------------------------------------------------------------------
double FX (double dx, double k)
{
double silaX;
        
        silaX=-k*dx;            
    
return silaX;           
}
//sila Y------------------------------------------------------------------------
double FY (double dy, double k)
{
double silaY;
        
        silaY=-k*dy;          
    
return silaY;           
}
//sila Z------------------------------------------------------------------------
double FZ (double dz, double k)
{
double silaZ;            
        
        silaZ=-k*dz;             
    
return silaZ;           
}
//druga pochodna potencjalu-----------------------------------------------------
double D2F (double k)
{
double Dsila;
                                  
        Dsila=k;                
    
return Dsila;           
}

        
};

//------------------------------------------------------------------------------      
//FUNKCJA GLOWNA
//------------------------------------------------------------------------------      
using namespace std;	
main()
{
//------------------------------------------------------------------------------
ofstream pozycje;
ofstream pedy;
ofstream energie;

pozycje.open("pozycje");
if(!pozycje)return 0;

pedy.open("pedy");
if(!pedy)return 0;

energie.open("energie");
if(!energie)return 0;      
//------------------------------------------------------------------------------      
Oscl oscylator;

double EnP, EnK, EnC, calka, ZmSter;

double m, k, dt, x0, y0, z0, vx0, vy0, vz0;    
double Q, g, kB, temp;
int petla=0, skok=100;

x0=2; y0=0; z0=0; 
vx0=0; vy0=0; vz0=0;      

dt=0.001;
temp=1;
m=1;      
k=1;
Q=2;
g=1;
kB=1;

oscylator.ustaw(x0, y0, z0, vx0, vy0, vz0, m, k);      

cout<<oscylator.x<<"   "<<oscylator.y<<"   "<<oscylator.z<<endl;   
oscylator.stan();
cout<<oscylator.energiaC()<<"   "<<oscylator.energiaP()<<"   "<<oscylator.energiaK()<<endl;    
for(int i=0;i<500000;i++)
{
//oscylator.NewtonRK4(dt);
oscylator.NoseHooverRK4(Q, temp, g, kB, dt);
//oscylator.TravisBragaRK4(Q, temp, kB, dt);
oscylator.stan();

ZmSter=oscylator.dzeta;

EnP=oscylator.energiaP();
EnK=oscylator.energiaK();
EnC=oscylator.energiaC();
//calka=EnC;//N
calka=EnC+0.5*Q*oscylator.dzeta*oscylator.dzeta+g*kB*temp*log(oscylator.s);//NH
//calka=EnC+0.5*Q*oscylator.dzeta*oscylator.dzeta+kB*temp*log(oscylator.s);//TB

cout.precision(16);
//cout<<calka<<"   "<<EnC<<"   "<<oscylator.energiaP()<<"   "<<oscylator.energiaK()<<endl;   
cout<<calka<<endl;
petla=i-1;
if(petla%skok==0)
{
pozycje.precision(16);
pozycje<<i+1<<"   "<<oscylator.x<<"   "<<oscylator.y<<"   "<<oscylator.z<<endl;   
pozycje.flush();

pedy.precision(16);
pedy<<i+1<<"   "<<oscylator.px<<"   "<<oscylator.py<<"   "<<oscylator.pz<<endl;   
pedy.flush();

energie.precision(16);
energie<<i+1<<"   "<<calka<<"   "<<EnC<<"   "<<EnP<<"   "<<EnK<<"   "<<ZmSter<<endl;   
energie.flush();
}
}
cout<<"KONIEC"<<endl;
   

getch();      
}

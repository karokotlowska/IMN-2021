#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void jawnaEulera(FILE *fp,double y,double lambda,double range[],double deltaT){
 
  for(double i=range[0];i<=range[1];i+=deltaT){
    y=y+deltaT*lambda*y;
    fprintf(fp,"%g %g %g\n",i + deltaT,y,y-exp(lambda*(i+deltaT))); 
  }
  fprintf(fp,"\n");
}


void jawnaRK2(FILE *fp,double y,double lambda,double range[],double deltaT){
  double k1,k2;
  for(double i=range[0];i<=range[1];i+=deltaT){
    k1=lambda*y;
    k2=lambda*(y+deltaT*k1);
    y=y+deltaT/2*(k1+k2);
    fprintf(fp,"%g %g %g\n",i + deltaT,y,y-exp(lambda*(i+deltaT)));
  }
    fprintf(fp,"\n");
}

void jawnaRK4(FILE *fp,double y,double lambda,double range[],double deltaT){
  double k1,k2,k3,k4;
  for(double i=range[0];i<=range[1];i+=deltaT){
    k1=lambda*y;
    k2=lambda*(y+deltaT/2*k1);
    k3=lambda*(y+deltaT/2*k2);
    k4=lambda*(y+deltaT*k3);
    y=y+deltaT/6*(k1+2*k2+2*k3+k4);
    fprintf(fp,"%g %g %g\n",i + deltaT,y,y-exp(lambda*(i+deltaT)));
  }
  fprintf(fp,"\n");
}


double V(double OmegaV,double t){ return 10.*sin(OmegaV * t);};
double f (double I,double k_i){return I + k_i;};
double g (double R,double L,double C,double omegaV,double t,double Q, double I){return V(omegaV,t)/L-(1.0/(L*C))*Q-(R/L)*I;};

void RRZ(FILE *fp,double R,double L,double C,double T0, double omega0,double range[],double Q0,double I0,double omegaV,double deltaT){
  double k1Q,k2Q,k3Q,k4Q;
  double k1I,k2I,k3I,k4I;
  for(double i=range[0];i<=range[1];i+=deltaT){
    k1Q=I0;
    k1I=g(R,L,C,omegaV,i,Q0,I0);
    k2Q=f(I0,deltaT/2*k1I);
    k2I=g(R,L,C,omegaV,i,Q0+deltaT/2*k1Q,I0+deltaT/2*k1Q);
    k3Q=f(I0,deltaT/2*k2I);
    k3I=g(R,L,C,omegaV,i,Q0+deltaT/2*k3Q,I0+deltaT/2*k3Q);
    k4Q=f(I0,deltaT/2*k3I);
    k4I=g(R,L,C,omegaV,i,Q0+deltaT*k3Q,I0+deltaT*k3Q);

    Q0=Q0+deltaT/6*(k1Q+2*k2Q+2*k3Q+k4Q);
    I0=I0+deltaT/6*(k1I+2*k2I+2*k3I+k4I);

    fprintf(fp,"%g %g %g\n",i + deltaT,I0,Q0);
  }
  fprintf(fp,"\n");
}


int main(void) {
  //////////pierwsza czesc
  double t[]={0.01,0.1,1.};
  double y0=1.;
  double lambda=-1.;
  double range[]={0.,5.};

  FILE *fp=fopen("out.dat", "w");
  FILE *fp2=fopen("out2.dat", "w");
  FILE *fp3=fopen("out3.dat", "w");

  for(int i=0;i<sizeof(t)/sizeof(t[0]);i++){
    jawnaEulera(fp,y0,lambda,range,t[i]);
    jawnaRK2(fp2,y0,lambda,range,t[i]);
    jawnaRK4(fp3,y0,lambda,range,t[i]);
  }
  fclose (fp);
  fclose (fp2);
  fclose(fp3);

  //////////druga czesc
  double deltaT = 0.0001;
  double R = 100.;
  double L = 0.1;
  double C = 0.001;
  double omega0 = 1./(sqrt(L*C));
  double T0 = (2. * M_PI)/omega0;
  double range2[]={0.,4.*T0};
  double Q0=0.;
  double I0=0.;
  double t2[]={0.5*omega0,0.8*omega0,1.*omega0,1.2*omega0};

  FILE *fp4=fopen("out4.dat", "w");

  for(int i=0;i<sizeof(t2)/sizeof(t2[0]);i++){
    RRZ(fp4,R,L,C,T0,omega0,range2,Q0,I0,t2[i],deltaT);
  }
  fclose(fp4);

}

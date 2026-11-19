/*  Two-component dark matter of the chiral U(1)_X model, arXiv:2506.04718.
    CDM[1] = ~X  (Dirac,    Z2  odd),
    CDM[2] = ~~M (Majorana, Z2' odd).
    Usage: ./main data1.par
*/

#define MASSES_INFO
#define OMEGA
#define CDM_NUCLEON
#define CDM_NUCLEUS
#define DECAYS

#define CLEAN

#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"
#include"lib/pmodel.h"

static void nucleonCS(int k, double *csSIp, double *csSDp, double *csSDn)
{ double pA0[2],pA5[2],nA0[2],nA5[2];
  double Nmass=0.939;
  double SCcoeff;

  nucleonAmplitudes(CDM[k], pA0,pA5,nA0,nA5);
  SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*McdmN[k]/(Nmass+McdmN[k]),2.);
  /* average over DM and anti-DM */
  *csSIp=  SCcoeff*(pA0[0]*pA0[0]+pA0[1]*pA0[1])/2;
  *csSDp=3*SCcoeff*(pA5[0]*pA5[0]+pA5[1]*pA5[1])/2;
  *csSDn=3*SCcoeff*(nA5[0]*nA5[0]+nA5[1]*nA5[1])/2;
  printf("%s[%s]-nucleon cross sections[pb]: proton SI %.3E SD %.3E, neutron SD %.3E\n",
         CDM[k],antiParticle(CDM[k]),*csSIp,*csSDp,*csSDn);
}

int main(int argc,char** argv)
{  int err;
   char cdmName[10];
   double Omega=0;

  ForceUG=1;  /* the model is written in the unitary gauge */
  VZdecay=0; VWdecay=0;  cleanDecayTable();

  if(argc==1)
  {
      printf(" Correct usage:  ./main  <file with parameters> \n");
      printf("Example: ./main data1.par\n");
      exit(1);
  }

  err=readVar(argv[1]);
  if(err==-1)     {printf("Can not open the file\n"); exit(1);}
  else if(err>0)  { printf("Wrong file contents at line %d\n",err);exit(1);}

  err=sortOddParticles(cdmName);
  if(err) { printf("Can't calculate %s\n",cdmName); return 1;}

  printf("\nDM candidates: %s (M=%.3E GeV), %s (M=%.3E GeV)\n",
         CDM[1],McdmN[1],CDM[2],McdmN[2]);
  printf("MZ=%.6E MZp=%.4E MX=%.4E vq=%.4E sin(chi)=%.3E epsDP=%.3E\n",
         findValW("MZ"),findValW("MZp"),findValW("MX"),findValW("vq"),
         findValW("schi"),findValW("epsDP"));
  printf("sin(alpha)=%.3E lamHX=%.3E\n",findValW("sa"),findValW("lamHX"));

#ifdef MASSES_INFO
  printf("\n=== MASSES OF HIGGS AND ODD PARTICLES: ===\n");
  printHiggs(stdout);
  printMasses(stdout,1);
#endif

#ifdef OMEGA
{ int fast=1;
  double Beps=1.E-5;

  printf("\n==== Calculation of relic density =====\n");
  Tend=1E-3;
  Omega=darkOmega2(fast,Beps,&err);
  if(err<16)
  { printf("Omega h^2 total=%.3E  %s: %.3E  %s: %.3E\n",
           Omega,CDM[1],Omega*fracCDM[1],CDM[2],Omega*fracCDM[2]);
  } else printf("darkOmega2 failed, err=%d\n",err);
}
#endif

#ifdef CDM_NUCLEON
{ double siD,sdpD,sdnD, siM,sdpM,sdnM;
  printf("\n==== CDM-nucleon cross sections =====\n");
  nucleonCS(1,&siD,&sdpD,&sdnD);
  nucleonCS(2,&siM,&sdpM,&sdnM);
  /* Eq. of Fig. 4 of the paper: sigma_SI,eff = sigma_SI,p * Omega_Dirac h^2/0.12 */
  if(Omega>0) printf("sigma_SI,eff(%s) = %.3E pb\n",CDM[1],siD*Omega*fracCDM[1]/0.12);
}
#endif

#ifdef CDM_NUCLEUS
{ char* expName;
  printf("\n===== Direct detection exclusion:======\n");
  double pval=DD_pval(AllDDexp, Maxwell, &expName);
  if(pval<0.1 )  printf("Excluded by %s  %.1f%%\n", expName, 100*(1-pval));
  else printf("Not excluded by DD experiments  at 90%% level \n");
}
#endif

#ifdef DECAYS
{ txtList L;
  char *names[3]={"Zp","hX","h"};
  printf("\n================= Decays ==============\n");
  for(int i=0;i<3;i++)
  { double width=pWidth(names[i],&L);
    printf("\n%s :   total width=%.3E \n",names[i],width);
    printTxtList(L,stdout);
  }
}
#endif

#ifdef CLEAN
  system("rm -f HB.* HS.* hb.* hs.*  debug_channels.txt debug_predratio.txt  Key.dat");
#endif

  killPlots();
  return 0;
}

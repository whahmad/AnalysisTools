#include "TauSolver.h"
#include "PDG_Var.h"

TauSolver::TauSolver(TVector3 Tau_, TLorentzVector a1_):
  Tau(Tau_),
  a1(a1_)
{
  if(Tau.Mag()!=0)Tau*=1/Tau.Mag();
}

TauSolver::~TauSolver(){

}

void TauSolver::SolvebyRotation(TLorentzVector &Tau1,TLorentzVector &Tau2,TLorentzVector &nu1,TLorentzVector &nu2, int mode){
  TLorentzVector a1rot=a1;
  double phi(Tau.Phi()),theta(Tau.Theta());
  a1rot.RotateZ(-phi);
  a1rot.RotateY(-theta);
  double Enu1(0), Enu2(0), Ea1(a1.E()),Pz(a1rot.Pz()),mtau(PDG_Var::Tau_mass()),Pt(a1rot.Pt()),ma1(a1.M());
  if(mode==PZ)PzSolver(Enu1,Enu2,Ea1,ma1,Pz,Pt);
  if(mode=E)ESolver(Enu1,Enu2,Ea1,ma1,Pz,Pt);
  TLorentzVector Neutrino1(-a1rot.Px(),-a1rot.Py(),sqrt(Enu1*Enu1-Pt*Pt),Enu1);
  Neutrino1.RotateY(theta);
  Neutrino1.RotateZ(phi);
  Tau1=a1+Neutrino1;
  nu1=Neutrino1;
  TLorentzVector Neutrino2(-a1rot.Px(),-a1rot.Py(),sqrt(Enu2*Enu2-Pt*Pt),Enu2);
  Neutrino2.RotateY(theta);
  Neutrino2.RotateZ(phi);
  Tau2=a1+Neutrino2;
  nu2=Neutrino2;
}

void TauSolver::SolvebyProjection(TLorentzVector &Tau1,TLorentzVector &Tau2,TLorentzVector &nu1,TLorentzVector &nu2, int mode){
  double Enu1(0), Enu2(0), Ea1(a1.E()),Pz(0),mtau(PDG_Var::Tau_mass()),Pt(0),ma1(a1.M());
  Pz=Tau.Dot(a1.Vect());
  TVector3 Pz_vec=Tau;
  Pz_vec*=Pz;
  TVector3 Pt_vec=a1.Vect();
  Pt_vec-=Pz_vec;
  Pt=Pt_vec.Mag();
  if(mode==PZ)PzSolver(Enu1,Enu2,Ea1,ma1,Pz,Pt);
  if(mode=E)ESolver(Enu1,Enu2,Ea1,ma1,Pz,Pt);
  double Pnuz1=0;if(Enu1-Pt>0)Pnuz1=sqrt(Enu1*Enu1-Pt*Pt);
  double Pnuz2=0;if(Enu2-Pt>0)Pnuz1=sqrt(Enu2*Enu2-Pt*Pt);
  TLorentzVector Neutrino1(Pnuz1*Tau.Px()-Pt_vec.Px(),Pnuz1*Tau.Py()-Pt_vec.Py(),Pnuz1*Tau.Pz()-Pt_vec.Pz(),Enu1);
  Tau1=a1+Neutrino1;
  nu1=Neutrino1;
  TLorentzVector Neutrino2(Pnuz2*Tau.Px()-Pt_vec.Px(),Pnuz2*Tau.Py()-Pt_vec.Py(),Pnuz2*Tau.Pz()-Pt_vec.Pz(),Enu2);
  Tau1=a1+Neutrino2;
  nu2=Neutrino1;
}

void TauSolver::quadratic(double &x_plus,double &x_minus,double a, double b, double c){
  double R=b*b-4*a*c;
  if(R<0){R=0;}
  x_minus=(-b-sqrt(R))/(2.0*a);
  x_plus=(-b+sqrt(R))/(2.0*a);
}

void TauSolver::ESolver(double &Enu1,double &Enu2,double Ea1,double ma1, double Pz, double Pt){
  double mtau(PDG_Var::Tau_mass());
  double a=1-Ea1*Ea1/(Pz*Pz);
  double K=(mtau*mtau-ma1*ma1-2*Pt*Pt);
  double b=K*Ea1/(Pz*Pz);
  double c=-(Pt*Pt+K*K/(4*Pz*Pz));
  quadratic(Enu1,Enu2,a,b,c);
}

void TauSolver::PzSolver(double &Enu1,double &Enu2,double Ea1,double ma1, double Pz, double Pt){
  double mtau(PDG_Var::Tau_mass());
  double a=1-Pz*Pz/(Ea1*Ea1);
  double K=(mtau*mtau-ma1*ma1-2.0*Pt*Pt);
  double b=-K*Pz/(Ea1*Ea1);
  double c=Pt*Pt-K*K/(4.0*Ea1*Ea1);
  double Pnuz1(0),Pnuz2(0);
  quadratic(Pnuz1,Pnuz2,a,b,c);
  Enu1=sqrt(Pnuz1*Pnuz1+Pt*Pt);
  Enu2=sqrt(Pnuz2*Pnuz2+Pt*Pt);
  if(Pnuz1<0)Enu1=Pt;
  if(Pnuz2<0)Enu2=Pt;
}

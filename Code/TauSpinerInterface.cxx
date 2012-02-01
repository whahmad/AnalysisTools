#include "TauSpinerInterface.h"

//#include "HepMC/IO_GenEvent.h"
#include "Tauola.h"
#include "LHAPDF/LHAPDF.h"
#include "tau_reweight_lib.h"
#include "read_particles_from_TAUOLA.h"
#include<iostream>
#include "TauDataFormat/TauNtuple/interface/PdtPdgMini.h"

TauSpinerInterface::TauSpinerInterface(){
  ////////////////////////////////////////////////////////////
  //
  // This class is based on TauSpiner Example
  Tauola::initialize();
  string name="MSTW2008nnlo90cl.LHgrid";
  LHAPDF::initPDFSetByName(name);
  double CMSENE = 7000.0; // center of mass system energy.
                          // used in PDF calculation. For pp collisions only
  bool Ipp = true;  // for pp collisions
  // Initialize TauSpinner
  initialize_spiner(Ipp, CMSENE);
}

TauSpinerInterface::~TauSpinerInterface(){

}

double TauSpinerInterface::Get(TauSpinerType type, SimpleParticle X, SimpleParticle tau, std::vector<SimpleParticle> tau_daughters,SimpleParticle tau2, std::vector<SimpleParticle> tau_daughters2){
  if(hminus==type || hplus==type){
    double HHState=tautauHelicityState(X, tau, tau2, tau_daughters,tau_daughters2);
    if(HHState>0.0 && hplus==type)  return 1.0;
    if(HHState<0.0 && hminus==type) return 1.0;
    return 0.0;
  }
  if(LPolarization==type) {
    double S=X.e()*X.e() - X.px()*X.px() - X.py()*X.py() - X.pz()*X.pz();
    if(tau2.pdgid()==PdtPdgMini::tau_minus)return getLongitudinalPolarization(S, tau, tau2);
    return getLongitudinalPolarization(S, tau2, tau);
  }
  double WT    = 1.0; // assume that there is 1 bosons decaying into                                                                                                                                                                        
  // Calculate weight for first boson
  if( abs(X.pdgid())==24 ||  abs(X.pdgid())==37 ){
    WT = calculateWeightFromParticlesWorHpn(X, tau, tau2, tau_daughters); // note that tau2 is tau neutrino 
  }
  else if( X.pdgid()==25 || X.pdgid()==36 || X.pdgid()==22 || X.pdgid()==23 ){
    WT = calculateWeightFromParticlesH(X, tau, tau2, tau_daughters,tau_daughters2);
  }
  else{
    cout<<"WARNING: Unexpected PDG for tau mother: "<<X.pdgid()<<endl;
  }
  if(Spin==type) return WT;
  if(UnSpin==type) return 1.0/WT;
  if(FlipSpin==type) return (2.0-WT)/(WT);
  std::cout << "Warning TauSpinerWeight TauSpinerType " << type << " is INVALID. Returning Spin WT=1.0." << std::endl;
  return 1.0;
}


double TauSpinerInterface::tautauHelicityState(SimpleParticle &sp_X, SimpleParticle &sp_tau1, SimpleParticle &sp_tau2, std::vector<SimpleParticle> &sp_tau1_daughters, std::vector<SimpleParticle> &sp_tau2_daughters)
{
  SimpleParticle         sp_tau;
  SimpleParticle         sp_nu_tau;
  std::vector<SimpleParticle> sp_tau_daughters;
  
  // First iteration is for tau plus, so the 'nu_tau' is tau minus
  if (sp_tau1.pdgid() == -15 )
  {
    sp_tau           = sp_tau1;
    sp_nu_tau        = sp_tau2;
    sp_tau_daughters = sp_tau1_daughters;
  }
  else
  {
    sp_tau           = sp_tau2;
    sp_nu_tau        = sp_tau1;
    sp_tau_daughters = sp_tau2_daughters;
  }

  double *HHp, *HHm;
  
  // We use this to separate namespace for tau+ and tau-
  if(true)
  {
    // Create Particles from SimpleParticles
    Particle X     (      sp_X.px(),      sp_X.py(),      sp_X.pz(),      sp_X.e(),      sp_X.pdgid() );
    Particle tau   (    sp_tau.px(),    sp_tau.py(),    sp_tau.pz(),    sp_tau.e(),    sp_tau.pdgid() );
    Particle nu_tau( sp_nu_tau.px(), sp_nu_tau.py(), sp_nu_tau.pz(), sp_nu_tau.e(), sp_nu_tau.pdgid() );

    vector<Particle> tau_daughters;

    // tau pdgid
    int tau_pdgid = sp_tau.pdgid();

    // Create list of tau daughters
    for(int i=0; i<sp_tau_daughters.size(); i++)
    {
      Particle pp(sp_tau_daughters[i].px(),
                  sp_tau_daughters[i].py(),
                  sp_tau_daughters[i].pz(),
                  sp_tau_daughters[i].e(),
                  sp_tau_daughters[i].pdgid() );

      tau_daughters.push_back(pp);
    }

    double phi2 = 0.0, theta2 = 0.0;


    //  Move decay kinematics first to tau rest frame  with z axis pointing along nu_tau direction
    //  later rotate again to have neutrino from tau decay along z axis: angles phi2, theta2
    prepareKinematicForHH   (tau, nu_tau, tau_daughters, &phi2, &theta2);


    //  Determine decay channel and then calculate polarimetric vector HH
    HHp = calculateHH(tau_pdgid, tau_daughters, phi2, theta2);

  } // end of tau+

  // Second iteration is for tau minus, so the 'nu_tau' is tau minus
  if(sp_tau1.pdgid() == 15 )
  {
    sp_tau           = sp_tau1;
    sp_nu_tau        = sp_tau2;
    sp_tau_daughters = sp_tau1_daughters;
  }
  else
  {
    sp_tau           = sp_tau2;
    sp_nu_tau        = sp_tau1;
    sp_tau_daughters = sp_tau2_daughters;
  }
  
  // We use this to separate namespace for tau+ and tau-
  if(true)
  {
    // Create Particles from SimpleParticles
    Particle X     (      sp_X.px(),      sp_X.py(),      sp_X.pz(),      sp_X.e(),      sp_X.pdgid() );
    Particle tau   (    sp_tau.px(),    sp_tau.py(),    sp_tau.pz(),    sp_tau.e(),    sp_tau.pdgid() );
    Particle nu_tau( sp_nu_tau.px(), sp_nu_tau.py(), sp_nu_tau.pz(), sp_nu_tau.e(), sp_nu_tau.pdgid() );

    vector<Particle> tau_daughters;

    // tau pdgid
    int tau_pdgid = sp_tau.pdgid();

    // Create list of tau daughters
    for(int i=0; i<sp_tau_daughters.size(); i++)
    {
      Particle pp(sp_tau_daughters[i].px(),
                  sp_tau_daughters[i].py(),
                  sp_tau_daughters[i].pz(),
                  sp_tau_daughters[i].e(),
                  sp_tau_daughters[i].pdgid() );

      tau_daughters.push_back(pp);
    }

    double phi2 = 0.0, theta2 = 0.0;


    //  Move decay kinematics first to tau rest frame  with z axis pointing along nu_tau direction
    //  later rotate again to have neutrino from tau decay along z axis: angles phi2, theta2
    prepareKinematicForHH   (tau, nu_tau, tau_daughters, &phi2, &theta2);


    //  Determine decay channel and then calculate polarimetric vector HH
    HHm = calculateHH(tau_pdgid, tau_daughters, phi2, theta2);

  } // end of tau-

  // polarization vector (unit vector) has been rotated such that the z axis is allon the tau direction 
  // (note: polarization is -s where s is the tau spin direction) 
  double Helicity_taup=0.5*HHp[2]/fabs(HHp[2]);
  double Helicity_taum=0.5*HHm[2]/fabs(HHm[2]);

  delete[] HHp;
  delete[] HHm;
  
  return Helicity_taup+Helicity_taum;
}


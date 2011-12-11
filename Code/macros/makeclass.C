void makeclass(){
   TChain c("t");
   c.Add("/user/inugent/MC_DY_SkimmedTauNtuple_1_1_1Yh.root");
   c.MakeClass("NtupleReader");
}

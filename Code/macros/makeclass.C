void makeclass(){
   TChain c("t");
   c.Add("/net/scratch_cms/institut_3b/cherepanov/MC_DY_SkimmedTauNtuple.root");
   c.MakeClass("NtupleReader");
}

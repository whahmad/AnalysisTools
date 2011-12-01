void makeclass(){
   TChain c("t");
   c.Add("~/Software/TauNtupleMaker/CMSSW_4_2_4/src/TauDataFormat/TauNtuple/output.root");
   c.MakeClass("NtupleReader");
}

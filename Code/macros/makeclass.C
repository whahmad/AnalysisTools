void makeclass(){
   TChain c("t");
   c.Add("../Set_1/TauNtuple_123_1_t4m.root");
   c.MakeClass("NtupleReader");
}

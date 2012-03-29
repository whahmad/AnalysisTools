void makeclass(){
   TChain c("t");
   c.Add("/user/scratch/nugent/TauNtuple_115_1_4P9.root");
   c.MakeClass("NtupleReader");
}

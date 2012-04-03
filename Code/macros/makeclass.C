void makeclass(){
   TChain c("t");
   c.Add("/net/scratch_cms/institut_3b/nugent/workdirV1/TauNtuple_27_1_I2n.root");
   c.MakeClass("NtupleReader");
}

void makeclass(TString file){
   TChain c("t");
   c.Add(file);
   c.MakeClass("NtupleReader");
}

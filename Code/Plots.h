#ifndef Plots_h
#define Plots_h

#include <string.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include<vector>


class Plots {
 public:
 Plots();
 ~Plots();

  void Set_Data_Type(TString Data_type){
    Data_type_=Data_type;
  }

  void SaveHistograms(TString File, std::vector<TString> HistogramNames);


 void PlotStyle();
 //void SetAtlasStyle ();
 void AtlasStyle();
 void ATLASLabel(Double_t x,Double_t y,bool Preliminary,Color_t color);

 void Plot1D(std::vector<TH1D> histo,float Lumi,std::vector<float> CrossSectionandAcceptance,std::vector<float> nevents,std::vector<int> colour,std::vector<TString> legend);
 void Plot1D(std::vector<std::vector<TH1D> > histo,float Lumi,std::vector<float> CrossSectionandAcceptance,std::vector<float> nevents,std::vector<int> colour,std::vector<TString> legend);
 void Plot2D(std::vector<TH2D>  histo,float Lumi,std::vector<float> CrossSectionandAcceptance,std::vector<float> nevents,std::vector<int> colour,std::vector<TString> legend);
 void Make_Figure(TString name,TString cap);
 void Plot1DSignificance(std::vector<TH1D> histo, bool gt,bool lt,float Lumi,std::vector<float> CrossSectionandAcceptance,std::vector<float> nevents,std::vector<int> colour,std::vector<TString> legend);
 void Plot1Dsigtobkg(std::vector<TH1D> histo, bool gt,bool lt,float Lumi,std::vector<float> CrossSectionandAcceptance,std::vector<float> nevents,std::vector<int> colour,std::vector<TString> legend);
 void Plot1D_DataMC_Compare(std::vector<TH1D> histo, float Lumi,std::vector<float> CrossSectionandAcceptance,std::vector<float> nevents,std::vector<int> colour,std::vector<TString> legend);

 private:
 static TString File_;
 static std::vector<TString> HistogramNames_;
 static TString Data_type_;
 bool doscale;
 bool doprofiles;
 bool verbose;
 bool dooneprofile;
};
#endif

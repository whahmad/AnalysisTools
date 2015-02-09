#ifndef SVFitObject_H
#define SVFitObject_H

class SVFitObject{
 public:
  SVFitObject(){};
  virtual ~SVFitObject(){};
  
  void Set_a(double a_){a=a_;}
  double Get_a(){return a;}

private:
  double a;

};

#endif 

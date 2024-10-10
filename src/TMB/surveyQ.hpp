template<class Type>
struct surveyQ_t{
public:

  
  Type get_Q(int length,int year){return 0;};

};


//Logistic shaped Q
template<class Type>
struct logisticQ_t : surveyQ_t<Type>{

  objective_function<Type> *objective;

  std::string prefix;


  //parameters
  Type Qmax;
  Type QL50;
  Type QL95;
  Type QK;

  vector<Type> QL;
  int L;
  Type bin_adjust;

  //Q prior related things
  Type Q_max_sd;
  vector<Type> Qpriors;

  //Constructors
  logisticQ_t(objective_function<Type> *obj,std::string pre,int L,Type bin_adjust,int start_length);

  Type get_Q(int length,int year);
  void Q_prior();
  void report();

};

template<class Type>
logisticQ_t<Type>::logisticQ_t(objective_function<Type> *obj, std::string pre, int L, Type bin_adjust,int start_length){
  objective = obj;
  prefix = pre;

  Qmax = Parameter("log_Qmax",obj,pre);
  Qmax = exp(Qmax);

  QL50 = Parameter("log_QL50",obj,pre);
  QL50 = exp(QL50);

  QL95 = Parameter("log_QL95",obj,pre);
  QL95 = exp(QL95);

  QK = -1*log(0.05/0.95)/(QL95-QL50);
  vector<Type> QLtemp(L);
  QL = QLtemp;


  for(int l = 0; l < L; ++l){
    Type len = Type(l)+bin_adjust+Type(start_length);
    QL(l) = Qmax/(1+exp(-QK*(len-QL50)));
  }

}

template<class Type> 
Type logisticQ_t<Type>::get_Q(int length,int year){
  return QL(length);
}


template<class Type>
void logisticQ_t<Type>::report(){
  Report("Qmax",Qmax,objective,prefix);
  Report("QL50",QL50,objective,prefix);
  Report("QL95",QL95,objective,prefix);
  Report("QK",QK,objective,prefix);

  Adreport("Qmax",Qmax,objective,prefix);
  Adreport("QL50",QL50,objective,prefix);
  Adreport("QL95",QL95,objective,prefix);
  Adreport("QK",QK,objective,prefix);
}


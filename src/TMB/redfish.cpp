#define TMB_LIB_INIT R_init_redfish
#include <TMB.hpp>
#include "data_parameter.hpp"
#include "transformations.hpp"
#include "sigmaGen.hpp"
#include "schnutevb.hpp"
#include "projections.hpp"
#include "help.hpp"
#include "surveyQ.hpp"
#include "pnorm4.hpp"

//' @name isNA
//' @title Check if Type is NA
//' @description Check if a Type variable is NA by R standards
//' @param x the variable to check
//' 
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

using namespace density;



//' @name fit_std_landings
//' @title Fit Landings
//' @description Fit reported landings assuming uncensored normally distributed data
//' @param nll the negative log likelihood being tracked by the model
//' @param log_expected_landings the expected landings predicted by the model at the given parameters
//' @param obj pointer to the objective function
template<class Type>
void fit_std_landings(Type &nll,vector<Type> &log_expected_landings, objective_function<Type> *obj){
  vector<Type> landing_nums = DataVector("landing_nums",obj);
  Type landings_sd = Parameter("log_landings_sd",obj);
  landings_sd = exp(landings_sd);

  vector<Type> log_landing_resids(landing_nums.size());
  vector<Type> std_landing_resids(landing_nums.size());

  if(check_simulate(obj)){
    for(int i = 0; i < landing_nums.size(); ++i){
      landing_nums(i) = exp(rnorm(log_expected_landings(i),landings_sd));
      Report("landing_nums",landing_nums,obj);
    }
  }
	
  for(int i = 0; i < landing_nums.size(); ++i){
    log_landing_resids(i) = log(landing_nums(i)) - log_expected_landings(i);
    std_landing_resids(i) = log_landing_resids(i)/landings_sd;
    nll -= dnorm(log(landing_nums(i)),log_expected_landings(i),landings_sd,true);
  }

  Adreport("landings_sd",landings_sd,obj);
  Adreport("log_expected_landings",log_expected_landings,obj);
  Report("log_landing_resids",log_landing_resids,obj);
  Report("std_landing_resids",std_landing_resids,obj);
}

//' @name fit_censored_landings
//' @title Fit censored landings
//' @description Fit reported landings assuming they are interval censored normally
//' @param nll the negative log likelihood being tracked by the model
//' @param log_expected_landings the expected landings predicted by the model at the given parameters
//' @param obj pointer to the objective function
template<class Type>
void fit_censored_landings(Type &nll,vector<Type> &log_expected_landings, objective_function<Type> *obj){
  vector<Type> landing_nums = DataVector("landing_nums",obj);
  Type landings_sd = Parameter("log_landings_sd",obj);
  landings_sd = exp(landings_sd);

  Rcout << "Using Cens. Bounds" << std::endl;

  vector<Type> log_landingsL = DataVector("log_landingsL",obj);
  vector<Type> log_landingsU = DataVector("log_landingsU",obj);
  vector<Type> log_landing_resids(landing_nums.size());
  vector<Type> std_landing_resids(landing_nums.size());


  if(check_simulate(obj)){
    for(int i = 0; i < landing_nums.size(); ++i){
      //change this
      landing_nums(i) = exp(rnorm(log_expected_landings(i),landings_sd));
      
    }
  }


  vector<Type> landnll(landing_nums.size());
  for(int i = 0; i < landing_nums.size(); ++i){
    log_landing_resids(i) = log(landing_nums(i)) - log_expected_landings(i);
    std_landing_resids(i) = log_landing_resids(i)/landings_sd;
      landnll(i) = censored_bounds(log(landing_nums(i)),log_expected_landings(i),landings_sd,-log_landingsL(i),log_landingsU(i));
      nll -= censored_bounds(log(landing_nums(i)),log_expected_landings(i),landings_sd,-log_landingsL(i),log_landingsU(i));
    
  }

  Adreport("landings_sd",landings_sd,obj);
  Adreport("log_expected_landings",log_expected_landings,obj);
  Report("log_landing_resids",log_landing_resids,obj);
  Report("std_landing_resids",std_landing_resids,obj);
  Report("landnll",landnll,obj);

}


//' @title Convert proportions to Continuation Ratio Logit
//' @param proportions The vector to convert to CRL
template<class Type>
vector<Type> make_CRL_vec(vector<Type> proportions){
  int S = proportions.size();
  vector<Type> crl(S-1);

  for(int s= 0; s < S-1; ++s){
    // Type denom = 0;
    // for(int j=s; j < S; ++j){
    //   denom += proportions(j);
    // }
    Type denom = proportions.segment(s,S-s).sum();
    Type num = proportions(s);
    Type pi = num/denom;
    crl(s) = log(pi/(1-pi));
  }

  return crl;
  
}

//' @title Convert vector to vector of proportions
//' @param something the vector to convert to proportions
template<class Type>
vector<Type> make_proportions(vector<Type> something){
  Type summy = something.sum();
  vector<Type> props(something.size());
  
  for(int i = 0; i < something.size(); ++i){
    props(i) = something(i)/summy;
  }
  return props;
}

//' @title Aggregate model catch to reflect observed sizes
//' @description Given a predetermined key, convert model catch to match up with observed freq.
//' @param catch_vec the catch vector to aggregate
//' @param catch_key the catch key for a given year
//' @param n_size the size of the aggregated catch vector
template<class Type>
vector<Type> aggregate_catch(vector<Type> catch_vec,vector<int> catch_key,int n_size){
  vector<Type> aggr_catch(n_size);

  for(int i = 0; i < n_size; ++i){
    aggr_catch(i) = 0;
  }

  for(int i = 0; i < catch_vec.size(); ++i){
    aggr_catch(catch_key(i)) += catch_vec(i);
  }

  return aggr_catch;
}

//' @name fit_survey
//' @description Fit observed survey indices, these are read in from R as a list per year
//' @param nll the negative log likelihood being tracked by the model
//' @param QLM the matrix containing the catchabilities by length
//' @param N the numbers at age matrix to use
//' @param L3 the value of L3
//' @param Y the total number of years in the model
//' @param start_length the start length of the model
//' @param end_length the end length of the model
//' @param obj pointer to the objective function
template<class Type>
void fit_survey(Type &nll,matrix<Type> &QLM,matrix<Type> &N,int L3,int Y, int start_length,int end_length,objective_function<Type>* obj){
  using namespace density;
  
  //Survey indices and lenghts are kept in a list because I really hate dealing with indexing...
  SEXP survey_list = getListElement(obj -> data,"survey_list");
  //number of years to process
  int n_years = LENGTH(survey_list);

  int survey_size = DataInteger("survey_size",obj);
  int survey_sigtype = DataInteger("survey_sigtype",obj);

  vector<Type> rhoS(1);

  switch(survey_sigtype){
  case 0:
    rhoS(0) = 0;
    break;
  case 1:
    rhoS(0) = Parameter("logit_rhoS",obj,"");
    break;
  default:
    error("Survey correlation mode not supported");
    break;
  }

  vector<Type> sd_survey = ParameterVector("log_sd_survey",obj);

  //vector<Type> log_exp_index;
  //vector<Type> exp_index;
  vector<Type> log_exp_index(survey_size);
  vector<Type> exp_index(survey_size);
  vector<Type> log_survey_resids(survey_size);
  vector<Type> std_log_survey_resids(survey_size);
  //vector<Type> total_survey_abundance;
  vector<Type> total_survey_abundance(n_years);

  vector< matrix<Type> > Ssigmas(n_years);

  int r_proj = DataInteger("r_proj",obj);

  
  for(int yy = 0; yy < n_years; ++yy){
    // get the current list
    SEXP curr_list = VECTOR_ELT(survey_list,yy);

    Type tyear = asVector<Type>(getListElement(curr_list,"year",&isNumericScalar))[0];
    int year = CppAD::Integer(tyear);

    Type ttype = asVector<Type>(getListElement(curr_list,"type",&isNumericScalar))[0];
    int stype = CppAD::Integer(ttype);
    
    Type tprojy = asVector<Type>(getListElement(curr_list,"projy",&isNumericScalar))[0];
    int projy = CppAD::Integer(tprojy);

    
    vector<Type> s_ind = asVector<Type>(getListElement(curr_list,"indices",&isNumeric));
    vector<int> s_len = asVector<int>(getListElement(curr_list,"lengths",&isNumeric));
    vector<int> s_place = asVector<int>(getListElement(curr_list,"place",&isNumeric));
    vector<int> s_map = asVector<int>(getListElement(curr_list,"mmap",&isNumeric));


    Type tnlen = asVector<Type>(getListElement(curr_list,"nlen",&isNumericScalar))[0];
    int nlen = CppAD::Integer(tnlen);
    
    
    vector<Type> exp_ind(s_ind.size());
    vector<Type> c_sd(s_ind.size());

    for(int i = 0; i < s_ind.size(); ++i){
      Type Q = QLM(s_len(i),stype);
      if(s_len(i) < start_length){
	Type NLminus = 0;
	for(int l = 0; l < start_length; ++l){
	  Type Qt = QLM(l,stype);
	  NLminus += Qt*N(l,year);
	}
	exp_ind(i) = NLminus;
      }else if(s_len(i) != end_length-1){
	exp_ind(i) = Q*N(s_len(i),year);
      }else{
	Type NLplus = 0;
	for(int l = end_length-1; l < L3; ++l){
	  Type QLt = QLM(l,stype); 
	  NLplus += QLt*N(l,year);
	}
	exp_ind(i) = NLplus;
      }

      exp_index(s_place(i)) = exp_ind(i);
      log_exp_index(s_place(i)) = log(exp_ind(i));
      if(projy == 0){
	c_sd(i) = sd_survey(s_map(i));
      } 
    }
    total_survey_abundance(yy) = sum(exp_ind);

    // vector<Type> c_sd = sd_survey.col(yy);
      
    
      vector<Type> log_diff = log(s_ind)-log(exp_ind);

      matrix<Type> Ssigma;
      switch(survey_sigtype){
      case 0:
	//Rcout << "c_sd:" << c_sd << std::endl;
	Ssigma = sigmaGen(nlen,c_sd);
	break;
      case 1:
	Ssigma = sigmaGen(nlen,3,c_sd,rhoS,0);
	break;
      default:
	error("survey sigma not supported");
      }

      nll += MVNORM(Ssigma)(log_diff);
    
      Ssigmas(yy) = Ssigma;

    
      Eigen::LLT< Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> > lltDiff(Ssigma);
      matrix<Type> LD = lltDiff.matrixL();
      matrix<Type> LinvD = LD.inverse();
      vector<Type> std_resids = LinvD*(log_diff);

    
      for(int i = 0; i < s_ind.size(); ++i){
	log_survey_resids(s_place(i)) = log_diff(i);
      
	std_log_survey_resids(s_place(i)) = std_resids(i);
      }
    
  }

  if(r_proj == 1){
  vector<Type> ttotal_survey_abundance(Y);
    for(int y = 0; y < n_years; ++y){
       ttotal_survey_abundance(y) = total_survey_abundance(y);
     }
    for(int y = n_years; y < Y; ++y){
       ttotal_survey_abundance(y) = 0;
     }
    
     int diffy = Y-n_years;
     int addspace = QLM.rows()*diffy;
     vector<Type> texp_index(exp_index.size()+addspace);
     for(int j = 0; j < exp_index.size(); ++j){
       texp_index(j) = exp_index(j);
     }	
     int qq = exp_index.size();
     for(int jj = n_years; jj < Y; ++jj){
        for(int zz = 0; zz < QLM.rows(); ++zz){
    // 	 //Rcout << "qq: " << qq << std::endl;
      	texp_index(qq) = QLM(zz,1)*N(zz,jj);
	ttotal_survey_abundance(jj) += texp_index(qq);
      	qq++;
        }
    //  
      }
     //Report("total_survey_abundance",ttotal_survey_abundance,obj);
    Report("exp_index",texp_index,obj);
    vector<Type> log_total_survey_abundance = ttotal_survey_abundance.log();
    Adreport("log_total_survey_abundance",log_total_survey_abundance,obj,"",true);
    
    Adreport("total_survey_abundance",ttotal_survey_abundance,obj,"",true);
    vector<Type> log_exp_index = log(exp_index);
    Adreport("log_exp_index",log_exp_index,obj);
    //total_survey_abundance = ttotal_survey_abundance;
    //exp_index = texp_index;
    // log_exp_index = log(texp_index);
   }
  

  if(r_proj == 0){
    vector<Type> log_total_survey_abundance = total_survey_abundance.log();
    Adreport("log_total_survey_abundance",log_total_survey_abundance,obj,"",true);
  
    Report("log_exp_index",log_exp_index,obj);
    Adreport("log_exp_index",log_exp_index,obj,"",true);
    Report("exp_index",exp_index,obj);
  }
  Report("log_survey_resids",log_survey_resids,obj);
  Report("std_log_survey_resids",std_log_survey_resids,obj);

  Report("Ssigmas",Ssigmas,obj);
  Adreport("sd_survey",exp(sd_survey),obj);
  if(survey_sigtype == 1){
    Type ross = rho_trans(rhoS(0));
    Adreport("rhoS",ross,obj);
  }
  
}

//' @name fit_catch_prop
//' @description Fit observed catch proportions using CRLs, these are read in from R as a list per year
//' @param nll the negative log likelihood being tracked by the model
//' @param catch_at_length the catch at length estimated by the model
//' @param obj pointer to the objective function
template<class Type>
void fit_catch_prop(Type &nll, matrix<Type> &catch_at_length, objective_function<Type>* obj){

  using namespace density;
  
  //We keep the catch proportion keys and stuff in a list
  //std::string cl_nam = "catch_list";
  //char* cl_cnam = new char[cl_nam.size()+1];
  SEXP catch_list = getListElement(obj -> data,"catch_list");
  //Number of years to process
  int n_years = LENGTH(catch_list);

  //Get parameters
  //Type log_sd_catch_prop = Parameter("log_sd_catch_prop",obj,"");
  //matrix<Type> log_sd_catch_prop_m = ParameterMatrix("log_sd_catch_prop_m",obj);
  vector<Type> log_sd_catch_prop_m = ParameterVector("log_sd_catch_prop_m",obj);

  
  vector< vector<Type> > obs_crls(n_years);
  vector< vector<Type> > exp_props(n_years);
  vector< vector<Type> > diffs(n_years);
  vector< matrix<Type> > sigmas(n_years);
  vector< vector<Type> > agg_catches(n_years);
  vector<int> years(n_years);
  vector<vector<int> > keys(n_years);
  vector<vector<Type> > props(n_years);
  vector<vector<Type> > agg_props(n_years);
  vector<vector<Type> > std_diffs(n_years);
  vector<vector<Type> > logit_exp_props(n_years);

  enum catch_sigma_t{
    independent_c = 0,
    ARone_c = 1,
    compound_c
  };

  int catch_sig_sw = DataInteger("catch_sig_sw",obj,"");

  vector<Type> rhoC(1);

  switch(catch_sig_sw){
  case independent_c:
    rhoC(0) = 0;
    break;
  case ARone_c:
    rhoC = ParameterVector("logit_rhoC",obj,"");
    break;
  case compound_c:
    rhoC = ParameterVector("logit_rhoC",obj,"");
    break;
  default:
    error("Catch correlation mode not supported");
    break;
  }

  Report("n_years_catch_prop",n_years,obj);
  
   for(int i = 0; i < n_years; ++i){
     //get the CURRENT list
     SEXP curr_list = VECTOR_ELT(catch_list,i);
     
     //get the needed items from the list
     Type tyear = asVector<Type>(getListElement(curr_list,"year",&isNumericScalar))[0];
     int year = CppAD::Integer(tyear);

     Type tsize = asVector<Type>(getListElement(curr_list,"ysize",&isNumericScalar))[0];
     int ysize = CppAD::Integer(tsize);
     
     vector<int> key = asVector<int>(getListElement(curr_list,"key",&isNumeric));
        
     vector<Type> prop = asVector<Type>(getListElement(curr_list,"prop",&isNumeric));

     vector<int> c_map = asVector<int>(getListElement(curr_list,"mmap",&isNumeric));



     years(i) = year;
     keys(i) = key;
     props(i) = prop;
     
     vector<Type> curr_c_at_l = catch_at_length.col(year);
     vector<Type> agg_prop = aggregate_catch(prop,key,ysize);
     vector<Type> obs_crl = make_CRL_vec(agg_prop);

     vector<Type> agg_c_at_l = aggregate_catch(curr_c_at_l,key,ysize);
     vector<Type> exp_prop = make_proportions(agg_c_at_l);
     vector<Type> logit_exp_prop = logit(exp_prop);
     vector<Type> exp_crl = make_CRL_vec(exp_prop);

     vector<Type> diff(obs_crl.size());
     for(int j = 0; j < obs_crl.size(); ++j){
       diff(j) = obs_crl(j)-exp_crl(j);
     }

     int nlen = diff.size();

     // vector<Type> log_sd_catch_prop = log_sd_catch_prop_m.col(i);
     vector<Type> log_sd_catch_prop(nlen);
     for(int j = 0; j < nlen; ++j){
       log_sd_catch_prop(j) = log_sd_catch_prop_m(c_map(j));
     }
     
      //Create the covariance matrix
     matrix<Type> sigma;
     switch(catch_sig_sw){
     case independent_c:
       sigma = sigmaGen(nlen,log_sd_catch_prop);
       break;
     case ARone_c:
       sigma = sigmaGen(nlen,3,log_sd_catch_prop,rhoC,0);
       break;
     case compound_c:
       sigma = sigmaGen(nlen,2,log_sd_catch_prop,rhoC,0);
       break;
     default:
       error("catch sigma not supported");
     }
     
// matrix<Type> sigma = sigmaGen(nlen,log_sd_catch_prop);

      obs_crls(i) = obs_crl;
      exp_props(i) = exp_prop;
      logit_exp_props(i) = logit_exp_prop;
      diffs(i) = diff;
      sigmas(i) = sigma;
      agg_catches(i) = agg_c_at_l;
      agg_props(i) = agg_prop;
      Eigen::LLT< Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> > lltDiff(sigma);
      matrix<Type> LD = lltDiff.matrixL();
      matrix<Type> LinvD = LD.inverse();
      std_diffs(i) = LinvD*(diff);
  
      nll += MVNORM(sigma)(diff);
      
   }
   matrix<Type> sd_catch_prop_m = log_sd_catch_prop_m.array().exp();
   
   Report("years",years,obj);
   Report("props",props,obj);
   Report("agg_props",agg_props,obj);
   Report("keys",keys,obj);
   Report("obs_crls",obs_crls,obj);
   Report("exp_props",exp_props,obj);
   Report("diffs",diffs,obj);
   Report("std_diffs",std_diffs,obj);
   Report("sigmas",sigmas,obj);
   Report("agg_catches",agg_catches,obj);
   Adreport("sd_catch_prop_m",sd_catch_prop_m,obj);
   AdreportList("std_diffs",std_diffs,obj);
   AdreportList("logit_exp_props",logit_exp_props,obj);
   vector<Type> irho_C = invlogit(rhoC);
   Adreport("rhoC",irho_C,obj);
   
   
}


enum N_t {
  fixedN = 0,
  randN = 1
};

enum F_t{
  sep_no_disc = 0,
  AR = 3,
  rwrw = 4,
  simp_len_blocked = 5
};

enum GF_t{
  sig_schnute = 0,
  tv_schnute = 1,
  fixed_L = 2
};

enum M_t{
  fixedM = 0,
  Mdev = 1,
  rwM = 2,
  LorMAP = 3,
  LorMAPdev = 4
};

enum Q_t{
  logisticQ = 0,
  doubleQ = 1,
  spline_double = 2,
  spline_rho = 3,
  AR_rho = 4
};

enum weight_t{
  fixedW = 0,
  rawW = 1
};

enum landing_t{
  std_landing = 0,
  cens_bound = 1,
  discards = 2,
  prop_discards = 3
};

enum catch_t{
  no_prop = 0,
  prop = 1
};

enum sel_t{
  at_age = 0,
  at_len_logistic = 1
};

template<class Type>
Type objective_function<Type>::operator() ()
{

  using newton::newton_config_t;
  using newton::Newton;

    DATA_STRUCT(cfg, newton_config_t);
  // //Functor<TMBad::ad_aug> F(m, x);
  // stupid<TMBad::ad_aug> F2(X,y,nll);
  // vector<Type> sol = Newton(F2, parms, cfg);

  int proj_years = DataInteger("proj_years",this);


    
  //Data inputs
  //These use the functions defined in data_parameter.hpp instead of the
  //TMB Macros. They are more flexible for using in functions and some other scenarios and use consistent C++ syntax.
  int start_length = DataInteger("start_length",this);
  int end_length = DataInteger("end_length",this);
  int start_year = DataInteger("start_year",this);
  int end_year = DataInteger("end_year",this);
  int Y = DataInteger("Y",this);
  int A = DataInteger("A",this);
  int L = DataInteger("L",this);
  int L3 = DataInteger("L3",this);
  matrix<Type> weightsF = DataMatrix("weightsF",this);
  matrix<Type> maturityF = DataMatrix("maturityF",this);
  matrix<Type> maturityM = DataMatrix("maturityM",this);
  vector<Type> survey_index = DataVector("survey_index",this);
  vector<int> survey_year = DataIVector("survey_year",this);
  vector<int> survey_length = DataIVector("survey_length",this);
  vector<int> landing_year = DataIVector("landing_year",this);
  Type bin_adjust = DataScalar("bin_adjust",this);
  vector<int> survey_type = DataIVector("survey_type",this);
  int N_type = DataInteger("N_type",this);
  int F_type = DataInteger("F_type",this);
  int GF_type = DataInteger("GF_type",this);
  int M_type = DataInteger("M_type",this);
  int catch_type = DataInteger("catch_type",this);
  int Q_prior = DataInteger("Q_prior",this);
  int landing_type = DataInteger("landing_type",this);
  int adr_flag = DataInteger("adr_flag",this);
  vector<Type> Qpriors = DataVector("Qpriors",this);
  

  Type nll = 0.0;

  //We need a bunch of parameters early for things like projections...
  matrix<Type> log_N_a;
  vector<Type> log_recruit;

  Type recruit_sd = Parameter("log_recruit_sd",this);
  recruit_sd = exp(recruit_sd);

  if(N_type == 1){
    log_N_a = ParameterMatrix("log_N_a",this,"");
  }else{
    log_recruit = ParameterVector("log_recruit",this);
  }

  Type Fy_sd;
  vector<Type> log_Fy;
  Type Fa_sd;
  vector<Type> log_Fa;

   switch(F_type){
  case sep_no_disc:
    Fy_sd = Parameter("log_Fy_sd",this);
    Fy_sd = exp(Fy_sd);
    log_Fy = ParameterVector("log_Fy",this);
    log_Fa = ParameterVector("log_Fa",this);
    Fa_sd = Parameter("log_Fa_sd",this);
    exp(Fa_sd);
    break;
  case AR:
  case rwrw:
  case simp_len_blocked:
    Fy_sd = Parameter("log_Fy_sd",this);
    Fy_sd = exp(Fy_sd);
    log_Fy = ParameterVector("log_Fy",this);
    break;
  default:
    error("Incorrect F choice");
    break;
  };


  int start_age = 1;
  Type fall_adjust = 10.5/12;
  schnutevb_t<Type> growth_fun(this,"",A,GF_type,N_type,log_N_a);
 




  
  //F and Z
  matrix<Type> M(A,Y);
  matrix<Type> F(A,Y);

  
  matrix<Type> Z(A,Y);

  //M stuff
  switch(M_type){
  case fixedM:
    scalar_M(M,this);
    break;
  case Mdev:
    M_dev(M,nll,this);
    break;
  case LorMAP:
    Lorenzen_M_int(M,nll,this,growth_fun,start_age);
    break;
  case LorMAPdev:
    Lorenzen_M_int_dev(M,nll,this,growth_fun,start_age);
    break;
  default:
    error("M type not implemented!");
    break;
  };


  //make length bins + ages
  vector<Type> len_bins(L3);
  vector<Type> ages_jan(A);
  vector<Type> ages_fall(A);
  vector<Type> ages_mid(A);

  matrix<Type> len_sds(A,Y);
  matrix<Type> len_mus(A,Y);

  for(int l = 0; l < L3; ++l){
    len_bins(l) = 1+l+bin_adjust;
    //len_bins(l) = start_length+l+bin_adjust;
  }


  
  for(int a = 0; a < A; ++a){
    Type age_temp = Type(a)+Type(start_age);
    ages_jan(a) = age_temp;
    ages_mid(a) = age_temp+0.5;
    ages_fall(a) = age_temp+fall_adjust;
    for(int y = 0; y < Y; ++y){
      len_sds(a,y) = growth_fun.sd_len_at_age(ages_jan(a),y);
      len_mus(a,y) = growth_fun.mean_len_at_age(ages_jan(a),y);
    }
  }

  

  //Hold PLAs for every year 
  vector<matrix <Type> > plas_jan(Y);
  vector<matrix <Type> > plas_mid(Y);
  vector<matrix <Type> > plas_fall(Y);

  for(int y = 0; y < Y; ++y){
    plas_jan(y) = prob_len_at_age(growth_fun,len_bins,ages_jan,y);
    plas_mid(y) = prob_len_at_age(growth_fun,len_bins,ages_mid,y);
    plas_fall(y) = prob_len_at_age(growth_fun,len_bins,ages_fall,y);
  }

  vector< vector<Type> > ages_combo(3);
  ages_combo(0) = ages_jan;
  ages_combo(1) = ages_mid;
  ages_combo(2) = ages_fall;
  
  matrix<Type> S_ly;
  matrix<Type> S_ay;

  //F stuff
  switch(F_type){
  case sep_no_disc:
    rw_sep_F(F,nll,this,log_Fy,Fy_sd,log_Fa,Fa_sd);
    break;
  case AR:
    S_ly = mv_AR_sel(nll,this,Y,L3,len_bins);
    S_ay = convert_Sly_to_Say(S_ly,plas_jan,A);
    aggregated_F(S_ay,F,nll,this,log_Fy,Fy_sd);
    Report("S_ly",S_ly,this);
    Report("S_ay",S_ay,this);
    break;
  case rwrw:
    S_ly = rwrw_sel(L3,Y,len_bins,nll,this);
    S_ay = convert_Sly_to_Say(S_ly,plas_jan,A);
    aggregated_F(S_ay,F,nll,this,log_Fy,Fy_sd);
    Report("S_ly",S_ly,this);
    Report("S_ay",S_ay,this);
    break;
  case simp_len_blocked:
    S_ly = blocked_super_simple_sel(L3,Y,len_bins,nll,this);
    S_ay = convert_Sly_to_Say(S_ly,plas_jan,A);
    aggregated_F(S_ay,F,nll,this,log_Fy,Fy_sd);
    Report("S_ly",S_ly,this);
    Report("S_ay",S_ay,this);
    break;
  default:
    error("Incorrect F choice");
    break;
  };

  Z = F.array()+M.array();


  //N stuff
  matrix<Type> N(A,Y);

  switch(N_type){
  case fixedN:
    fixed_N(N,Z,nll,this,log_recruit);
    break;
  case randN:
    neo_cohort_effects_N(N,Z,nll,this,log_N_a,recruit_sd,F,S_ay,M,cfg);
    break;
  default:
    error("Incorrect N choice");
    break;
  };


  matrix<Type> N_mid(A,Y);
  matrix<Type> N_fall(A,Y);



  N_mid = N_adjust(N,Z,Type(0.5));
  N_fall = N_adjust(N,Z,fall_adjust);



  matrix<Type> NL(L3,Y);
  matrix<Type> NL_mid(L3,Y);
  matrix<Type> NL_fall(L3,Y);

  for(int y = 0; y < Y; ++y){
      NL.col(y) = plas_jan(y)*N.col(y);
      NL_mid.col(y) = plas_mid(y)*N_mid.col(y);
      NL_fall.col(y) = plas_fall(y)*N_fall.col(y);
  }


  //observation equations

  //Q stuff
  vector<Type> QL(L3);
  matrix<Type> QLM(L3,2);

 
      
	
  logisticQ_t<Type> camp_Q(this,"",L3,bin_adjust,1);
  logisticQ_t<Type> engel_Q(this,"engel.",L3,bin_adjust,1);
  QLM.col(1) = camp_Q.QL;
  QLM.col(0) = engel_Q.QL;
  vector<Type> Q_rho(L3);
  for(int l = 0; l < L3; ++l){
    Q_rho(l) = log(QLM(l,0)/QLM(l,1));
  }
      
  if(Q_prior == true){
    Type Q_max_sd = Parameter("log_Q_max_sd",this,"");
    Q_max_sd = exp(Q_max_sd);
    if(check_simulate(this)){
      for(int l = 0; l < L3; ++l){
      }
    }
    for(int l = 0; l < L3; ++l){
      nll -= dnorm(Q_rho(l),log(Qpriors(l)),Q_max_sd,true);
    }
    Report("QLM",QLM,this);
    Adreport("Q_max_sd",Q_max_sd,this,"",adr_flag);

	
  }
  camp_Q.report();
  engel_Q.report();
  Adreport("Q_rho",Q_rho,this,"",adr_flag);
  Report("Q_rho",Q_rho,this);

  fit_survey(nll,QLM,NL_fall,L3,Y,start_length,end_length,this);
  
  matrix<Type> catch_at_age(A,Y);

  for(int a = 0; a < A; ++a){
    for(int y = 0; y < Y; ++y){
      catch_at_age(a,y) = N_mid(a,y)*(1-exp(-Z(a,y)))*(F(a,y)/Z(a,y));
    }
  }

  matrix<Type> catch_at_length(L3,Y);
  for(int y = 0; y < Y; ++y){
    catch_at_length.col(y) = plas_mid(y)*catch_at_age.col(y);
  }
  
  matrix<Type> catch_biomass_mat = (weightsF.array()*catch_at_length.array());
  vector<Type> expected_landings = catch_biomass_mat.colwise().sum();
  vector<Type> log_expected_landings = log(expected_landings);



 

  switch(landing_type){
  case std_landing:
    fit_std_landings(nll,log_expected_landings,this);
    break;
  case cens_bound:
    fit_censored_landings(nll,log_expected_landings,this);
    break;
  default:
    error("Landing Type Not Implemented!");
  };
      
  switch(catch_type){
  case no_prop:
    break;
  case prop:
    fit_catch_prop(nll,catch_at_length,this);
    break;
  default:
    error("Catch Type not implemented");
  };
  
  
  
  matrix<Type> biomass_mat = weightsF.array()*NL.array();
  matrix<Type> ssb_mat = maturityF.array()*biomass_mat.array();


  Report("len_bins",len_bins,this);
  Report("ages_jan",ages_jan,this);
  Report("ages_fall",ages_fall,this);
  Report("ages_mid",ages_mid,this);

  //Adreport("sd_survey",sd_survey,this,"",adr_flag);
  Adreport("QLM",QLM,this,"",adr_flag);
  
  Report("plas_jan",plas_jan,this);
  Report("plas_mid",plas_mid,this);
  Report("plas_fall",plas_fall,this);
  Report("len_sds",len_sds,this);
  Report("len_mus",len_mus,this);
  
  Report("F",F,this);
  Adreport("F",F,this,"",adr_flag);
  Report("Z",Z,this);
  Adreport("Z",Z,this,"",adr_flag);
  Report("N",N,this);
  Report("N_mid",N_mid,this);
  Report("N_fall",N_fall,this);
  
  vector<Type> Fbar(Y);
  for(int y = 0; y < Y; ++y){
    Fbar(y) = 0;
    for(int a = 0; a < A; ++a){
      Fbar(y) += F(a,y);
    }
    Fbar(y) /= A;
  }
  Adreport("Fbar",Fbar,this,"",adr_flag);
  Report("Fbar",Fbar,this);

  vector<Type> log_Fbar = log(Fbar);
  Adreport("log_Fbar",log_Fbar,this,"",adr_flag);
  Report("log_Fbar",log_Fbar,this);
  
  Report("NL",NL,this);
  Report("NL_mid",NL_mid,this);
  Report("NL_fall",NL_fall,this);


  Report("log_expected_landings",log_expected_landings,this);
  Adreport("log_expected_landings",log_expected_landings,this,"",adr_flag);

  matrix<Type> log_N = N.array().log();
  matrix<Type> log_NL = NL.array().log();

  //They don't REALLY need to be reported
  Adreport("log_N",log_N,this,"",false);
  Adreport("log_NL",log_NL,this,"",false);
  Report("M",M,this);

  Report("catch_biomass_mat",catch_biomass_mat,this);
  Report("catch_at_length",catch_at_length,this);
  Report("catch_at_age",catch_at_age,this);
  Adreport("catch_at_age",catch_at_age,this,"",false);
  Adreport("catch_at_length",catch_at_length,this,"",false);
  Report("biomass_mat",biomass_mat,this);
  Adreport("biomass_mat",biomass_mat,this,"",false);
  Report("ssb_mat",ssb_mat,this);
  Adreport("ssb_mat",ssb_mat,this,"",false);

  vector<Type> tot_ssb = ssb_mat.colwise().sum();
  vector<Type> log_tot_ssb = ssb_mat.array().colwise().sum().log();
  Adreport("log_tot_ssb",log_tot_ssb,this,"",adr_flag);
  Adreport("tot_ssb",tot_ssb,this,"",adr_flag);

  vector<Type> log_tot_N = N.array().colwise().sum().log();
  Adreport("log_tot_N",log_tot_N,this,"",adr_flag);


  growth_fun.report();

  if(N_type == 1){
    log_recruit = log_N_a.row(0);
  }

  

  return nll;

}

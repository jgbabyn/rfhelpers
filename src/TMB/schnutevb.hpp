template<class Type>
Type weighted_mean(vector<Type> x,vector<Type> w){
  Type numsum = 0.0;
  Type denomsum = 0.0;
  for(int i = 0; i < x.size(); ++i){
    numsum += x(i)*w(i);
    denomsum += w(i);
  }
  Type xmean = numsum/denomsum;
  return xmean;
}

namespace my_expecteds {

  ///SD
template<class Float>
struct expsd_t {
  typedef Float Scalar;
  Float a_pg;
  Float lil_s;
  Float big_s;
  Float lil_l;
  Float big_L;
  Float M;

  Float operator() (Float x) {
    Float t = x+M;
    Float sd = lil_s+(big_s-lil_s)*((t-lil_l)/(big_L-lil_l));
    Float lambda = 1/a_pg;
    Float dexpp = lambda*exp(-lambda*x);
    Float ans = dexpp*sd;
    return ans;
  }

  Float expected() {
    using gauss_kronrod::integrate;
    Float ans = integrate(*this,0,INFINITY);
    return ans;
  }
  
};

  template<class Float>
  Float eval(Float a_pg,Float lil_s,Float big_s, Float lil_l,Float big_L,Float M){
    expsd_t<Float> f = {a_pg,lil_s,big_s,lil_l,big_L,M};
    return f.expected();
  }

  TMB_BIND_ATOMIC(func,111110,eval(x[0],x[1],x[2],x[3],x[4],x[5]));

  template<class Type>
  Type expected_pg_sd(Type a_pg,Type s,Type S, Type l, Type L, Type M){
    vector<Type> args(7);
    args << a_pg, s, S, l, L, M, 0;
    return my_expecteds::func(CppAD::vector<Type>(args))[0];
  }


  //Mu
  template<class Float>
  struct expmu_t {
    typedef Float Scalar;
    Float a_pg;
    Float lil_l;
    Float big_L;
    Float k;
    Float M;

  Float operator() (Float x) {
    Float t = x+M;
    Float mu = lil_l+(big_L-lil_l)*((1-pow(k,t-1))/(1-pow(k,M-1)));
    Float lambda = 1/a_pg;
    Float dexpp = lambda*exp(-lambda*x);
    Float ans = dexpp*mu;
    return ans;
  }

  Float expected() {
    using gauss_kronrod::integrate;
    Float ans = integrate(*this,0,INFINITY);
    return ans;
  }
  
};

  template<class Float>
  Float evalmu(Float a_pg,Float lil_l,Float big_L, Float k,Float M){
    expsd_t<Float> f = {a_pg,lil_l,big_L,k,M};
    return f.expected();
  }

  TMB_BIND_ATOMIC(funcmu,11110,evalmu(x[0],x[1],x[2],x[3],x[4]));

    template<class Type>
  Type expected_pg_mu(Type a_pg,Type l,Type L, Type k, Type M){
    vector<Type> args(6);
    args << a_pg, l, L, k, M, 0;
    return my_expecteds::funcmu(CppAD::vector<Type>(args))[0];
  }

  
}
  
template<class Type>
struct schnutevb_t{

  //Needed to make sure parameters get attached to obj fun.
  objective_function<Type> *objective;

  //needed to make sure things get put where they should be
  std::string prefix;

  //Max Age
  Type M;

  //Holds the type of growth function
  int gf_type;

  //Parameters for single time version
  //Parameters
  //First component length
  Type ell;
  //Last component length
  Type L;
  //Growth fraction
  Type k;

  //SD parameters
  Type s;
  Type S;

  //Plus group Z to account for survival
  Type pgZ;
  int pg_ext;
  
  vector<Type> pg_surv;

  Type init_a_pg;
  vector<Type> as_pg;

  matrix<Type> log_N_a;

  //Parameters for the time varying by years version
  //2 x Y matrix, first row is ell, L is 2nd
  matrix<Type> ells_Ls;
  matrix<Type> ss_Ss;
  vector<Type> ks;

  Type phi_ells;
  Type phi_Ls;
  Type phi_ks;
  Type phi_ss;
  Type phi_Ss;

  Type ells_sd;
  Type Ls_sd;
  Type ks_sd;
  Type ss_sd;
  Type Ss_sd;

  //Fixed L infinity
  Type f_L_inf;


  // constructors
  schnutevb_t(objective_function<Type> *obj, std::string pre,int A,int gftype,int N_type,matrix<Type> &log_N);
  //This empty one is used for the projections
  schnutevb_t();
  
  //Functions
  //Stuff to convert to classical VB
  Type get_K();
  Type get_L_infinity();
  Type get_t_zero(Type a1);
  void report();

  //Mean and sd
  Type mean_len_at_age(Type t);
  Type sd_len_at_age(Type t);
  Type mean_len_at_ageP(Type t);
  Type sd_len_at_ageP(Type t);

  Type sd_len_at_age2(Type t);
  Type mean_len_at_age(Type t,int y);
  Type sd_len_at_age(Type t,int y);

  void L_from_f_L_inf();

  Type log_lik(objective_function<Type> * obj);

  Type inverse_gf(Type mu);
  Type calc_x(Type j,Type k_dev);

};


#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<class Type>
schnutevb_t<Type>::schnutevb_t(objective_function<Type> *obj, std::string pre,int A,int gftype,int N_type,matrix<Type> &log_N){
  objective = obj;
  M = Type(A);
  prefix = pre;
  gf_type = gftype;

  switch(gf_type){
  case 0:
  case 3:
    {
      vector<Type> ell_L(2);

      //Keeps ell and L ordered
      Type log_L = Parameter("log_L",obj,pre);
      Type log_ell = Parameter("log_ell",obj,pre);
      ell_L(0) = log_ell;
      ell_L(1) = log_L;
      ell_L = exp(ordered_inv_transform(ell_L));

      //Keep s and S ordered
      vector<Type> s_S(2);
      //s_S(0) = Parameter("log_s",objective,pre);
      //s_S(1) = Parameter("log_S",objective,pre);
      Type log_s = Parameter("log_s",obj,pre);
      Type log_S = Parameter("log_S",obj,pre);
      s_S(0) = exp(log_s);
      s_S(1) = exp(log_S);

      //Parameters and transforms
      ell = ell_L(0);

      L = ell_L(1);

      Type log_k = Parameter("log_k",obj,pre);
      
      k = exp(log_k);

      S = s_S(1);

      s = s_S(0);
    }
    break;
  case 1:
  {
    matrix<Type> log_ells_Ls = ParameterMatrix("log_ells_Ls",obj,pre);
    ells_Ls = log_ells_Ls;
    matrix<Type> log_ss_Ss = ParameterMatrix("log_ss_Ss",obj,pre);
    ss_Ss = log_ss_Ss;
    vector<Type> log_ks = ParameterVector("log_ks",obj,pre);

    Type t_phi_ells = Parameter("t_phi_ells",obj,pre);
    Type t_phi_Ls = Parameter("t_phi_Ls",obj,pre);
    Type t_phi_ks = Parameter("t_phi_ks",obj,pre);
    Type t_phi_ss = Parameter("t_phi_ss",obj,pre);
    Type t_phi_Ss = Parameter("t_phi_Ss",obj,pre);


    Type log_ells_sd = Parameter("log_ells_sd",obj,pre);
    Type log_Ls_sd = Parameter("log_Ls_sd",obj,pre);
    Type log_ks_sd = Parameter("log_ks_sd",obj,pre);
    Type log_ss_sd = Parameter("log_ss_sd",obj,pre);
    Type log_Ss_sd = Parameter("log_Ss_sd",obj,pre);

    //Each column is ordered
    for(int i = 0; i < ells_Ls.cols(); ++i){
      vector<Type> temp(2);
      temp = log_ells_Ls.col(i);
      ells_Ls.col(i) = exp(ordered_inv_transform(temp));
      temp = log_ss_Ss.col(i);
      ss_Ss.col(i) = exp(ordered_inv_transform(temp));
    }

    ks = exp(log_ks);

    phi_ells = invlogit(t_phi_ells);
    phi_Ls = invlogit(t_phi_Ls);
    phi_ks = invlogit(t_phi_ks);
    phi_ss = invlogit(t_phi_ss);
    phi_Ss = invlogit(t_phi_Ss);

    ells_sd = exp(log_ells_sd);
    Ls_sd = exp(log_Ls_sd);
    ks_sd = exp(log_ks_sd);
    ss_sd = exp(log_ss_sd);
    Ss_sd = exp(log_Ss_sd);
  }
  break;
    case 2:
    {
      vector<Type> ell_L(2);

      //Keeps ell and L ordered
      //PARAMETER(log_L);
      

      
      PARAMETER(log_ell);
      ell_L(0) = log_ell;
      ell_L(1) = log_ell+0.001;
      //ell_L(1) = log_L;

      // ell_L(0) = Parameter("log_ell",objective,pre);
      // ell_L(1) = Parameter("log_L",objective,pre);
      ell_L = exp(ordered_inv_transform(ell_L));

      //Keep s and S ordered
      vector<Type> s_S(2);
      //s_S(0) = Parameter("log_s",objective,pre);
      //s_S(1) = Parameter("log_S",objective,pre);
      PARAMETER(log_s);
      PARAMETER(log_S);
      s_S(0) = exp(log_s);
      s_S(1) = exp(log_S);

      //s_S = exp(ordered_inv_transform(s_S));


      //Parameters and transforms
      ell = ell_L(0);

      DATA_SCALAR(fixed_L_inf);
      f_L_inf = fixed_L_inf;

      //L = ell_L(1);

      PARAMETER(log_k);
      //k = Parameter("log_k",objective,pre);
      k = exp(log_k);

      S = s_S(1);

      s = s_S(0);
      //Calculate the L from the fixed L_inf.
      L_from_f_L_inf();

    }
    break;
  case 4:
    {
      vector<Type> ell_L(2);

      //Keeps ell and L ordered

      Type log_L = Parameter("log_L",obj,pre);
      Type log_ell = Parameter("log_ell",obj,pre);
  
      ell_L(0) = log_ell;
      ell_L(1) = log_L;

      // ell_L(0) = Parameter("log_ell",objective,pre);
      // ell_L(1) = Parameter("log_L",objective,pre);
      ell_L = exp(ordered_inv_transform(ell_L));

      //Keep s and S ordered
      vector<Type> s_S(2);
      //s_S(0) = Parameter("log_s",objective,pre);
      //s_S(1) = Parameter("log_S",objective,pre);
      Type log_s = Parameter("log_s",obj,pre);
      Type log_S = Parameter("log_S",obj,pre);

      s_S(0) = exp(log_s);
      s_S(1) = exp(log_S);

      pg_ext = DataInteger("pg_ext",obj,"");

      Type log_init_a_pg = Parameter("log_init_a_pg",obj,pre);
      //Upper bounded by the max age.
      init_a_pg = ub_inv(log_init_a_pg,Type(pg_ext));
  

      if(N_type == 1){
	log_N_a = log_N;
	//matrix<Type> log_N_a = glanceParameterMatrix("log_N_a",obj,"");
	vector<Type> _as_pg(log_N_a.cols());
	pgZ = -(log_N_a(A-1,0)-log_N_a(A-2,0));

	vector<Type> pg_s_t(pg_ext);
	vector<Type> pg_ag(pg_ext);
	pg_surv = pg_s_t;
	pg_surv(0) = exp(-pgZ);
	pg_ag(0) = 0;
	for(int i = 1; i < pg_surv.size(); ++i){
	  pg_ag(i) = i;
	  pg_surv(i) = exp(-pgZ)*pg_surv(i-1);
	}
	//init_a_pg = weighted_mean(pg_ag,pg_surv);
	
	as_pg = _as_pg;
	as_pg(0) = init_a_pg;
	for(int y = 1; y < log_N_a.cols(); ++y){
	  as_pg(y) = (exp(log_N_a(A-1,y-1))*1+exp(log_N_a(A,y-1))*(as_pg(y-1)+1))/(exp(log_N_a(A-1,y-1))+exp(log_N_a(A,y-1)));
	 }
      }else{
	error("This doesn't work lol!");
      }

      
      //s_S = exp(ordered_inv_transform(s_S));


      //Parameters and transforms
      ell = ell_L(0);

      L = ell_L(1);

      //PARAMETER(log_k);
      Type log_k = Parameter("log_k",objective,pre);
      k = exp(log_k);

      S = s_S(1);

      s = s_S(0);
    }
    break;
  case 5:
    {
      vector<Type> ell_L(2);
      
      //Keeps ell and L ordered

      Type log_L = Parameter("log_L",obj,pre);
      Type log_ell = Parameter("log_ell",obj,pre);
  
      ell_L(0) = log_ell;
      ell_L(1) = log_L;

      // ell_L(0) = Parameter("log_ell",objective,pre);
      // ell_L(1) = Parameter("log_L",objective,pre);
      ell_L = exp(ordered_inv_transform(ell_L));

      //Keep s and S ordered
      vector<Type> s_S(2);
      //s_S(0) = Parameter("log_s",objective,pre);
      //s_S(1) = Parameter("log_S",objective,pre);
      Type log_s = Parameter("log_s",obj,pre);
      Type log_S = Parameter("log_S",obj,pre);

      s_S(0) = exp(log_s);
      s_S(1) = exp(log_S);

      pg_ext = DataInteger("pg_ext",obj,"");

      //Type log_init_a_pg = Parameter("log_init_a_pg",obj,pre);
      //Upper bounded by the max age.
      //init_a_pg = ub_inv(log_init_a_pg,Type(pg_ext));
      
  

      if(N_type == 1){
	log_N_a = log_N;
	//matrix<Type> log_N_a = glanceParameterMatrix("log_N_a",obj,"");
	vector<Type> _as_pg(log_N_a.cols());
	pgZ = -(log_N_a(A-1,0)-log_N_a(A-2,0));

	vector<Type> pg_s_t(pg_ext);
	vector<Type> pg_ag(pg_ext);
	pg_surv = pg_s_t;
	pg_surv(0) = exp(-pgZ);
	pg_ag(0) = 0;
	for(int i = 1; i < pg_surv.size(); ++i){
	  pg_ag(i) = i;
	  pg_surv(i) = exp(-pgZ)*pg_surv(i-1);
	}
	init_a_pg = weighted_mean(pg_ag,pg_surv);
	
	as_pg = _as_pg;
	as_pg(0) = init_a_pg;
	for(int y = 1; y < log_N_a.cols(); ++y){
	  as_pg(y) = (exp(log_N_a(A-1,y-1))*1+exp(log_N_a(A,y-1))*(as_pg(y-1)+1))/(exp(log_N_a(A-1,y-1))+exp(log_N_a(A,y-1)));
	 }
      }else{
	error("This doesn't work lol!");
      }

      
      //s_S = exp(ordered_inv_transform(s_S));


      //Parameters and transforms
      ell = ell_L(0);

      L = ell_L(1);

      //PARAMETER(log_k);
      Type log_k = Parameter("log_k",objective,pre);
      k = exp(log_k);

      S = s_S(1);

      s = s_S(0);
    }
    break;
 
  default:
    error("GF type not implemented?");
  };


}
  #undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this


  template<class Type>
  Type schnutevb_t<Type>::get_K(){
    Type K = -log(k);
    return K;
  }

  template<class Type>
  Type schnutevb_t<Type>::get_L_infinity(){
    Type L_inf_top = L-ell*pow(k,M-1);
    Type L_inf_bot = 1-pow(k,M-1);
    Type L_inf = L_inf_top/L_inf_bot;
    return L_inf;
  }

  template<class Type>
  Type schnutevb_t<Type>::get_t_zero(Type a1){
    Type part = (L-ell)/(L-ell*pow(k,M-1));
    Type t_zero = a1-(1/log(k))*log(part);
    return t_zero;
  }

 template<class Type>
 Type schnutevb_t<Type>::mean_len_at_ageP(Type t){
   Type mu_t = ell+(L-ell)*((1-pow(k,t-1))/(1-pow(k,M-1)));
   return mu_t;
 }

  template<class Type>
  Type schnutevb_t<Type>::mean_len_at_age(Type t){
    Type mu_t = ell+(L-ell)*((1-pow(k,t-1))/(1-pow(k,M-1)));
    if(t >= M){
      vector<Type> pg_mus(pg_ext);
      for(int b = 0; b < pg_ext; ++b){
	Type age_p_temp = t+Type(b);
	pg_mus(b) = ell+(L-ell)*((1-pow(k,age_p_temp-1))/(1-pow(k,M-1)));
      }
      mu_t = weighted_mean(pg_mus,pg_surv);
    }
    return mu_t;
  }

template<class Type>
Type schnutevb_t<Type>::mean_len_at_age(Type t,int y){
  switch(gf_type){
  case 0:
  case 2:  
    return mean_len_at_ageP(t);
    break;
  case 1:
    {
    Type ellt = ells_Ls(0,y);
    Type Lt = ells_Ls(1,y);
    Type kt = ks(y);
    Type mu_t = ellt+(Lt-ellt)*((1-pow(kt,t-1))/(1-pow(kt,M-1)));
    return mu_t;
    }
    break;
  case 3:
    return mean_len_at_ageP(t);
    break;
  case 4:
  case 5:
   if(t >= M){
     vector<Type> weights(pg_ext);
     vector<Type> mus(pg_ext);
    Type tdff = t-M;
    for(int i = 0; i < pg_ext; ++i){
      weights(i) = dexp(i+tdff,1/as_pg(y));
      mus(i) = mean_len_at_ageP(t+Type(i)); 
    }
    Type wmmu = weighted_mean(mus,weights);
    //Type intmu = my_expecteds::expected_pg_mu(as_pg(y),ell,L,k,M);
    //Rcout << "wmmu:" << wmmu << std::endl;
    //Rcout << "intmu:" << intmu << std::endl;
    return wmmu;
    //return intmu;
  
    }else{
      return mean_len_at_ageP(t);
     }
    break;
  default:
    error("Not implemented GF function");
    break;
  };
}

  // template<class Type>
  // Type schnutevb_t<Type>::sd_len_at_age(Type t){
  //   Type mu_t = mean_len_at_age(t);
  //   Type sd_t = s+(S-s)*((mu_t-ell)/(L-ell));
  //   return sd_t;
  // }


// template<class Type>
// Type schnutevb_t<Type>::sd_len_at_ageP(Type t){
//   // Type mu_t = mean_len_at_age(t);
//   Type age = t+1;
//   Type sd_t = s+(S-s)*((t-ell)/(L-ell));
//   return sd_t;
// }

template<class Type>
Type schnutevb_t<Type>::sd_len_at_ageP(Type t){
  // Type mu_t = mean_len_at_age(t);
  Type age = t+1;
  Type sd_t = s+(S-s)*((t-ell)/(L-ell));
  return sd_t;
}


//Based on age rather than length
  template<class Type>
  Type schnutevb_t<Type>::sd_len_at_age(Type t){
    Type mu_t = mean_len_at_age(t);
    Type age = t+1;
    Type sd_t = s+(S-s)*((t-ell)/(L-ell));
    if(t >= M){
      vector<Type> pg_sds(pg_ext);
      for(int b = 0; b < pg_ext; ++b){
	Type age_p_temp = t+Type(b);
	pg_sds(b) = s+(S-s)*((age_p_temp-ell)/(L-ell));
      }
      sd_t = weighted_mean(pg_sds,pg_surv);
    }
    return sd_t;
  }


//Increase the standard deviation for the plus group..
template<class Type>
Type schnutevb_t<Type>::sd_len_at_age2(Type t){
  Type sd_t = 0;
  if(t < M){
    sd_t = s+(S-s)*((t-ell)/(L-ell));
  }else{
    //FoR NOW
    sd_t = 2*S;
  }
  return sd_t;
}
 
template<class Type>
Type schnutevb_t<Type>::sd_len_at_age(Type t,int y){
  switch(gf_type){
  case 0:
  case 2:  
    return sd_len_at_ageP(t);
    break;
  case 1:
    {
    Type ellt = ells_Ls(0,y);
    Type Lt = ells_Ls(1,y);
    Type kt = ks(y);
    Type mu_t = mean_len_at_age(t,y);
    Type st = ss_Ss(0,y);
    Type St = ss_Ss(1,y);
    Type sd_t = st+(St-st)*((mu_t-ellt)/(Lt-ellt));
    return sd_t;
    }
    break;
  case 3:
    return sd_len_at_ageP(t);
    break;
  case 4:
  case 5:  
     if(t >= M){
       vector<Type> weights(pg_ext);
       vector<Type> sds(pg_ext);
       Type tdff = t-M;
       for(int i = 0; i < pg_ext; ++i){
	 weights(i) = dexp(i+tdff,1/as_pg(y));
	 sds(i) = sd_len_at_ageP(t+Type(i)); 
       }
       Type wmsd = weighted_mean(sds,weights);
       //Type intsd = my_expecteds::expected_pg_sd(as_pg(y),s,S,ell,L,M);
       //Rcout << "wmsd:" << wmsd << std::endl;
       //Rcout << "intsd:" << intsd << std::endl;
       return wmsd;
       //return intsd;
    }else{
      return sd_len_at_ageP(t);
     }
  default:
    error("Not implemented GF function");
    break;
  };
}


template<class Type>
schnutevb_t<Type>::schnutevb_t(){
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<class Type>
Type schnutevb_t<Type>::log_lik(objective_function<Type> * obj){
  switch(gf_type){
  case 0:
  case 2:
  case 3:
  case 4:  
    return 0.0;
    break;
  case 1:
    {
      Type pnll = 0.0;
      vector<Type> tempv = ells_Ls.array().row(0).log();
      //We want neg cause we subtract
      pnll -= density::SCALE(density::AR1(phi_ells),ells_sd)(tempv);
      tempv = ells_Ls.array().row(0).log();
      pnll -= density::SCALE(density::AR1(phi_Ls),Ls_sd)(tempv);
      pnll -= density::SCALE(density::AR1(phi_ks),ks_sd)(log(ks));
      return pnll;

    }
    break;
  default:
    error("Not implemented GF function");
    break;
  };
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this


  template<class Type>
  Type schnutevb_t<Type>::inverse_gf(Type mu){
    Type Z1 = (1-(((mu-ell)/(L-ell))*(1-pow(k,(M-1)))));
    return log(Z1)/log(k)+1;
  }

  template<class Type>
  Type schnutevb_t<Type>::calc_x(Type j, Type k_dev){
    Type mean = mean_len_at_age(j);
    Type sd = sd_len_at_age(j);
    Type inv = inverse_gf(mean+k_dev*sd);
    //Type x = inv - j;
    //return x;
    return inv;
  }


#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR objective
template<class Type>
  void schnutevb_t<Type>::report(){
  REPORT(ell);
  REPORT(L);
  REPORT(k);
  REPORT(S);
  REPORT(s);
  REPORT(pgZ);
  ADREPORT(ell);
  ADREPORT(L);
  ADREPORT(k);
  ADREPORT(S);
  ADREPORT(s);
  if(gf_type == 4){
    Type init_a_pg_t = init_a_pg+M;
    Adreport("init_a_pg",init_a_pg_t,objective,"");
    Report("init_a_pg",init_a_pg_t,objective,"");
    REPORT(as_pg);
    REPORT(pg_surv);
    REPORT(pgZ);
  }
  Type K = get_K();
  Type Linf = get_L_infinity();
  Type tzero = get_t_zero(Type(1.0));
  Adreport("vbK",K,objective,"growth.");
  Adreport("vbLinf",Linf,objective,"growth.");
  Adreport("vbtzero",tzero,objective,"growth.");
    // Report("ell",ell,objective,prefix);
    // Report("L",L,objective,prefix);
    // Report("k",k,objective,prefix);
    // Report("S",S,objective,prefix);
    // Report("s",s,objective,prefix);
  }
  #undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

template<class Type>
void schnutevb_t<Type>::L_from_f_L_inf(){
  L = f_L_inf*(1-pow(k,M-1))+ell*pow(k,M-1);

}

//' @name prob_len_at_age
//' @description Probability of Length at Age assuming a VB growth function and normal distribution assumptions
//' @param growth_fun the VB growth function
//' @param length_bins the length bins used by the model
//' @param ages the vector of ages to calculate the probability of length at age of
//' @param y the current year to use for the growth model
template<class Type>
matrix<Type> prob_len_at_age(schnutevb_t<Type> &growth_fun,vector<Type> &length_bins,vector<Type> &ages,int y){
  int A = ages.size();
  int L = length_bins.size();

  matrix<Type> pla(L,A);

  for(int a = 0; a < A; ++a){
    Type mu_a = growth_fun.mean_len_at_age(ages(a),y);
    Type sd_a = growth_fun.sd_len_at_age(ages(a),y);


    for(int l = 0; l < L; ++l){
      if(l == 0){
	pla(l,a) = pnorm(length_bins(l),mu_a,sd_a);
      }else if(l == L-1){
	pla(l,a) = 1-pnorm(length_bins(l-1),mu_a,sd_a);
      }else{
	pla(l,a) = pnorm(length_bins(l),mu_a,sd_a)-pnorm(length_bins(l-1),mu_a,sd_a);
      }
    }


  }


  return pla;
}


//This is so gross
template<class Type>
void lognormal_recruits(vector<Type> &log_recruit,Type &nll, objective_function<Type> *obj){
  Type recruit_sd = Parameter("log_recruit_sd",obj,"");
  int adr_flag = DataInteger("adr_flag",obj);
  recruit_sd = exp(recruit_sd);
  for(int i = 1; i < log_recruit.size(); ++i){
    nll -= dnorm(log_recruit(i),log_recruit(i-1),recruit_sd,true);
  }
  
  Adreport("recruit_sd",recruit_sd,obj,"",adr_flag);
}


  
template<class Type>
matrix<Type> cohort_dynamics(matrix<Type> &N, matrix<Type> &Z){
  int A = N.rows();
  int Y = N.cols();
  for(int a = 1; a < A; ++a){
    for(int y = 1; y < Y; ++y){
      N(a,y) = N(a-1,y-1)*exp(-Z(a-1,y-1));
      if(a == A-1){
	N(a,y) += N(a,y-1)*exp(-Z(a,y-1));
      }
    }
  }
  return N;
}

template<class Type>
matrix<Type> N_adjust(matrix<Type> &N, matrix<Type> &Z, Type adjust_time){
  int A = N.rows();
  int Y = N.cols();

  matrix<Type> adj_N = N;
  for(int a = 0; a < A; ++a){
    for(int y = 0; y < Y; ++y){
      adj_N(a,y) = N(a,y)*exp(-Z(a,y)*adjust_time);
    }
  }
  return adj_N;
}

template<class Type>
void fixed_N(matrix<Type> &N, matrix<Type> &Z, Type &nll, objective_function<Type> *obj,vector<Type> log_recruit){
  int A = N.rows();
  int Y = N.cols();

  int adr_flag = DataInteger("adr_flag",obj);
 

  vector<Type> log_N0 = ParameterVector("log_N0",obj);
  Type log_N0_sd = Parameter("log_N0_sd",obj);
  Type N0_sd = exp(log_N0_sd);
  //currently no good
  //int rec_type = 0;
  
  //Simulate recruits
  // if(check_simulate(obj)){
  //   log_recruit(0) = rnorm(log_recruit(0),recruit_sd);

  //   for(int y = 1; y < log_recruit.size(); ++y){
  //     log_recruit(y) = rnorm(log_recruit(y-1),recruit_sd);
  //   }
  //   Report("log_recruit",log_recruit,obj);
  // }

  vector<Type> goose_recruit(log_recruit.size());
  
  lognormal_recruits(log_recruit,nll,obj);
  goose_recruit = log_recruit;
  
  N.row(0) = exp(goose_recruit);

  //Simulate intial numbers at age
  if(check_simulate(obj)){
    vector<Type> ratios(A-1);
    for(int a = 1; a < A; ++a){
      ratios(a-1) = rnorm(-Z(a-1,0),N0_sd);
      if(a == 1){
	log_N0(a-1) = goose_recruit(0)+ratios(a-1);
      }else{
	log_N0(a-1) = log_N0(a-2)+ratios(a-1);
      }
    }
    Report("log_N0",log_N0,obj);
  }
  
  for(int a = 1; a < A; ++a){
    N(a,0) = exp(log_N0(a-1));
  }



  N = cohort_dynamics(N,Z);

  for(int i = 1; i < log_N0.size()-1; ++i){
    nll -= dnorm(log_N0(i)-log_N0(i-1),-Z(i-1,0),N0_sd,true);
  }
  int At = log_N0.size()-1;
  nll -= dnorm(log_N0(At)-log_N0(At-1),-log(1-exp(-Z(At-1,0))),N0_sd,true);
  



  


  vector<Type> recruit = exp(goose_recruit);
  Report("recruit",recruit,obj);
  vector<Type> N0 = exp(log_N0);
  Report("N0",N0,obj);
  Adreport("log_recruit",goose_recruit,obj,"",adr_flag);
  //REPORT(log_recruit);
  Adreport("N0_sd",N0_sd,obj,"",adr_flag);

}

template<class Type>
void rw_sep_F(matrix<Type> &F,Type &nll, objective_function<Type> *obj,vector<Type> &log_Fy,Type &Fy_sd,vector<Type> &log_Fa,Type &Fa_sd){
  int A = F.rows();
  int Y = F.cols();
  int adr_flag = DataInteger("adr_flag",obj);
 
  for(int i = 1; i < log_Fa.size(); ++i){
    nll -= dnorm(log_Fa(i),log_Fa(i-1),Fa_sd,true);
  }

  //Simulate F_a
  if(check_simulate(obj)){
    log_Fa(0) = rnorm(log_Fa(0),Fa_sd);

    for(int i = 1; i < log_Fa.size(); ++i){
      log_Fa(i) = rnorm(log_Fa(i-1),Fa_sd);
    }
    Report("log_Fa",log_Fa,obj);
  }
  
  for(int i = 1; i < log_Fy.size(); ++i){
    nll -= dnorm(log_Fy(i),log_Fy(i-1),Fy_sd,true);
  }

  //Simulate log_Fy
  if(check_simulate(obj)){
    log_Fy(0) = rnorm(log_Fy(0),Fy_sd);

    for(int i = 1; i < log_Fy.size(); ++i){
      log_Fy(i) = rnorm(log_Fy(i-1),Fy_sd);
    }
    Report("log_Fy",log_Fy,obj);
  }


  Adreport("Fa_sd",Fa_sd,obj,"",adr_flag);
  Adreport("Fy_sd",Fy_sd,obj,"",adr_flag);
  vector<Type> Fa = exp(log_Fa);
  vector<Type> Fy = exp(log_Fy);
  Adreport("Fa",Fa,obj,"",adr_flag);
  Adreport("Fy",Fy,obj,"",adr_flag);


  for(int a = 0; a < A; ++a){
    for(int y = 0; y < Y; ++y){
      if( a > 1){
	F(a,y) = exp(log_Fa(a))+exp(log_Fy(y));
	}else{
	F(a,y) = 0;
	}
    }
  } 
}

//Doesn't set F for ages 1 and 2 to be zero
template<class Type>
void rw_sep_F2(matrix<Type> &F,Type &nll, objective_function<Type> *obj){

  int A = F.rows();
  int Y = F.cols();

  vector<Type> log_Fa = ParameterVector("log_Fa",obj,"");
  Type Fa_sd = Parameter("log_Fa_sd",obj,"");
  Fa_sd = exp(Fa_sd);

  vector<Type> log_Fy = ParameterVector("log_Fy",obj,"");
  Type Fy_sd = Parameter("log_Fy_sd",obj,"");
  Fy_sd = exp(Fy_sd);

  for(int i = 1; i < log_Fy.size(); ++i){
    nll -= dnorm(log_Fy(i),log_Fy(i-1), Fy_sd,true);

  }

  for(int a = 1; a < log_Fa.size(); ++a){
    Type temp = exp(log_Fa(a-1));
    nll -= dnorm(log_Fa(a),log(temp),Fa_sd,true);
  }

  
  vector<Type> Fa = exp(log_Fa);
  vector<Type> Fy = exp(log_Fy);

  
  for(int a = 0; a < A; ++a){
    for(int y = 0; y < Y; ++y){
      F(a,y) = exp(log_Fa(a))+exp(log_Fy(y));
    }
  }

  Adreport("Fa",Fa,obj,"");
  Adreport("Fy",Fy,obj,"");
  Adreport("Fy_sd",Fy_sd,obj,"");
  Adreport("Fa_sd",Fa_sd,obj,"");

  
}

template<class Type>
void mv_AR_F(matrix<Type> &F,Type &nll, objective_function<Type> *obj){
  using namespace density;
  
  array<Type> array_F = ParameterArray("array_F",obj,"");
  Type phi_F_age = Parameter("t_phi_F_age",obj,"");
  phi_F_age = invlogit(phi_F_age);
  Type phi_F_year = Parameter("t_phi_F_year",obj,"");
  phi_F_year = invlogit(phi_F_year);
  Type sd_F = Parameter("log_sd_F",obj,"");
  sd_F = exp(sd_F);

  int A = F.rows();
  int Y = F.cols();

  for(int a = 0; a < A; ++a){
    for(int y = 0; y < Y; ++y){
      F(a,y) = exp(array_F(a,y));
    }
  }

  nll += SCALE(SEPARABLE(AR1(phi_F_year),AR1(phi_F_age)),sd_F)(array_F);

  
}

template<class Type>
matrix<Type> mv_AR_sel(Type &nll, objective_function<Type> *obj,int Y,int L, vector<Type> len_bins){
    using namespace density;
  
  array<Type> array_F = ParameterArray("array_F",obj,"");
  Type phi_F_age = Parameter("t_phi_F_age",obj,"");
  phi_F_age = invlogit(phi_F_age);
  Type phi_F_year = Parameter("t_phi_F_year",obj,"");
  phi_F_year = invlogit(phi_F_year);
  Type sd_F = Parameter("log_sd_F",obj,"");
  sd_F = exp(sd_F);


  matrix<Type> S_ly(L,Y);
  
  for(int l = 0; l < L; ++l){
    for(int y = 0; y < Y; ++y){
      S_ly(l,y) = invlogit(array_F(l,y));
    }
  }

  nll += SCALE(SEPARABLE(AR1(phi_F_year),AR1(phi_F_age)),sd_F)(array_F);

  return S_ly;
}
  
template<class Type>
void scalar_M(matrix<Type> &M,objective_function<Type> *obj){
  Type base_M = DataScalar("base_M",obj);

  int A = M.rows();
  int Y = M.cols();

  for(int a = 0; a < A; ++a){
    for(int y = 0; y < Y; ++y){
      M(a,y) = base_M;
    }
  }
}



template<class Type>
void M_dev(matrix<Type> &M, Type &nll, objective_function<Type> *obj){
  using namespace density;

  array<Type> process_e = ParameterArray("process_e",obj,"");
  Type phi_pe_year = Parameter("t_phi_pe_year",obj,"");
  phi_pe_year = invlogit(phi_pe_year);
  Type phi_pe_age = Parameter("t_phi_pe_age",obj,"");
  phi_pe_age = invlogit(phi_pe_age);
  Type sd_pe = Parameter("log_sd_pe",obj,"");
  sd_pe = exp(sd_pe);
  Type base_M = DataScalar("base_M",obj,"");

  
  int A = M.rows();
  int Y = M.cols();
  if(check_simulate(obj)){
    //SEPARABLE(AR1(phi_pe_year),AR1(phi_pe_age)).simulate(process_e);
    N01<Type> farts_age;
    N01<Type> farts_year;
    AR1_t<N01<Type> > pe_age(phi_pe_age,farts_age);
    AR1_t<N01<Type> > pe_year(phi_pe_year,farts_year);
    //You can't simulate from the scaled version... (https://github.com/kaskr/adcomp/issues/314)
    SEPARABLE_t< AR1_t<N01<Type> >, AR1_t<N01<Type> > > sep_age_year(pe_year,pe_age);
    sep_age_year.simulate(process_e);
    process_e = process_e/sd_pe;
    Report("process_e",process_e,obj);
  }
  
  
  for(int a = 0; a < A; ++a){
    for(int y = 0; y < Y; ++y){
      M(a,y) = exp(process_e(a,y));
    }
  }

  nll += SCALE(SEPARABLE(AR1(phi_pe_year),AR1(phi_pe_age)),sd_pe)(process_e);

  

  Adreport("phi_pe_year",phi_pe_year,obj,"");
  Adreport("phi_pe_age",phi_pe_age,obj,"");
  Adreport("sd_pe",sd_pe,obj,"");
  Adreport("process_e",process_e,obj,"");


}


template<class Type>
void Lorenzen_M_int(matrix<Type> &M, Type &nll, objective_function<Type> *obj, schnutevb_t<Type> &growth_fun,int start_age){
  Type base_M = DataScalar("base_M",obj,"");

  int A = M.rows();
  int gft = DataInteger("GF_type",obj,"");

  Type f_age = Type(A-1)+Type(start_age)+0.5;
  Type mu_f = 0.0;
  if(gft == 0){
    mu_f = growth_fun.mean_len_at_age(f_age);
  }else{
    mu_f = growth_fun.mean_len_at_ageP(f_age);
  } 
  int Y = M.cols();
  for(int a = 0; a < A; ++a){
    for(int y = 0; y < Y; ++y){
      //Mid year age
      Type age = Type(a)+Type(start_age) + 0.5;
      Type mu_len = 0.0;
      if(gft == 0){
	mu_len = growth_fun.mean_len_at_age(age);
      }else{
	mu_len = growth_fun.mean_len_at_ageP(age);
      }
      M(a,y) = base_M*(pow(mu_len,-1)/pow(mu_f,-1));
    }
  }

}

template<class Type>
void Lorenzen_M_int_dev(matrix<Type> &M, Type &nll, objective_function<Type> *obj, schnutevb_t<Type> &growth_fun,int start_age){
  using namespace density;

  array<Type> process_e = ParameterArray("process_e",obj,"");
  Type phi_pe_year = Parameter("t_phi_pe_year",obj,"");
  phi_pe_year = invlogit(phi_pe_year);
  Type phi_pe_age = Parameter("t_phi_pe_age",obj,"");
  phi_pe_age = invlogit(phi_pe_age);
  Type sd_pe = Parameter("log_sd_pe",obj,"");
  sd_pe = exp(sd_pe);
  Type base_M = DataScalar("base_M",obj,"");

  
  int A = M.rows();
  int Y = M.cols();
  if(check_simulate(obj)){
    //SEPARABLE(AR1(phi_pe_year),AR1(phi_pe_age)).simulate(process_e);
    N01<Type> farts_age;
    N01<Type> farts_year;
    AR1_t<N01<Type> > pe_age(phi_pe_age,farts_age);
    AR1_t<N01<Type> > pe_year(phi_pe_year,farts_year);
    //You can't simulate from the scaled version... (https://github.com/kaskr/adcomp/issues/314)
    SEPARABLE_t< AR1_t<N01<Type> >, AR1_t<N01<Type> > > sep_age_year(pe_year,pe_age);
    sep_age_year.simulate(process_e);
    process_e = process_e/sd_pe;
    Report("process_e",process_e,obj);
  }
  
  Type f_age = Type(A-1)+Type(start_age)+0.5;
  Type mu_f = growth_fun.mean_len_at_age(f_age);
  for(int a = 0; a < A; ++a){
    for(int y = 0; y < Y; ++y){
      Type age = Type(a)+Type(start_age) + 0.5;
      Type mu_len = growth_fun.mean_len_at_age(age);
      M(a,y) = exp(process_e(a,y))*(pow(mu_len,-1)/pow(mu_f,-1));
    }
  }

  nll += SCALE(SEPARABLE(AR1(phi_pe_year),AR1(phi_pe_age)),sd_pe)(process_e);

  Adreport("phi_pe_year",phi_pe_year,obj,"");
  Adreport("phi_pe_age",phi_pe_age,obj,"");
  Adreport("sd_pe",sd_pe,obj,"");
  Adreport("process_e",process_e,obj,"");


}



//Construct an AYxAY sigma matrix with ind. entries
template<class Type>
matrix<Type> ind_N_sigma(int A,int Y,Type recruit_sd, Type N0_sd, Type surv_sd,objective_function<Type> *obj){

  vector<int> vagesindices = DataIVector("vagesindices",obj,"");
  vector<int> vyearsindices = DataIVector("vyearsindices",obj,"");

  
  int AY = A*Y;
  matrix<Type> N_sigma(AY,AY);

  
  
  for(int i = 0; i < AY; ++i){
    for(int j = 0; j <= i; ++j){
      if(i != j){
	N_sigma(i,j) = 0;
	N_sigma(j,i) = 0;
      }else{
	if(vagesindices[i] == 0){
	  N_sigma(i,j) = pow(recruit_sd,2);
	}else if(vyearsindices[i] == 0){
	  N_sigma(i,j) = pow(N0_sd,2);
	}else{
	  N_sigma(i,j) = pow(surv_sd,2);
	}
      }
    }
  }

  return N_sigma;

}


//Build beta for logistic selectivity at length
template<class Type>
Type build_beta(Type lower_b,Type upper_b,Type logit_mean,
				Type eta_y){
  Type top = upper_b-lower_b;
  Type bot = 1+exp(-(logit_mean+eta_y));
  Type beta = lower_b+top/bot;
  return beta;
}

//Build the logistic length selectivity matrix elements
template<class Type>
Type make_Sly(Type l,Type beta_1y,Type beta_2y){
  //Type inty = (-l-beta_1y)/beta_2y;
  //Type Sly = 1/(1+exp(inty));
  //Type interior = -beta_2y*(l-beta_1y);
  Type interior = -(1/beta_2y)*(l-beta_1y);
  Type bot = (1+exp(interior));
  Type Sly = 1/bot;
  
  return Sly;
}

template<class Type>
matrix<Type> build_logistic_len_sel(Type &nll, objective_function<Type> *obj,int Y,int L, vector<Type> len_bins){
  using namespace density;

  Type logit_m_1 = Parameter("logit_m_1",obj,"");
  Type logit_m_2 = Parameter("logit_m_2",obj,"");

  Type lower_b1 = DataScalar("lower_b1",obj,"");
  Type upper_b1 = DataScalar("upper_b1",obj,"");

  Type lower_b2 = DataScalar("lower_b2",obj,"");
  Type upper_b2 = DataScalar("upper_b2",obj,"");

  //Random effects
  Type sd_sel = Parameter("log_sd_sel",obj,"");
  sd_sel = exp(sd_sel);
  Type rho_par = Parameter("logit_rho_par",obj,"");
  rho_par = invlogit(rho_par);
  Type rho_year = Parameter("logit_rho_year",obj,"");
  rho_year = invlogit(rho_year);
  
  //vector<Type> eta_y1 = ParameterVector("eta_y1",obj,"");
  //vector<Type> eta_y2 = ParameterVector("eta_y2",obj,"");
  array<Type> eta_mat = ParameterArray("eta_mat",obj,"");
  vector<Type> eta_y1 = eta_mat.matrix().row(0);
  vector<Type> eta_y2 = eta_mat.matrix().row(1);

  
  nll += SCALE(SEPARABLE(AR1(rho_par),AR1(rho_year)),sd_sel)(eta_mat);

  vector<Type> beta_1(Y);
  vector<Type> beta_2(Y);

  for(int y = 0; y < Y; ++y){
    beta_1(y) = build_beta(lower_b1,upper_b1,logit_m_1,eta_y1(y));
    beta_2(y) = build_beta(lower_b1,upper_b1,logit_m_1,eta_y1(y));
  }

  matrix<Type> S_ly(L,Y);
  for(int y = 0; y < Y; ++y){
    for(int l = 0; l < L; ++l){
      Type l_bin = len_bins(l);
      S_ly(l,y) = make_Sly(l_bin,beta_1(y),beta_2(y));
    }
  }

  Report("sd_sel",sd_sel,obj,"");
  Report("rho_par",rho_par,obj,"");
  Report("rho_year",rho_year,obj,"");
  Report("eta_mat",eta_mat,obj,"");
  Report("beta_1",beta_1,obj,"");
  Report("beta_2",beta_2,obj,"");

  return S_ly;
}

template<class Type>
matrix<Type> convert_Sly_to_Say(matrix<Type> &S_ly,vector<matrix <Type>> &plas,int A){
  int Y = S_ly.cols();
  int L = S_ly.rows();
  matrix<Type> S_ay(A,Y);

  for(int y = 0; y < Y; ++y){
    matrix<Type> pla_y = plas(y);
    for(int a = 0; a < A; ++a){
      S_ay(a,y) = 0;
      for(int l = 0; l < L; ++l){
	S_ay(a,y) += pla_y(l,a)*S_ly(l,y);
      }
    }
  }

  return S_ay;
}

template<class Type>
vector<Type> convert_Sly_to_Say(vector<Type> &S_ly,matrix <Type> &plas,int A){
  int L = S_ly.size();
  vector<Type> S_ay(A);

    for(int a = 0; a < A; ++a){
      S_ay(a) = 0;
      for(int l = 0; l < L; ++l){
	S_ay(a) += plas(l,a)*S_ly(l);
      }
    }
  

  return S_ay;
}


template<class Type>
matrix<Type> rw_sel(int L,int Y,vector<Type> len_bins, Type &nll,objective_function<Type> *obj){
  using namespace density;
  
  matrix<Type> S_ly(L,Y);

  int mora_year = DataInteger("mora_year",obj,"");
  Type phi_age = Parameter("logit_phi_age",obj,"");
  phi_age = invlogit(phi_age);
  Type sd_sel = Parameter("log_sd_sel",obj,"");
  sd_sel = exp(sd_sel);
  matrix<Type> logit_sly = ParameterMatrix("logit_sly",obj,"");
  N01<Type> sel_age;
  AR1_t<N01<Type> > sel_rw(phi_age,sel_age);
  SCALE_t< AR1_t<N01<Type > > > sc_sel_rw(sel_rw,sd_sel);
  
  for(int y = 0; y < Y; ++y){
    vector<Type> sel_col = logit_sly.col(y);
    nll += sc_sel_rw(sel_col);
  }

  for(int y = 0; y < Y; ++y){
    for(int l = 0; l < L; ++l){
      //vector<Type> sel_col = logit_sly.col(y);
      S_ly(l,y) = invlogit(logit_sly(l,y));
    }
  }

  return S_ly;
  
}


template<class Type>
matrix<Type> rwrw_sel(int L,int Y,vector<Type> len_bins,Type &nll,objective_function<Type> *obj){
  matrix<Type> S_ly(L,Y);
  //the year when moratorium started when we'll shift gears
  int mora_year = DataInteger("mora_year",obj,"");

  vector<Type> rw1 = ParameterVector("rw1",obj);
  vector<Type> rw2 = ParameterVector("rw2",obj);
  Type sd_rw1 = Parameter("log_sd_rw1",obj);
  Type sd_rw2 = Parameter("log_sd_rw2",obj);
  sd_rw1 = exp(sd_rw1);
  sd_rw2 = exp(sd_rw2);

  for(int i = 1; i < rw1.size(); ++i){
    nll -= dnorm(rw1(i-1),rw1(i),sd_rw1,true);
    nll -= dnorm(rw2(i-1),rw2(i),sd_rw2,true);
  }

  for(int y = 0; y < mora_year; ++y){
      for(int l = 0; l < L; ++l){
	S_ly(l,y) = invlogit(rw1(l));
      }
  }

  for(int y = mora_year; y < Y; ++y){
      for(int l = 0; l < L; ++l){
	S_ly(l,y) = invlogit(rw2(l));
      }
  }

  return S_ly;
}
 
  

template<class Type>
matrix<Type> blocked_super_simple_sel(int L,int Y,vector<Type> len_bins,Type &nll,objective_function<Type> *obj){
    matrix<Type> S_ly(L,Y);
  //the year when moratorium started when we'll shift gears
  int mora_year = DataInteger("mora_year",obj,"");

  //First block is directed fishery, assuming logistic selectivity
  Type b_beta1 = Parameter("log_b_beta1",obj,"");
  b_beta1 = exp(b_beta1);
  //Beta 2 is now L_95
  //Type b_beta2 = Parameter("log_b_beta2",obj,"");
  //b_beta2 = exp(b_beta2);
  Type rounding_bit = DataScalar("rounding_bit",obj,"");
  Type b_beta2 = b_beta1+rounding_bit*b_beta1;
  //Type r_b_beta2 = -1*log(0.05/0.95)/(b_beta2-b_beta1);
  Type r_b_beta2 = -1*(b_beta2-b_beta1)/(log(0.05/0.95));
  
  for(int y = 0; y < mora_year; ++y){
    for(int l = 0; l < L; ++l){
      Type l_bin = len_bins(l);
      S_ly(l,y) = make_Sly(l_bin,b_beta1,r_b_beta2);
    }
  }
  
  //Assuming Gamma distributed Selectivity
  Type sel_shape = Parameter("log_sel_shape",obj,"");
  sel_shape = exp(sel_shape);
  Type sel_scale = Parameter("log_sel_scale",obj,"");
  sel_scale = exp(sel_scale);

  for(int y = mora_year; y < Y; ++y){
    for(int l = 0; l < L; ++l){
      Type l_bin = len_bins(l);
      S_ly(l,y) = dgamma(l_bin,sel_shape,sel_scale,0);
    }
  }

  //Adjust to be between 0 and 1?
  for(int y = mora_year; y < Y; ++y){
    Type maxy = S_ly.col(y).maxCoeff();
    S_ly.col(y) = S_ly.col(y)/maxy;
  }

  Adreport("sel_shape",sel_shape,obj);
  Adreport("sel_scale",sel_scale,obj);
  Adreport("b_beta1",b_beta1,obj);
  Adreport("b_beta2",b_beta2,obj);
  
  return S_ly;
  
}


  
template<class Type>
void aggregated_F(matrix<Type> &S_ay,matrix<Type> &F_ay,Type &nll,objective_function<Type> *obj,vector<Type> &log_Fy,Type &Fy_sd){

  int Y = S_ay.cols();
  int A = S_ay.rows();

  //Adding bits for projections...
   int proj_type = DataInteger("proj_type",obj);
   int og_Y = DataInteger("og_Y",obj);
   int r_proj = DataInteger("r_proj",obj);
  
 
  for(int y = 1; y < log_Fy.size(); ++y){
    nll -= dnorm(log_Fy(y),log_Fy(y-1),Fy_sd,true);
  }

  vector<Type> Fy = exp(log_Fy);

  for(int a = 0; a < A; ++a){
    for(int y = 0; y < og_Y; ++y){
      F_ay(a,y) = Fy(y)*S_ay(a,y);
    }
  }

  vector<Type> proj_Fy;
  if(proj_type == 0){
    proj_Fy = DataVector("supplied_F",obj);
  }

  for(int a = 0; a < A; ++a){
    for(int y = og_Y; y < Y; ++y){
      if(proj_type == 0){
	int yc = y-(og_Y);
	F_ay(a,y) = proj_Fy(yc)*S_ay(a,y);
      }else{
	F_ay(a,y) = 0;
      }
    }
  }
  
  

  Report("log_Fy",log_Fy,obj,"");
  Adreport("Fy_sd",Fy_sd,obj,"");

}

template<class Type>
void neo_cohort_effects_N(matrix<Type> &N,matrix<Type> &Z, Type &nll,objective_function<Type> *obj,matrix<Type> &log_N_a,Type recruit_sd,matrix<Type> &F,matrix<Type> &S_ay,matrix<Type> &M,newton::newton_config_t<Type> cfg){
  int A = N.rows();
  int Y = N.cols();

  //Adding bits for projections...
   int proj_type = DataInteger("proj_type",obj);
   int og_Y = DataInteger("og_Y",obj);
   int r_proj = DataInteger("r_proj",obj);
  
  Type surv_sd = Parameter("log_surv_sd",obj,"");
  //Type surv_plus_sd = Parameter("log_surv_plus_sd",obj,"");
  //surv_plus_sd = exp(surv_plus_sd);
  surv_sd = exp(surv_sd);

  Type N0_sd = Parameter("log_N0_sd",obj,"");
  N0_sd = exp(N0_sd);

  //The predicted mean
  matrix<Type> mpredN(A,Y);
  //The dynamics one
  matrix<Type> mpredN2(A,Y);

  //The "obsevered"
  matrix<Type> mcomp(A,Y);

  mpredN2(0,0) = log_N_a(0,0);
  mpredN(0,0) = log_N_a(0,0);
  mcomp(0,0) = log_N_a(0,0);

  //recruits
  for(int y = 1; y < Y; ++y){
    mpredN(0,y) = log_N_a(0,y-1);
    mpredN2(0,y) = log_N_a(0,y-1);
    mcomp(0,y) = log_N_a(0,y);
  }

  //initial numbers
  for(int a = 1; a < A; ++a){
    mpredN(a,0) = -Z(a-1,0);
    mpredN2(a,0) = log_N_a(a,0);
    mcomp(a,0) = log_N_a(a,0)-log_N_a(a-1,0);
    if(a == A-1){
      mpredN(a,0) = -log(1-exp(-Z(a-1,0)));
    }
  }

  vector<Type> NN(A);
  vector<Type> MM(A);
  vector<Type> SS(A);
  vector<Type> proj_F(A);

  vector<Type> given_catch;
  vector<Type> proj_Fys;
  vector<Type> solved_Fs;
  vector<vector<Type> > MMs;
  vector<vector<Type> > SSs;
  vector<vector<Type> > NNs;
  vector<vector<Type> > pFs;
  

  if(proj_type == 1){
    given_catch = DataVector("given_catch",obj);
    //proj_Fys = ParameterVector("proj_Fys",obj);
    vector<Type> ssF(given_catch.size());
    vector<vector<Type>> goose(given_catch.size());
    MMs = goose;
    SSs = goose;
    NNs = goose;
    pFs = goose;
    solved_Fs = ssF;
  }
 
  //Cohort dynamics
  for(int y = 1; y < Y; ++y){


    //Update F/Z for catch projections if they're being used
      //Need to make it so we run the functor to get the new F and hence Z before it's used

       Type cur_catch;
      //Are we in a projection year?
      if(y >= og_Y){
	if(proj_type == 1){
	  for(int a = 0; a < A; ++a){
	    MM(a) = M(a,y-1);
	    SS(a) = S_ay(a,y-1);
	    NN(a) = exp(log_N_a(a,y-1));
	  }
	  
	  cur_catch = given_catch(y-og_Y);
	  Functor<TMBad::ad_aug> Fcalc(cur_catch,MM,SS,NN);
	  vector<Type> aFy(1);
	  aFy(0) = log(0.1);
	  vector<Type> nFy = newton::Newton(Fcalc,aFy,cfg);
	  proj_F = exp(nFy(0))*SS;
	  pFs(y-og_Y) = proj_F;
	  solved_Fs(y-og_Y) = nFy(0);
	  //update this years Z?
	  for(int a = 0; a < A; ++a){
	    Z(a,y) = proj_F(a)+M(a,y);
	    F(a,y) = proj_F(a);
	  }
	  SSs(y-og_Y) = SS;
	  MMs(y-og_Y) = MM;
	  NNs(y-og_Y) = NN;

	}
      }
	

    
    for(int a = 1; a < A; ++a){


      //The actual dynamics
      mcomp(a,y) = log_N_a(a,y);

      mpredN2(a,y) = log_N_a(a-1,y-1)-Z(a-1,y-1);
      if(a == A-1){
	mpredN2(a,y) = logspace_add(mpredN2(a,y),log_N_a(a,y-1)-Z(a,y-1));
      }
      mpredN(a,y) = mpredN2(a,y);
      

      
    }
  }

  //Build sigma for the recruits
  vector<Type> log_rec_sds(Y);
  for(int y = 0; y < Y; ++y){
    log_rec_sds(y) = log(recruit_sd);
  }
  vector<Type> dummy_rho(0);
  
  matrix<Type> rec_sigma = sigmaGen(Y,0,log_rec_sds,dummy_rho,0);

  //Build sigma for every other age class
  vector<Type> log_ev_sds(Y);
  log_ev_sds(0) = log(N0_sd);
  for(int y = 1; y < Y; ++y){
    log_ev_sds(y) = log(surv_sd);
  }

  //JUST FOR THE PLUS GROUP
  vector<Type> log_plus_sds(Y);
  log_plus_sds(0) = log(N0_sd);
  for(int y = 1; y < Y; ++y){
    log_plus_sds(y) = log(surv_sd);
  }

 
  matrix<Type> ev_sigma = sigmaGen(Y,0,log_ev_sds,dummy_rho,0);
  matrix<Type> plus_sigma = sigmaGen(Y,0,log_plus_sds,dummy_rho,0);
  matrix<Type> resNraw(A,Y);
  matrix<Type> resN(A,Y);
  matrix<Type> resN2(A-1,Y);

  Eigen::LLT< Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> > lltCovNrec(rec_sigma);
  matrix<Type> LN = lltCovNrec.matrixL();
  matrix<Type> LinvN = LN.inverse();

  Eigen::LLT< Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> > lltCovNEV(ev_sigma);
  matrix<Type> LNEV = lltCovNEV.matrixL();
  matrix<Type> LinvNEV = LN.inverse();

  Eigen::LLT< Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> > lltCovNPLUS(plus_sigma);
  matrix<Type> LNPLUS = lltCovNPLUS.matrixL();
  matrix<Type> LinvNPLUS = LN.inverse();

  
  density::MVNORM_t<Type> recdens(rec_sigma);
  vector<Type> dum_rel(Y);
  for(int y = 0; y < Y;++y){
    dum_rel(y) = mcomp(0,y)-mpredN(0,y);
  }
  nll += recdens(dum_rel);
  resNraw.row(0) = mcomp.row(0)-mpredN.row(0);
  resN.row(0) = LinvN*(vector<Type>(mcomp.row(0)-mpredN.row(0)));
  
  
  density::MVNORM_t<Type> evdens(ev_sigma);
  density::MVNORM_t<Type> plusdens(plus_sigma);

  
  for(int a = 1; a < A; ++a){
    for(int y = 0; y < Y;++y){
      dum_rel(y) = mcomp(a,y)-mpredN(a,y);
    }
      if(a < A-1){
	nll += evdens(dum_rel);
      }else{
	nll += plusdens(dum_rel);
      }
    
    resNraw.row(a) = mcomp.row(a)-mpredN.row(a);
    if(a < A-1){
      resN.row(a) = LinvNEV*(vector<Type>(mcomp.row(a)-mpredN.row(a)));
      resN2.row(a-1) = LinvNEV*(vector<Type>(mcomp.row(a)-mpredN.row(a)));
    }else{
      resN.row(a) = LinvNPLUS*(vector<Type>(mcomp.row(a)-mpredN.row(a)));
      resN2.row(a-1) = LinvNPLUS*(vector<Type>(mcomp.row(a)-mpredN.row(a)));
    }
  }

  for(int a = 0; a < A; ++a){
    for(int y = 0; y < Y; ++y){
      N(a,y) = exp(log_N_a(a,y));
    }
  }

  vector<Type> log_recruit(Y);
  for(int y = 0; y < Y; ++y){
    log_recruit(y) = log_N_a(0,y);
  }

  Report("ssF",solved_Fs,obj,"");
  Report("NNs",NNs,obj,"");
  Report("SSs",SSs,obj,"");
  Report("MMs",MMs,obj,"");
  

  
  Report("mcomp",mcomp,obj,"");
  Report("mpredN",mpredN,obj,"");
  Report("mpredN2",mpredN2,obj,"");
  Report("rec_sigma",rec_sigma,obj,"");
  Report("ev_sigma",ev_sigma,obj,"");
  Adreport("recruit_sd",recruit_sd,obj,"");
  Adreport("surv_sd",surv_sd,obj,"");
  Adreport("N0_sd",N0_sd,obj,"");
  Adreport("log_recruit",log_recruit,obj,"");
  //Adreport("log_N_a",log_N_a,obj,"");
  Report("resN",resN,obj,"");
  Adreport("resN",resN,obj,"");
  //Adreport("resN2",resN2,obj,"");
  Report("resNraw",resNraw,obj,"");
}



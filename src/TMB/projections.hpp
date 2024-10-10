template<class Type>
struct Functor{
  Type cur_catch;
  vector<Type> M;
  vector<Type> S;
  vector<Type> N;

  Functor(const Type &cur_catch, const vector<Type> &M, const vector<Type> &S, const vector<Type> N) : cur_catch(cur_catch), M(M), S(S), N(N) {}

  Type operator()(const vector<Type> &lF){
    Type F = exp(lF(0));
    vector<Type> FS(S.size());
    vector<Type> ecatch(S.size());
    vector<Type> Z(S.size());
    for(int i = 0; i < S.size(); ++i){
      FS(i) = F*S(i);
      Z(i) = FS(i)+M(i);
      ecatch(i) = N(i)*(1-exp(-Z(i)))*(FS(i)/Z(i));
    }

    Type secatch = ecatch.sum();
    Type ret = (secatch-cur_catch)*(secatch-cur_catch);
    return ret;
				   
  }
  
};

// //This structs holds all the crud needed for projections
// template<class Type>
// struct proj_helper_t{

//   schnutevb_t<Type> *gf;
//   vector<Type> len_bins;
//   vector< vector<Type> > ages_combo;
//   matrix<Type> *QLM;

//   Type recruit_sd;
//   Type Fy_sd;
//   Type Fa_sd;
//   vector<Type> log_Fy;
//   vector<Type> log_Fa;

//   matrix<Type> *N;
//   matrix<Type> *M;
//   matrix<Type> *NL;
//   matrix<Type> *SA;
//   matrix<Type> *F;
//   matrix<Type> *SL;

  

// };


// template<class Type>
// struct projections_t{


  
//   objective_function<Type> *obj;

//   int proj_years;
//   vector<Type> supplied_F;
//   vector<Type> given_catch;
  
//   matrix<Type> N_p;
//   matrix<Type> NL_p;
//   matrix<Type> NL_fall_p;
//   matrix<Type> SA_p;
//   matrix<Type> SL_p;
//   matrix<Type> F_p;
//   matrix<Type> M_p;
//   matrix<Type> Z_p;
//   matrix<Type> N_fall_p;
//   vector<matrix <Type> > PLA_p;
//   vector<matrix <Type> > PLA_fall_p;
//   proj_helper_t<Type> proj_stuff;

//   int A;
//   int Y;
//   int L;
  
//   matrix<Type> sur_N;
//   vector<Type> tot_sur_N;
//   newton::newton_config_t<Type> cfg;

  
//   int proj_type;

//   vector<Type> neo_Fy;
//   vector<Type> neo_Fa;
//   vector<Type> neo_l_rec;

//   projections_t(objective_function<Type> *objective,proj_helper_t<Type> proj_stuff,newton::newton_config_t<Type> tcfg);
//   void log_lik(Type &nll,matrix<Type> &N);

//   void report();
  
 
// };


// template<class Type>
// projections_t<Type>:: projections_t(objective_function<Type> *objective,proj_helper_t<Type> tproj_stuff,newton::newton_config_t<Type> tcfg){



  
//    obj = objective;
//    proj_stuff = tproj_stuff;
//   // // recruit_sd = trecruit_sd;
//   // // Fy_sd = tFy_sd;
//   // // log_Fy = tlog_Fy;
//   // // QLM = tQLM;
//    cfg = tcfg;
  
//   proj_type = DataInteger("proj_type",obj);

//   proj_years = DataInteger("proj_years",obj);
//   A = proj_stuff.N ->rows();
//   Y = proj_stuff.N -> cols();
//   L = proj_stuff.NL -> rows();

//   matrix<Type> tN(A,proj_years);
//   matrix<Type> tNL(L,proj_years);
//   matrix<Type> tF(A,proj_years);
//   matrix<Type> tM(A,proj_years);
//   matrix<Type> tSA(A,proj_years);
//   matrix<Type> tSL(L,proj_years);
//   matrix<Type> tZ(A,proj_years);
//   vector<matrix <Type> > tPLA(proj_years);
//   vector<matrix <Type> > tPLA_fall(proj_years);

//   for(int a = 0; a < A; ++a){
//     for(int yy = 0; yy < proj_years; ++yy){
//       tN(a,yy) = 0;
//       tF(a,yy) = 0;
//       tM(a,yy) = 0;
//       tSA(a,yy) = 0;
//       tZ(a,yy) = 0;
//     }
//   }
  
//   for(int y = 0; y < proj_years; ++y){
//     tM.col(y) = proj_stuff.M -> col(Y);
//     tSA.col(y) = proj_stuff.SA -> col(Y);
//     tSL.col(y) = proj_stuff.SL -> col(Y);
//   }

//   M_p = tM;
//   SA_p = tSA;
//   SL_p = tSL;
//   Z_p = tZ;
//   N_p = tN;
//   N_fall_p = tN;
//   NL_p = tNL;
//   NL_fall_p = tNL;

//   //WE HAVE TO CUSTOM RUN THE GF GRRR
//   schnutevb_t<Type> proj_gf;

//   //Fill it up with the stuff we need
//   switch(proj_stuff.gf -> gf_type){
//   case 3:
//     proj_gf.objective = proj_stuff.gf -> objective;
//     proj_gf.M = proj_stuff.gf -> M;
//     proj_gf.prefix = "projections.";
//     proj_gf.gf_type = proj_stuff.gf -> gf_type;
//     proj_gf.ell = proj_stuff.gf -> ell;
//     proj_gf.L = proj_stuff.gf -> L;
//     proj_gf.k = proj_stuff.gf -> k;
//     proj_gf.s = proj_stuff.gf -> s;
//     proj_gf.S = proj_stuff.gf -> S;
//     //We don't break because we want all the above in case 4
//   case 4:
//     {
//       vector<Type> as_pg(proj_years);
//       proj_gf.as_pg = as_pg;
//       proj_gf.pg_ext = proj_stuff.gf -> pg_ext;
//     break;
//     }
//   default:
//     error("This GF type is not supported for projections sorry!");
//   };

//   // Rcout << "gf_tpype?" << proj_gf.gf_type << std::endl;

//   neo_l_rec = ParameterVector("neo_l_rec",objective);
  

//   if(proj_type == 0){
//       supplied_F = DataVector("supplied_F",obj);
//   }else{
//     given_catch = DataVector("given_catch",obj);
//     neo_Fy = ParameterVector("neo_Fy",obj);
//   }
      

  

  
  // for(int y = 0; y < proj_years; ++y){


  //   Rcout << "y: " << y << std::endl;
    
  //   vector<Type> SS(A);
  //   vector<Type> MM(A);
  //   vector<Type> NN(A);
  //   vector<Type> FF(A);
  //   vector<Type> SSL(L);
    
  //   if(y == 0){
  //     SS = proj_stuff.SA -> col(Y);
  //     MM = proj_stuff.M -> col(Y);
  //     NN = proj_stuff.N -> col(Y);
  //     SSL = proj_stuff.SL -> col(Y); 
  //   }else{
  //     //SS = SA_p.col(y-1);
  //     MM = M_p.col(y-1);
  //     NN = N_p.col(y-1);
  //     SSL = SL_p.col(y-1);
  //   }

    

  //   Type Fy;
  //   Type cur_catch;
    
  //   if(proj_type == 0){
  //     Fy = supplied_F(y);
  //     FF = Fy*SS;
  //   }else{
  //     cur_catch = given_catch(y);
  //     Functor<TMBad::ad_aug> Fcalc(cur_catch,MM,SS,NN);
  //     vector<Type> nFy(1);
  //     nFy(0) = neo_Fy(y);
  //     vector<Type> F_sol = newton::Newton(Fcalc,nFy,cfg);
  //     //FF = neo_Fy(y)*SS;
  //   }

  //   vector<Type> ZZ = FF+MM;
  //   Z_p.col(y) = ZZ;

  //   for(int a = 1; a < A; ++a){
  //     N_p(a,y) = NN(a-1)*exp(-ZZ(a-1));
  //     if(a == A-1){
  // 	N_p(a,y) += NN(a)*exp(-ZZ(a));
  //     }
  //   }

  //   Rcout << "We got here threepo??" << std::endl;
 
  //   Type a_pg;
  //   if(proj_stuff.gf -> gf_type == 4){
  //     if(y == 0){
  // 	a_pg = proj_stuff.gf -> as_pg(Y);
  //     }else{
  // 	a_pg = proj_gf.as_pg(y-1);
  //     }
  //     proj_gf.as_pg(y) = (NN(A-1)*1+NN(A)*(a_pg+1))/(NN(A-1)+NN(A));
  //   }
  //   else{
  //     a_pg = 0;
  //   }
  

  //   Rcout << "We got here fourry??" << std::endl;


  //   Type dum = proj_gf.mean_len_at_age(5,y);
  //   Rcout << "test: " << dum << std::endl;

  //   Type dum2 = proj_gf.sd_len_at_age(5,y);
  //   Rcout << "test2: " << dum2 << std::endl;

  //   Type dum3 = proj_gf.mean_len_at_age(proj_stuff.ages_combo(0)(1),y);
  //   Rcout << "test3: " << dum3 << std::endl;

  //   matrix<Type> PLA_p = prob_len_at_age(proj_gf,proj_stuff.len_bins,proj_stuff.ages_combo(0),y);
  //   matrix<Type> PLA_fall_p = prob_len_at_age(proj_gf,proj_stuff.len_bins,proj_stuff.ages_combo(2),y);
    
  //   //PLA_p(y) = prob_len_at_age(proj_gf,proj_stuff.len_bins,proj_stuff.ages_combo(0),y);
  //    //PLA_fall_p(y) = prob_len_at_age(proj_gf,proj_stuff.len_bins,proj_stuff.ages_combo(2),y);

  //   Rcout << "We got here five??" << std::endl;
    
    
  //    SL_p.col(y) = SSL;
  //    SA_p.col(y) = convert_Sly_to_Say(SSL,PLA_p,A);
  //    F_p.col(y) = FF;

  //    Rcout << "We got here six??" << std::endl;
    
    
    
  //    Type fall_adjust = 10.5/12;

  //    for(int a = 0; a < A; ++a){
  //      N_fall_p(a,y) = N_p(a,y)*exp(-Z_p(a,y)*fall_adjust);
  //    }

  //    Rcout << "We are here?" << std::endl;

  //    NL_p.col(y) = PLA_p*N_p.col(y);
  //    NL_fall_p.col(y) = PLA_fall_p*N_fall_p.col(y);

  //    Rcout << "We made it here?" << std::endl;

  //  }
  
// }

// template<class Type>
// void projections_t<Type>::log_lik(Type &nll,matrix<Type> &N){

  
//   for(int y = 0; y < proj_years; ++y){
//     if(y == 0){
//       vector<Type> Nn = N.col(Y);
//       Rcout << "WHAT: " << Nn << std::endl;
//       // nll -= dnorm(neo_l_rec(y),log(Nn(0)),proj_stuff.recruit_sd,true);
//     }else{
//       //nll -= dnorm(neo_l_rec(y),neo_l_rec(y-1),proj_stuff.recruit_sd,true);
//     }
//     //N_p(0,y) = exp(neo_l_rec(y));
//   }

// }

// template<class Type>
// void projections_t<Type>::report(){
//   if(proj_years > 0){
//     Report("N_p",N_p,obj);
//     Report("F_p",F_p,obj);
//     Report("Z_p",Z_p,obj);
//     Report("N_fall_p",N_fall_p,obj);
    
//     Adreport("N_p",N_p,obj);
//     Adreport("F_p",F_p,obj);
//     Adreport("Z_p",Z_p,obj);
//     Adreport("N_fall_p",N_fall_p,obj);

//   }
// }

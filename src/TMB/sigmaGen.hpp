template<class Type>
matrix<Type> sigmaGen(int nage,int corrFlag, vector<Type> &logsdF,vector<Type> &tRhoF,int ageFree){

  enum corrType{
		independent = 0,
		parallel = 1,
		compound = 2,
		ARone = 3,
		custom = 4,
		ARsome = 5
  };

  matrix<Type> sigma(nage,nage);
  vector<Type> sdF = exp(logsdF);
  vector<Type> rhoF(tRhoF.size());
  for(int i = 0; i < tRhoF.size(); ++i){
    rhoF(i) = rho_trans(tRhoF(i));
  }
  sigma.setZero();

  switch(corrFlag){
  case independent:
    sigma.diagonal() = sdF*sdF;
    break;
  case parallel:
    sigma.diagonal() = sdF*sdF;
    for(int i=0; i < nage; i++){
      for(int j=0; j < i; j++){
	//Anders explained this is needed to get parallel to coverge. VERY PICKY ON NUMBER OF 0.999
	sigma(i,j) = 0.999*sdF(i)*sdF(j);
	sigma(j,i) = sigma(i,j);
      }
    }
    break;
  case compound:
    sigma.diagonal() = sdF*sdF;
    for(int i = 0;i < nage; i++){
      for(int j = 0; j < i; j++){
	sigma(i,j) = rhoF(0)*sdF(i)*sdF(j);
	sigma(j,i) = sigma(i,j);
      }
    }
    break;
  case ARone:
    sigma.diagonal() = sdF*sdF;
    for(int i = 0;i < nage; i++){
      for(int j = 0; j < i; j++){
	sigma(i,j) = pow(rhoF(0),Type(i-j))*sdF(i)*sdF(j);
	sigma(j,i) = sigma(i,j);
      }
    }
    break;
  case custom:
    sigma.diagonal() = sdF*sdF;
    for(int i = 0;i < nage; i++){
      for(int j = 0; j < i; j++){
	sigma(i,j) = rhoF(i)*rhoF(j)*sdF(i)*sdF(j);
	sigma(j,i) = sigma(i,j);
      }
    }
    break;
  case ARsome:
    sigma.diagonal() = sdF*sdF;
    for(int i = 0; i < nage; i++){
      for(int j = 0; j < i; j++){
	if(i <= ageFree || j <= ageFree){
	  sigma(i,j) = 0;
	}else{
	  sigma(i,j) = pow(rhoF(0),Type(i-j))*sdF(i)*sdF(j);
	}
	sigma(j,i) = sigma(i,j);
      }
    }
    break;
  default:
    error("F correlation mode not supported");
    break;
  }

  return sigma;
}


template<class Type>
matrix<Type> sigmaGen(int nage,vector<Type> logsdF,vector<Type> tRhoF){
  return sigmaGen(nage,3,logsdF,tRhoF,0);
}
  

template<class Type>
matrix<Type> sigmaGen(int nage,vector<Type> logsdF,Type tRhoF){
  vector<Type> tRhoFv(1);
  tRhoFv(0) = tRhoF;
  return sigmaGen(nage,3,logsdF,tRhoFv,0);
}

template<class Type>
matrix<Type> sigmaGen(int nage,int corrFlag,vector<Type> logsdF,Type tRhoF){
  vector<Type> tRhoFv(1);
  tRhoFv(0) = tRhoF;
  return sigmaGen(nage,corrFlag,logsdF,tRhoFv,0);
}

template<class Type>
matrix<Type> sigmaGen(int nage,int corrFlag,vector<Type> logsdF,Type tRhoF,int ageFree){
  vector<Type> tRhoFv(1);
  tRhoFv(0) = tRhoF;
  return sigmaGen(nage,corrFlag,logsdF,tRhoFv,ageFree);
}

template<class Type>
matrix<Type> sigmaGen(int nage,vector<Type> log_sd){
  vector<Type> tRhoFv(0);
  return sigmaGen(nage,0,log_sd,tRhoFv,0);
}

template<class Type>
matrix<Type> sigmaGen(int nage,Type log_sd){
  vector<Type> tRhoFv(0);
  vector<Type> log_sdv(nage);
  for(int i = 0; i < nage; ++i){
    log_sdv(i) = log_sd;
  }
  return sigmaGen(nage,0,log_sdv,tRhoFv,0);
}

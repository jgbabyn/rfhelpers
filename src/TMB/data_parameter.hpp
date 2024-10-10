template<class Type>
inline Type Parameter(const char *nam,objective_function<Type> *obj){
  Type ret(obj -> fillShape(asVector<Type>(obj -> getShape(nam,&isNumericScalar)),nam)[0]);
  return ret;
}

template<class Type>
inline Type Parameter(std::string name,objective_function<Type> *obj,std::string prefix){
  std::string full = prefix + name;
  char * nam = new char[full.size()+1];
  std::strcpy(nam,full.c_str());
  Type ret(obj -> fillShape(asVector<Type>(obj -> getShape(nam,&isNumericScalar)),nam)[0]);
  return ret;
}

template<class Type>
inline vector<Type> ParameterVector(const char *nam,objective_function<Type> *obj){
  vector<Type> ret(obj ->fillShape(asVector<Type>(obj -> getShape(nam, &Rf_isNumeric)),nam));
  return ret;
}

template<class Type>
inline vector<Type> ParameterVector(std::string name,objective_function<Type> *obj,std::string prefix){
  std::string full = prefix + name;
  char * nam = new char[full.size()+1];
  std::strcpy(nam,full.c_str());
  vector<Type> ret(obj ->fillShape(asVector<Type>(obj -> getShape(nam, &Rf_isNumeric)),nam));
  return ret;
}

template<class Type>
inline matrix<Type> ParameterMatrix(const char *nam,objective_function<Type> *obj){
  tmbutils::matrix<Type> ret(obj -> fillShape(asMatrix<Type> ( obj -> getShape(nam, &Rf_isMatrix) ),nam));
  return ret;
}

template<class Type>
inline matrix<Type> ParameterMatrix(std::string name,objective_function<Type> *obj,std::string prefix){
  std::string full = prefix + name;
  char * nam = new char[full.size()+1];
  std::strcpy(nam,full.c_str());
  tmbutils::matrix<Type> ret(obj -> fillShape(asMatrix<Type> ( obj -> getShape(nam, &Rf_isMatrix) ),nam));
  return ret;
}

template<class Type>
inline tmbutils::array<Type> ParameterArray(const char *nam,objective_function<Type> *obj){
  tmbutils::array<Type> ret(obj -> fillShape(tmbutils::asArray<Type>(obj -> getShape(nam, &Rf_isArray)),nam));
  return ret;
}

template<class Type>
inline tmbutils::array<Type> ParameterArray(std::string name,objective_function<Type> *obj,std::string prefix){
  std::string full = prefix + name;
  char * nam = new char[full.size()+1];
  std::strcpy(nam,full.c_str());
  tmbutils::array<Type> ret(obj -> fillShape(tmbutils::asArray<Type>(obj -> getShape(nam, &Rf_isArray)),nam));
  return ret;
}



template<class Type>
inline vector<Type> DataVector(const char *nam, objective_function<Type> *obj){
  vector<Type> ret = asVector<Type>(getListElement(obj -> data,nam,&Rf_isNumeric));
  return ret;
}

template<class Type>
inline vector<Type> DataVector(std::string name, objective_function<Type> *obj,std::string prefix){
  std::string full = prefix + name;
  char * nam = new char[full.size()+1];
  std::strcpy(nam,full.c_str());
  vector<Type> ret = asVector<Type>(getListElement(obj -> data,nam,&Rf_isNumeric));
  return ret;
}


template<class Type>
inline matrix<Type> DataMatrix(const char *nam, objective_function<Type> *obj){
  matrix<Type> ret(asMatrix<Type>(getListElement(obj -> data,nam, &Rf_isMatrix)));
  return ret;
}

template<class Type>
inline matrix<Type> DataMatrix(std::string name, objective_function<Type> *obj,std::string prefix){
  std::string full = prefix + name;
  char * nam = new char[full.size()+1];
  std::strcpy(nam,full.c_str());
  matrix<Type> ret(asMatrix<Type>(getListElement(obj -> data,nam, &Rf_isMatrix)));
  return ret;
}

//Do not use
// template<class Type>
// inline matrix<Type> glanceParameterMatrix(std::string name, objective_function<Type> *obj,std::string prefix){
//   std::string full = prefix + name;
//   char * nam = new char[full.size()+1];
//   std::strcpy(nam,full.c_str());
//   matrix<Type> ret(asMatrix<Type>(getListElement(obj -> parameters,nam, &Rf_isMatrix)));
//   return ret;
// }

template<class Type>
inline Type DataScalar(const char *nam, objective_function<Type> *obj){
  Type ret(asVector<Type>(getListElement(obj ->data,nam,&isNumericScalar))[0]);
  return ret;
}

template<class Type>
inline Type DataScalar(std::string name, objective_function<Type> *obj,std::string prefix){
  std::string full = prefix + name;
  char * nam = new char[full.size()+1];
  std::strcpy(nam,full.c_str());
  Type ret(asVector<Type>(getListElement(obj ->data,nam,&isNumericScalar))[0]);
  return ret;
}

template<class Type>
inline int DataInteger(const char *nam, objective_function<Type> *obj){
  int ret(CppAD::Integer(asVector<Type>(getListElement(obj -> data,nam, &isNumericScalar))[0]));
  return ret;
}


template<class Type>
inline int DataInteger(std::string name, objective_function<Type> *obj,std::string prefix){
  std::string full = prefix + name;
  char * nam = new char[full.size()+1];
  std::strcpy(nam,full.c_str());
  int ret(CppAD::Integer(asVector<Type>(getListElement(obj -> data,nam, &isNumericScalar))[0]));
  return ret;
}


template<class Type>
inline vector<int> DataFactor(const char *nam, objective_function<Type> *obj){
  vector<int> ret(asVector<int>(getListElement(obj -> data, nam, &Rf_isNumeric)));
  return ret;
}

template<class Type>
inline vector<int> DataFactor(std::string name, objective_function<Type> *obj,std::string prefix){
  std::string full = prefix + name;
  char * nam = new char[full.size()+1];
  std::strcpy(nam,full.c_str());
  vector<int> ret(asVector<int>(getListElement(obj -> data, nam, &Rf_isNumeric)));
  return ret;
}

template<class Type>
inline vector<int> DataIVector(const char *nam, objective_function<Type> *obj){
  vector<int> ret(asVector<int>(getListElement(obj -> data, nam, &Rf_isNumeric)));
  return ret;
}

template<class Type>
inline vector<int> DataIVector(std::string name, objective_function<Type> *obj,std::string prefix){
  std::string full = prefix + name;
  char * nam = new char[full.size()+1];
  std::strcpy(nam,full.c_str());
  vector<int> ret(asVector<int>(getListElement(obj -> data, nam, &Rf_isNumeric)));
  return ret;
}

template<class Type>
inline vector<std::string> DataStrings(std::string name, objective_function<Type> *obj,std::string prefix){
  std::string full = prefix + name;
  char * nam = new char[full.size()+1];
  std::strcpy(nam,full.c_str());
  SEXP slist = getListElement(obj -> data,nam);
  vector<std::string> ret(LENGTH(slist));
  for(int i = 0; i < LENGTH(slist); ++i){
    SEXP v_element = VECTOR_ELT(slist,i);
    ret[i] = CHAR(STRING_ELT(v_element,0));
  }

  return ret;  
}

template<class Type>
inline vector<matrix<Type> > DataMatrices(std::string name, objective_function<Type> *obj, std::string prefix){
  std::string full = prefix + name;
  char * nam = new char[full.size()+1];
  std::strcpy(nam,full.c_str());
  SEXP mlist = getListElement(obj -> data,nam);
  vector< matrix<Type> > ret(LENGTH(mlist));
  for(int i = 0; i < LENGTH(mlist); ++i){
    SEXP v_element = VECTOR_ELT(mlist,i);
    ret[i] = asMatrix<Type>(v_element);
  }
  return ret;
}

template<class Type>
inline Eigen::SparseMatrix<Type> DataSparseMatrix(const char *nam, objective_function<Type> *obj){
  Eigen::SparseMatrix<Type> ret(tmbutils::asSparseMatrix<Type>(getListElement(obj -> data,nam, &isValidSparseMatrix)));
  return ret;
}

template<class Type>
inline Eigen::SparseMatrix<Type> DataSparseMatrix(std::string name, objective_function<Type> *obj,std::string prefix){
  std::string full = prefix + name;
  char * nam = new char[full.size()+1];
  std::strcpy(nam,full.c_str());
  Eigen::SparseMatrix<Type> ret(tmbutils::asSparseMatrix<Type>(getListElement(obj -> data,nam, &isValidSparseMatrix)));
  return ret;
}


template<class Type>
inline array<Type> DataArray(const char *nam, objective_function<Type> *obj){
  tmbutils::array<Type> ret;                                             
  if (!Rf_isNull(getListElement(obj -> parameters,nam))){ 
    ret = obj -> fillShape(tmbutils::asArray<Type>(obj -> getShape(nam, &Rf_isArray)), nam);    
  } else {                                                                
    ret = tmbutils::asArray<Type>(getListElement(obj -> data, nam, &Rf_isArray));               
  }
  return ret;
}

template<class Type>
inline array<Type> DataArray(std::string name, objective_function<Type> *obj,std::string prefix){
  std::string full = prefix + name;
  char * nam = new char[full.size()+1];
  std::strcpy(nam,full.c_str());
  tmbutils::array<Type> ret;                                             
  if (!Rf_isNull(getListElement(obj -> parameters,nam))){ 
    ret = obj -> fillShape(tmbutils::asArray<Type>(        
						   obj -> getShape(nam, &Rf_isArray)), nam);    
  } else {                                                                
    ret = tmbutils::asArray<Type>(getListElement(                        
						 obj -> data, nam, &Rf_isArray));               
  }
  return ret;

}

template<class Type, class T>
inline void Report(const char *nam,T item, objective_function<Type> *obj){
  if( isDouble<Type>::value &&                                            
      obj -> current_parallel_region<0 )                    
    {                                                                       
      SEXP _TMB_temporary_sexp_;                                          
      PROTECT( _TMB_temporary_sexp_ = asSEXP(item) );                     
      Rf_defineVar(Rf_install(nam),                                     
		   _TMB_temporary_sexp_, obj -> report);    
      UNPROTECT(1);                                                       
    }
}

template<class Type, class T>
inline void Report(std::string name,T item, objective_function<Type> *obj,std::string prefix){
  std::string full = prefix + name;
  char * nam = new char[full.size()+1];
  std::strcpy(nam,full.c_str());
  if( isDouble<Type>::value &&                                            
      obj -> current_parallel_region<0 )                    
    {                                                                       
      SEXP _TMB_temporary_sexp_;                                          
      PROTECT( _TMB_temporary_sexp_ = asSEXP(item) );                     
      Rf_defineVar(Rf_install(nam),                                     
		   _TMB_temporary_sexp_, obj -> report);    
      UNPROTECT(1);                                                       
    }
}

template<class Type, class T>
inline void Adreport(const char *nam, T item, objective_function<Type> *obj){
  obj -> reportvector.push(item, nam);
}

template<class Type, class T>
inline void Adreport(std::string name, T item, objective_function<Type> *obj,std::string prefix){
  std::string full = prefix + name;
  char * nam = new char[full.size()+1];
  std::strcpy(nam,full.c_str());
  obj -> reportvector.push(item, nam);
}

template<class Type, class T>
inline void Adreport(std::string name, T item, objective_function<Type> *obj,std::string prefix,bool do_adreport){
  if(do_adreport == true){
    std::string full = prefix + name;
    char * nam = new char[full.size()+1];
    std::strcpy(nam,full.c_str());
    obj -> reportvector.push(item, nam);
  }
}

//Kind of replicates SIMULATE macro
template<class Type>
inline bool check_simulate(objective_function<Type> *obj){
  return isDouble<Type>::value && obj -> do_simulate;
}

template<class Type, class T>
inline void Report(const char *nam,vector<vector<T> > item,objective_function<Type> *obj){
  if(isDouble<Type>::value && obj -> current_parallel_region < 0){
      SEXP _TMB_temporary_sexp_ = PROTECT(allocVector(VECSXP,item.size()));
      for(int i = 0; i < item.size(); ++i){
	SEXP temp = asSEXP(item(i));
	SET_VECTOR_ELT(_TMB_temporary_sexp_,i,temp);
      }
      Rf_defineVar(Rf_install(nam),_TMB_temporary_sexp_,
		   obj -> report);
      UNPROTECT(1);
    }

}

template<class Type, class T>
inline void Report(const char *nam,vector<matrix<T> > item,objective_function<Type> *obj){
  if(isDouble<Type>::value && obj -> current_parallel_region < 0){
      SEXP _TMB_temporary_sexp_ = PROTECT(allocVector(VECSXP,item.size()));
      for(int i = 0; i < item.size(); ++i){
	SEXP temp = asSEXP(item(i));
	SET_VECTOR_ELT(_TMB_temporary_sexp_,i,temp);
      }
      Rf_defineVar(Rf_install(nam),_TMB_temporary_sexp_,
		   obj -> report);
      UNPROTECT(1);
    }

}

template<class Type, class T>
inline void AdreportList(std::string nam,vector<vector<T> > item, objective_function<Type> *obj){
  for(int i = 0; i < item.size(); ++i){
    int i_adj = i+1;
    std::string istring = std::to_string(i_adj);
    std::string name_temp = nam+ "_" + istring;
    char* nammy = new char[name_temp.size()+1];
    std::strcpy(nammy,name_temp.c_str());
    vector<T> element = item[i];
    obj -> reportvector.push(element,nammy);
  }
}

template<class Type, class T>
inline void AdreportList(std::string nam,vector<vector<T> > item, objective_function<Type> *obj,bool do_adreport){
  if(do_adreport == true){
    for(int i = 0; i < item.size(); ++i){
      int i_adj = i+1;
      std::string istring = std::to_string(i_adj);
      std::string name_temp = nam+ "_" + istring;
      char* nammy = new char[name_temp.size()+1];
      std::strcpy(nammy,name_temp.c_str());
      vector<T> element = item[i];
      obj -> reportvector.push(element,nammy);
    }
  }
}

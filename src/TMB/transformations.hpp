// template<class Type>
// Type logit(Type p){
//   return log(p/(1-p));
// }

// template<class Type>
// vector<Type> logit(vector<Type> p){
//   s = p.size();
//   vector<Type> ret(s);
//   for(int i = 0; i < s; ++i){
//     ret(i) = logit(p(i));
//   }
//   return ret;
// }
    

//Assume we have a vector where we want x_k < x_{k+1}
template<class Type>
vector<Type> ordered_inv_transform(vector<Type> y){
  vector<Type> x(y.size());
  x(0) = y(0);
  for(int i = 1; i < x.size(); ++i){
    x(i) = x(i-1)+exp(y(i));
  }
  return x;
}

//Logistic but with a different max
template<class Type>
Type glogit(Type x, Type L){
  return L/(1+exp(-x));
}

template<class Type>
Type ub_inv(Type y,Type b){
  Type x = b -exp(y);
  return x;
}

template<class Type>
Type iglogit(Type y, Type L){
  return -log(L/y-1);
}

//Simplex transform with a generic max
template<class Type>
vector<Type> simplex_transform(vector<Type> x,Type L){
  vector<Type> zk(x.size()-1);
  vector<Type> y(x.size()-1);
  zk(0) = x(0);
  Type K = x.size();

  Type sumx = 0;
  for(int i = 1; i < zk.size(); ++i){
    sumx += x(i-1);
    zk(i) = x(i)/(L-sumx);
  }

  Type sumk = 0;
  for(int i = 0; i < y.size(); ++i){
    sumk += 1;
    y(i) = glogit(zk(i),L)-log(L/(K-sumk));
  }

  return y;
}

template<class Type>
vector<Type> inv_simplex_transform(vector<Type> y, Type L){
  vector<Type> zk(y.size());
  vector<Type> x(y.size()+1);
  Type K = x.size();
  Type sumk = 0;
  for(int i = 0; i < zk.size(); ++i){
    sumk += 1;
    Type what = y(i)+log(L/(K-sumk));
    zk(i) = iglogit(y(i)+log(L/(K-sumk)),L);
  }
  x(0) = zk(0);
  
 
  Type sumx = 0;
  for(int i = 1; i < zk.size(); ++i){
    sumx += x(i-1);
    x(i) = (L-sumx)*zk(i);
  }
  sumx += x(y.size()-1);
  x(y.size()) = L -sumx;
  return x;
}

template<class Type>
vector<Type> inv_soft_max(vector<Type> y){
  vector<Type> z = log(y);
  vector<Type> x = z-z.mean();
  return x;
}

//' @name rho_trans
//' @description Convert from optimizer space to -1 to 1 for correlation parameters
//' @param x the variable to convert
template <class Type>
Type rho_trans(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}

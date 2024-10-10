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


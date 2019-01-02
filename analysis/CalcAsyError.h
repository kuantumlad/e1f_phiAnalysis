#ifndef calcasyerror_h
#define calcasyerror_h


double CalcAsyError(std::vector<double> v_err ){

  double nom_asy = v_err[0];
  double sum = 0.0;
  for( int i = 0; i < v_err.size(); i++ ){
    sum = sum + ( v_err[i]*v_err[i] - nom_asy*nom_asy );    
  }
  
  return sum/(v_err.size() - 1); 
}

#endif

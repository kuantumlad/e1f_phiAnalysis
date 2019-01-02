#ifndef calcbsaerr_h
#define calcbsaerr_h

double CalcBSAErr( int numP, int numN, double avg_pol ){
  //double bsa_err = sqrt( 4.0/(avg_pol*avg_pol) * (1.0/ std::pow(numP + numN, 4)) * ( std::pow(numP,2)*numN + std::pow(numN,2)*numP )) ;
  double temp_bsa = (1.0/avg_pol) * (( numP - numN )/(numP + numN));
  double bsa_err = sqrt(( (1.0 - temp_bsa*temp_bsa) / ( numP + numN ) ));
  std::cout << " error: " << bsa_err << std::endl;
  return bsa_err;
  
}
#endif

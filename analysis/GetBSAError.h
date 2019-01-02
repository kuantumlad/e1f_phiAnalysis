#ifndef getbsaerror_h
#define getbsaerror_h

std::vector<double> GetBSAError( int n_bins, TH1D *h_P, TH1D *h_N, double avg_pol ){

  std::vector<double> v_bsa_err;
  for(int i = 0; i < n_bins; i++ ){
    int numP = h_P->GetBinContent(i+1);   
    int numN = h_N->GetBinContent(i+1);  

    double bsa_err = sqrt( 4.0/(avg_pol*avg_pol) * (1.0/ std::pow(numP + numN, 4)) * ( std::pow(numP,2)*numN + std::pow(numN,2)*numP )) ;
    std::cout << " bin " << i << " error: " << bsa_err << std::endl;
    v_bsa_err.push_back(bsa_err);  
  }
  return v_bsa_err;
}
#endif

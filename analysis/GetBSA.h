#ifndef GetBSA_h
#define GetBSA_h

#include "GetBSAError.h"
#include "CalcBSAErr.h"

TGraphErrors* GetBSA(int n_bins, TH1D *h_asy_phi_P, TH1D *h_asy_phi_N, double avg_pol ){
   
  //TH1D *h_asy_phi_top = new TH1D("h_asy_phi_top","h_asy_phi_top",n_bins, -180.0, 180.0);// = h_asy_phi_P;
  //TH1D *h_asy_phi_bot = new TH1D("h_asy_phi_bot","h_asy_phi_bot",n_bins, -180.0, 180.0);;// = h_asy_phi_P;
  //TH1D *h_asy_phi_final =  new TH1D("h_asy_phi_final","h_asy_phi_final",n_bins,-180.0, 180.0);
   
   //h_asy_phi_top->Add(h_asy_phi_P, h_asy_phi_N, 1.0, -1.0);
   //h_asy_phi_bot->Add(h_asy_phi_P, h_asy_phi_N, 1.0, 1.0);

   std::vector<double> asy_values_to_graph;
   std::vector<double> v_bsa_err;
   std::vector<double> v_bin_err;
   std::vector<double> asy_bin_center;  
   
   for( int k = 0; k < n_bins; k++ ){

     int numP = h_asy_phi_P->GetBinContent(k+1);
     int numN = h_asy_phi_N->GetBinContent(k+1);
     double bin_center = h_asy_phi_P->GetBinCenter(k+1);
     
     double bsa = ((double)numP - (double)numN) / ( (double)numP + (double)numN );
     std::cout << "P " << numP << " N " << numN << " bsa calc " << bsa << std::endl;
     double bsa_err = CalcBSAErr( numP, numN, avg_pol );

     //h_asy_phi_top->SetBinContent(k, numP);
     //h_asy_phi_bot->SetBinContent(k, numN);     
     //h_asy_phi_final->SetBinContent(k, bsa );

     asy_values_to_graph.push_back( bsa );
     v_bsa_err.push_back(bsa_err);
     v_bin_err.push_back(0.0);
     asy_bin_center.push_back(bin_center);
         
   }

   
   //h_asy_phi_final->Divide(h_asy_phi_top, h_asy_phi_bot, 1.0, 0.75);

   //std::vector<double> bsa_err;
   //bsa_err = GetBSAError( n_bins, h_asy_phi_P, h_asy_phi_N, avg_pol );
   
   TGraphErrors *g_temp = new TGraphErrors(asy_bin_center.size(),  &(asy_bin_center[0]), &(asy_values_to_graph[0]), &(v_bin_err[0]), &(v_bsa_err[0]));
   //for( int i = 0; i < n_bins; i++ ){
   //  g_temp->SetPointError( i, 0.0, bsa_err[i] );    
   //}  
   return g_temp;
}
#endif

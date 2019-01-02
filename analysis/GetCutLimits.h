#ifndef getcutlimits_h
#define getcutlimits_h

void GetCutLimits( std::map< std::string, std::vector< std::vector< double> > >& m_cutlim ){

  std::vector< double > cc_1, cc_2, cc_0;// 1 == loose, 2 == tight
  cc_0 = {-1.0, 1.0};
  cc_1 = {-1.1, 1.1 };
  cc_2 = {-1.0, 1.0 };  
  m_cutlim["dist_cc"] = {cc_0, cc_1, cc_2};

  std::vector< double > dcr1_1, dcr1_2, dcr1_0;// 1 == loose, 2 == tight
  dcr1_0 = {-1.0, 1.0};
  dcr1_1 = {-1.1, 1.0 };
  dcr1_2 = {-0.9, 1.0 };
  m_cutlim["dist_dcr1"] = { dcr1_0, dcr1_1, dcr1_2 };

  std::vector< double > dcr3_1, dcr3_2, dcr3_0;// 1 == loose, 2 == tight
  dcr3_0 = {-1.0, 1.0};
  dcr3_1 = {-1.1, 1.0 };
  dcr3_2 = {-0.95, 1.0 };
  m_cutlim["dist_dcr3"] = { dcr3_0, dcr3_1, dcr3_2 };

  std::vector< double > ecsf_1, ecsf_2, ecsf_0;// 1 == loose, 2 == tight
  ecsf_0 = {-1.0, 1.0};
  ecsf_1 = {-1.1, 1.1 };
  ecsf_2 = {-0.9, 0.9 };
  m_cutlim["dist_ecsf"] =  { ecsf_0, ecsf_1, ecsf_2 };

  std::vector< double > edep_1, edep_2, edep_0;// 1 == loose, 2 == tight
  edep_0 = {-1.0, 1.0};
  edep_1 = {-1.1, 1.0 };
  edep_2 = {-0.9, 1.0 };
  m_cutlim["dist_edep"] = { edep_0, edep_1, edep_2 };

  std::vector< double > ecu_1, ecu_2, ecu_0;// 1 == loose, 2 == tight
  ecu_0 = {-1.0,1.0};
  ecu_1 = {-1.1, 1.1 };
  ecu_2 = {-0.9, 0.9 };
  m_cutlim["dist_ecu"] = { ecu_0, ecu_1, ecu_2 };

  std::vector< double > ecv_1, ecv_2, ecv_0;// 1 == loose, 2 == tight
  ecv_0 = {-1.0, 1.0};
  ecv_1 = {-1.0, 1.1 };
  ecv_2 = {-1.0, 0.98 };
  m_cutlim["dist_ecv"] = { ecv_0, ecv_1, ecv_2 };

  std::vector< double > ecw_1, ecw_2, ecw_0;// 1 == loose, 2 == tight
  ecw_0 = {-1.0, 1.0};
  ecw_1 = {-1.0, 1.1 };
  ecw_2 = {-1.0, 0.98 };
  m_cutlim["dist_ecw"] = { ecw_0, ecw_1, ecw_2 };

  std::vector< double > vz_1, vz_2, vz_0;// 1 == loose, 2 == tight
  vz_0 = {-1.0, 1.0};
  vz_1 = {-1.1, 1.1 };
  vz_2 = {-0.95, 0.95 };
  m_cutlim["dist_vz"] = { vz_0, vz_1, vz_2 };
    

}
#endif

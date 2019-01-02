#ifndef physicshistofiller_hh
#define physicshistofiller_hh

class PhysicsHistoFiller{

 public:
  PhysicsHistoFiller();
  PhysicsHistoFiller(const char*);
  virtual ~PhysicsHistoFiller();


 public:
  const char* temp_hname;
  PhysicsHistoPlotter h_epX;
  PhysicsHistoPlotter h_epkpX;
  PhysicsHistoPlotter h_ekpX;

 public:
  void InitAllHistos();

  void FillPhysics(std::vector<TLorentzVector> );
  //void FillPhysics(TLorentzVector);
  void FillPhysicsH1(TLorentzVector);
  void FillPhysicsH2(TLorentzVector);
  void PrintAllHisto();
  


};
#endif

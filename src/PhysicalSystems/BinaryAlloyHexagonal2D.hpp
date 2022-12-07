#ifndef BIN_ALLOY_HEXAGONAL2D_HPP
#define BIN_ALLOY_HEXAGONAL2D_HPP

#include <filesystem>
#include <vector>
#include <string>
#include "PhysicalSystemBase.hpp"
#include "Main/Globals.hpp"

class BinaryAlloyHexagonal2D : public PhysicalSystem {

public :

  BinaryAlloyHexagonal2D(const char* alloyConfigFile = "config_initial.dat", int = 0); 
  ~BinaryAlloyHexagonal2D();

  //void readCommandLineOptions()                         override;
  void writeConfiguration(int = 0, const char* = NULL)  override;
  void getObservables()                                 override;
  void doMCMove()                                       override;
  void acceptMCMove()                                   override;
  void rejectMCMove()                                   override;

  //void buildMPIConfigurationType()                      override;

private :

  unsigned int Size;

  std::string species[2]; // name of species A and B
  double concentration[2]; //concentration of species A and B
  unsigned int numPerSpecies[2];
  
  double singleSiteEnergy[2]; // one atom cluster energies for A and B
  double nnInteraction[3]; // interaction energies for A-A, A-B and B-B nearest neighbour bonds
  double nnInteractionMatrix[2][2];
  double nnnInteraction[3]; // interaction energies for A-A, A-B and B-B next nearest neighbour bonds
  double nnnInteractionMatrix[2][2];
  double threeSiteInteraction[4]; // interaction energies for A-A-A, A-A-B, A-B-B and B-B-B nearest neighbour triangle clusters
  double triangleInteractionTensor[2][2][2];
  
  typedef int SiteOccupation;  

  // Old configuration
  unsigned int CurX[2], CurY[2];
  SiteOccupation CurType[2];

  // New configuration
  SiteOccupation** siteOccupation;          // 2D array because it is a 2D model
  
  bool firstTimeGetMeasures;

  // Private functions
  ObservableType                                                             getSingleSiteEnergy();
  ObservableType                                                             getNNEnergy();
  ObservableType                                                             getNNNEnergy();
  ObservableType                                                             getTriangleEnergy();
  
  ObservableType                                                             getDifferenceInNNEnergy();
  ObservableType                                                             getDifferenceInNNNEnergy();
  ObservableType                                                             getDifferenceInTriangleEnergy();

  ObservableType getNNEnergyDelta(unsigned int i, unsigned int j, int s00);
  ObservableType getNNNEnergyDelta(unsigned int i, unsigned int j, int s00);
  ObservableType getTriangleEnergyDelta(unsigned int i, unsigned int j, int s00);
  
  void readHamiltonian(const char* inputFile);
  
  void readAlloyConfigFile(const std::filesystem::path& alloyConfigFile);
  void initializeAlloyConfiguration(int initial);

};

#endif

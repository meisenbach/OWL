#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <sstream>
#include "BinaryAlloyHexagonal2D.hpp"
#include "Utilities/RandomNumberGenerator.hpp"

BinaryAlloyHexagonal2D::BinaryAlloyHexagonal2D(const char* alloyConfigFile, int initial)
{

  printf("Simulation for 2D Binary Alloy on a hexagonal lattice: %dx%d \n", simInfo.modelLatticeSize, simInfo.modelLatticeSize);

  Size = simInfo.modelLatticeSize;
  setSystemSize(Size * Size);

  readHamiltonian(simInfo.MCInputFile);

  // printf("start allocate\n");
  
  siteOccupation = new SiteOccupation*[Size];

  siteOccupation[0] = new SiteOccupation[Size*Size];
  for (unsigned int i = 1; i < Size; i++) 
    siteOccupation[i] = &siteOccupation[0][i * Size];

  // printf("allocated\n");
  
  if (std::filesystem::exists(alloyConfigFile))
    readAlloyConfigFile(alloyConfigFile);
  else if (simInfo.restartFlag && std::filesystem::exists("configurations/config_checkpoint.dat"))
    readAlloyConfigFile("configurations/config_checkpoint.dat");
  else
    {
      if(initial == 0)
	initializeAlloyConfiguration(simInfo.configInitMethod);
      else
	initializeAlloyConfiguration(initial);
    }
 
  initializeObservables(1);
  observableName.push_back("Total energy, E");                            // observables[0] : total energy
  // observableName.push_back("Magnetization in x-direction, M_x");          // observables[1] : magnetization in x-direction
  // observableName.push_back("Magnetization in y-direction, M_y");          // observables[2] : magnetization in y-direction
  // observableName.push_back("Magnetization in z-direction, M_z");          // observables[3] : magnetization in z-direction
  // observableName.push_back("Total magnetization, M");                     // observables[4] : total magnetization

  firstTimeGetMeasures = true;
  getObservables();

}

void BinaryAlloyHexagonal2D::readHamiltonian(const char* mainInputFile)
{
   //if (GlobalComm.thisMPIrank == 0)
  std::cout << "\n   BinaryAlloyHexagonal2D class reading input file: " << mainInputFile << "\n\n";

  std::ifstream inputFile(mainInputFile);
  std::string line, key;

  species[0] = "A";
  species[1] = "B";
  concentration[0] = concentration[1] = 1.0;
  
  singleSiteEnergy[0] = singleSiteEnergy[1] = 0.0;
  nnInteraction[0] = nnInteraction[1] = nnInteraction[2] = 0.0;
  nnnInteraction[0] = nnnInteraction[1] = nnnInteraction[2] = 0.0;
  threeSiteInteraction[0] = threeSiteInteraction[1] = threeSiteInteraction[2] = threeSiteInteraction[3] = 0.0;

  if (inputFile.is_open()) {

    while (std::getline(inputFile, line)) {

      if (!line.empty()) {
	std::istringstream lineStream(line);
	lineStream >> key;

	if (key.compare(0, 1, "#") != 0) {
	  if (key == "Concentration") {
	    lineStream >> concentration[0] >> concentration[1];
            continue;
	  }
	  else if (key == "AtomSpecies") {
	    lineStream >> species[0] >> species[1];
            continue;
	  }
	  else if (key == "SingleSiteEnergy") {
	    lineStream >> singleSiteEnergy[0] >> singleSiteEnergy[1];
            continue;
	  }
	  else if (key == "NNInteraction") {
	    lineStream >> nnInteraction[0] >> nnInteraction[1] >> nnInteraction[2];
            continue;
	  }
	  else if (key == "NNNInteraction") {
	    lineStream >> nnnInteraction[0] >> nnnInteraction[1] >> nnnInteraction[2];
            continue;
	  }
	  else if (key == "ThreeSiteInteraction") {
	    lineStream >> threeSiteInteraction[0] >> threeSiteInteraction[1] >> threeSiteInteraction[2] >> threeSiteInteraction[3];
            continue;
	  }
	}
      }
    }
    
    inputFile.close();
  }

  nnInteractionMatrix[0][0] = nnInteraction[0];
  nnInteractionMatrix[0][1] = nnInteractionMatrix[1][0] = nnInteraction[1];
  nnInteractionMatrix[1][1] = nnInteraction[2];

  nnnInteractionMatrix[0][0] = nnnInteraction[0];
  nnnInteractionMatrix[0][1] = nnnInteractionMatrix[1][0] = nnnInteraction[1];
  nnnInteractionMatrix[1][1] = nnnInteraction[2];

  triangleInteractionTensor[0][0][0] = threeSiteInteraction[0];
  triangleInteractionTensor[1][0][0] = triangleInteractionTensor[0][1][0] = triangleInteractionTensor[0][0][1] = threeSiteInteraction[1];
  triangleInteractionTensor[0][1][1] = triangleInteractionTensor[1][0][1] = triangleInteractionTensor[1][1][0] = threeSiteInteraction[2];
  triangleInteractionTensor[1][1][1] = threeSiteInteraction[3];
  

  printf("Species: %s %s\n", species[0].c_str(), species[1].c_str());
  printf("Concentrations: %f %f\n", concentration[0], concentration[1]);
  printf("Single Site Energies : %f %f\n", singleSiteEnergy[0], singleSiteEnergy[1]);
  printf("Nearest Neighbor Interactions :  %f %f %f\n", nnInteraction[0], nnInteraction[1], nnInteraction[2]);
  printf("Next Nearest Neighbor Interactions :  %f %f %f\n", nnnInteraction[0], nnnInteraction[1], nnnInteraction[2]);
  printf("Three Site Interactions : %f %f %f %f\n\n", threeSiteInteraction[0], threeSiteInteraction[1], threeSiteInteraction[2], threeSiteInteraction[3]);
  
}

BinaryAlloyHexagonal2D::~BinaryAlloyHexagonal2D()
{
  delete[] siteOccupation[0];
  delete[] siteOccupation;

  printf("BinaryAlloyHexagonal2D finished\n");
}


//void HeisenbergHexagonal2D::readCommandLineOptions()
//{ };


void BinaryAlloyHexagonal2D::writeConfiguration(int format, const char* filename)
{

  FILE* f;
  if (filename != NULL) f = fopen(filename, "w");
  else f = stdout;

  switch (format) {

  case 10 : {
    fprintf(f, "# 2D Hexagonal Binary Alloy Model : %u x %u \n\n", Size, Size);
    fprintf(f, "# TotalNumberOfSites %u\n", systemSize);
    fprintf(f, "# Atoms per Species %u %s, %u %s\n",numPerSpecies[0], species[0].c_str(), numPerSpecies[1], species[1].c_str());
    fprintf(f, "# Energy: %10.5f\n\n", observables[0]);
    
    // fprintf(f, "Observables ");
    //
    // for (unsigned int i = 0; i < numObservables; i++)
    //   fprintf(f, " %10.5f", observables[i]);
    // fprintf(f, "\n");

    // fprintf(f, "\nAlloyConfiguration\n");
    for (unsigned int i = 0; i < Size; i++) {
      for (unsigned int j = 0; j < Size; j++) {
	double x = double(i) - 0.5 * double(j);
	double y = 0.8660 * double(j);
        fprintf(f, "%4d %4d  %10.5f %10.5f  %2d\n", i, j, x, y, siteOccupation[i][j]);
      }
    }
    break;
  }

  default : {
    fprintf(f, "# 2D Hexagonal Binary Alloy Model : %u x %u \n\n", Size, Size);
    fprintf(f, "TotalNumberOfSites %u\n", systemSize);
    fprintf(f, "Observables ");

    for (unsigned int i = 0; i < numObservables; i++)
      fprintf(f, " %10.5f", observables[i]);
    fprintf(f, "\n");

    fprintf(f, "\nAlloyConfiguration\n");
    for (unsigned int i = 0; i < Size; i++) {
      for (unsigned int j = 0; j < Size; j++)
        fprintf(f, "%2d\n", siteOccupation[i][j]);
    }

  }

  }

  if (filename != NULL) fclose(f);

}


void BinaryAlloyHexagonal2D::getObservables()
{

  if (firstTimeGetMeasures) {
    //resetObservables();
    observables[0] = getSingleSiteEnergy() + getNNEnergy() + getNNNEnergy() + getTriangleEnergy();

    firstTimeGetMeasures = false;
    //printf("First time getObservables. \n");
  }
  else {
    observables[0] = getSingleSiteEnergy() + getNNEnergy() + getNNNEnergy() + getTriangleEnergy();
    // observables[0] += getDifferenceInNNEnergy() + getDifferenceInNNNEnergy() + getDifferenceInTriangleEnergy();

    //printf("observables = %10.5f %10.5f %10.5f %10.5f %10.5f\n", observables[0], observables[1], observables[2], observables[3], observables[4]);
  }

}

ObservableType BinaryAlloyHexagonal2D::getSingleSiteEnergy()
{
  ObservableType energy {0.0};
  for (unsigned int i = 0; i < Size; i++)
    {
      for (unsigned int j = 0; j < Size; j++)
	{
	  energy += singleSiteEnergy[siteOccupation[i][j]];
	}
    }
  return energy;
}

ObservableType BinaryAlloyHexagonal2D::getNNEnergyDelta(unsigned int i, unsigned int j, int s00)
{
  unsigned int xPlus, xMinus, yPlus, yMinus;
  ObservableType dEnergy;

  // nearest neighbor
  
  if (i != 0) xMinus = i - 1; else xMinus = Size - 1;
  if (i != Size - 1) xPlus = i + 1; else xPlus = 0;
  
  if (j != 0) yMinus = j - 1; else yMinus = Size - 1;
  if (j != Size - 1) yPlus = j + 1; else yPlus = 0;

  // (-1, -1)
  dEnergy = nnInteraction[s00 + siteOccupation[xMinus][yMinus]];
  // (-1, 0)
  dEnergy += nnInteraction[s00 + siteOccupation[xMinus][j]];
  // (0, -1)
  dEnergy += nnInteraction[s00 + siteOccupation[i][yMinus]];
  // (1, 0)
  dEnergy += nnInteraction[s00 + siteOccupation[xPlus][j]];
  // (0, 1)
  dEnergy += nnInteraction[s00 + siteOccupation[i][yPlus]];
  // (1, 1)
  dEnergy += nnInteraction[s00 + siteOccupation[xPlus][yPlus]];
  
  return dEnergy;
}

ObservableType BinaryAlloyHexagonal2D::getNNEnergy()
{
  unsigned int xPlus, xMinus, yPlus, yMinus;
  ObservableType energy {0.0};
  ObservableType dEnergy;

  // nearest neighbor
  
  for (unsigned int i = 0; i < Size; i++)
    {
    if (i != 0) xMinus = i - 1; else xMinus = Size - 1;
    if (i != Size - 1) xPlus = i + 1; else xPlus = 0;
    for (unsigned int j = 0; j < Size; j++)
      {
	if (j != 0) yMinus = j - 1; else yMinus = Size - 1;
	if (j != Size - 1) yPlus = j + 1; else yPlus = 0;

	int s00 = siteOccupation[i][j];

	// (-1, -1)
	dEnergy = nnInteraction[s00 + siteOccupation[xMinus][yMinus]];
	// (-1, 0)
	dEnergy += nnInteraction[s00 + siteOccupation[xMinus][j]];
	// (0, -1)
	dEnergy += nnInteraction[s00 + siteOccupation[i][yMinus]];
	// (1, 0)
	dEnergy += nnInteraction[s00 + siteOccupation[xPlus][j]];
	// (0, 1)
	dEnergy += nnInteraction[s00 + siteOccupation[i][yPlus]];
	// (1, 1)
	dEnergy += nnInteraction[s00 + siteOccupation[xPlus][yPlus]];
	
        energy += 0.5*dEnergy;
    }  
  }
  
  return energy;
}

ObservableType BinaryAlloyHexagonal2D::getNNNEnergyDelta(unsigned int i, unsigned int j, int s00)
{
  unsigned int xPlus, xMinus, yPlus, yMinus;
  unsigned int xPlus2, xMinus2, yPlus2, yMinus2;
  ObservableType dEnergy;

  // next nearest neighbor

  if (i != 0) xMinus = i - 1; else xMinus = Size - 1;
  if (i != Size - 1) xPlus = i + 1; else xPlus = 0;
    
  if (i > 1) xMinus2 = i - 2;
  else if (Size > 1) xMinus2 = Size + i - 2; else xMinus2 = 0;
  if (i < Size - 2) xPlus2 = i + 2;
  else if (Size > 1) xPlus2 = i - Size + 2; else xPlus2 = 0;
    
  if (j != 0) yMinus = j - 1; else yMinus = Size - 1;
  if (j != Size - 1) yPlus = j + 1; else yPlus = 0;
  
  if (j > 1) yMinus2 = j - 2;
  else if (Size > 1) yMinus2 = Size + j - 2; else yMinus2 = 0;
  if (j < Size - 2) yPlus2 = j + 2;
  else if (Size > 1) yPlus2 = j - Size + 2; else yPlus2 = 0;

	// (-2, -1)
  dEnergy = nnnInteraction[s00 + siteOccupation[xMinus2][yMinus]];
	// (-1, -2)
  dEnergy += nnnInteraction[s00 + siteOccupation[xMinus][yMinus2]];
	// (-1, 1)
  dEnergy += nnnInteraction[s00 + siteOccupation[xMinus][yPlus]];
	// (1, -1)
  dEnergy += nnnInteraction[s00 + siteOccupation[xPlus][yMinus]];
	// (1, 2)
  dEnergy += nnnInteraction[s00 + siteOccupation[xPlus][yPlus2]];
	// (2, 1)
  dEnergy += nnnInteraction[s00 + siteOccupation[xPlus2][yPlus]];
  
  return dEnergy;
}

ObservableType BinaryAlloyHexagonal2D::getNNNEnergy()
{
  unsigned int xPlus, xMinus, yPlus, yMinus;
  unsigned int xPlus2, xMinus2, yPlus2, yMinus2;
  ObservableType energy {0.0};
  ObservableType dEnergy;

  // next nearest neighbor
  
  for (unsigned int i = 0; i < Size; i++)
    {
      if (i != 0) xMinus = i - 1; else xMinus = Size - 1;
      if (i != Size - 1) xPlus = i + 1; else xPlus = 0;
    
      if (i > 1) xMinus2 = i - 2;
      else if (Size > 1) xMinus2 = Size + i - 2; else xMinus2 = 0;
      if (i < Size - 2) xPlus2 = i + 2;
      else if (Size > 1) xPlus2 = i - Size + 2; else xPlus2 = 0;
    
      for (unsigned int j = 0; j < Size; j++)
	{
	  if (j != 0) yMinus = j - 1; else yMinus = Size - 1;
	  if (j != Size - 1) yPlus = j + 1; else yPlus = 0;
	  
	  if (j > 1) yMinus2 = j - 2;
	  else if (Size > 1) yMinus2 = Size + j - 2; else yMinus2 = 0;
	  if (j < Size - 2) yPlus2 = j + 2;
	  else if (Size > 1) yPlus2 = j - Size + 2; else yPlus2 = 0;

	  int s00 = siteOccupation[i][j];
	// (-2, -1)
	  dEnergy = nnnInteraction[s00 + siteOccupation[xMinus2][yMinus]];
	// (-1, -2)
	  dEnergy += nnnInteraction[s00 + siteOccupation[xMinus][yMinus2]];
	// (-1, 1)
	  dEnergy += nnnInteraction[s00 + siteOccupation[xMinus][yPlus]];
	// (1, -1)
	  dEnergy += nnnInteraction[s00 + siteOccupation[xPlus][yMinus]];
	// (1, 2)
	  dEnergy += nnnInteraction[s00 + siteOccupation[xPlus][yPlus2]];
	// (2, 1)
	  dEnergy += nnnInteraction[s00 + siteOccupation[xPlus2][yPlus]];
	
	  energy += 0.5*dEnergy;
	}  
    }
  
  return energy;
}

ObservableType BinaryAlloyHexagonal2D::getTriangleEnergyDelta(unsigned int i, unsigned int j, int s00)
{
  unsigned int xPlus, xMinus, yPlus, yMinus;
  ObservableType dEnergy;

  // nearest neighbor
  
  if (i != 0) xMinus = i - 1; else xMinus = Size - 1;
  if (i != Size - 1) xPlus = i + 1; else xPlus = 0;

  if (j != 0) yMinus = j - 1; else yMinus = Size - 1;
  if (j != Size - 1) yPlus = j + 1; else yPlus = 0;

  int sMM = siteOccupation[xMinus][yMinus];
  int s0M = siteOccupation[i][yMinus];
  int sP0 = siteOccupation[xPlus][j];
  int sPP = siteOccupation[xPlus][yPlus];
  int s0P = siteOccupation[i][yPlus];
  int sM0 = siteOccupation[xMinus][j];
	
	// (-1, -1) (0, -1)
  dEnergy = threeSiteInteraction[s00 + sMM + s0M];
	// (0, -1) (1, 0)
  dEnergy += threeSiteInteraction[s00 + s0M + sP0];
	// (1, 0) (1, 1)
  dEnergy += threeSiteInteraction[s00 + sP0 + sPP];
	// (1, 1) (0, 1)
  dEnergy += threeSiteInteraction[s00 + sPP + s0P];
	// (0, 1) (-1, 0)
  dEnergy += threeSiteInteraction[s00 + s0P + sM0];
	// (-1, 0) (-1, -1)
  dEnergy += threeSiteInteraction[s00 + sM0 + sMM];
  
  return dEnergy;
}

ObservableType BinaryAlloyHexagonal2D::getTriangleEnergy()
{
    unsigned int xPlus, xMinus, yPlus, yMinus;
  ObservableType energy {0.0};
  ObservableType dEnergy;

  // nearest neighbor
  
  for (unsigned int i = 0; i < Size; i++)
    {
    if (i != 0) xMinus = i - 1; else xMinus = Size - 1;
    if (i != Size - 1) xPlus = i + 1; else xPlus = 0;
    for (unsigned int j = 0; j < Size; j++)
      {
	if (j != 0) yMinus = j - 1; else yMinus = Size - 1;
	if (j != Size - 1) yPlus = j + 1; else yPlus = 0;

	int s00 = siteOccupation[i][j];
	int sMM = siteOccupation[xMinus][yMinus];
	int s0M = siteOccupation[i][yMinus];
	int sP0 = siteOccupation[xPlus][j];
	int sPP = siteOccupation[xPlus][yPlus];
	int s0P = siteOccupation[i][yPlus];
	int sM0 = siteOccupation[xMinus][j];
	
	// (-1, -1) (0, -1)
	dEnergy = threeSiteInteraction[s00 + sMM + s0M];
	// (0, -1) (1, 0)
	dEnergy += threeSiteInteraction[s00 + s0M + sP0];
	// (1, 0) (1, 1)
	dEnergy += threeSiteInteraction[s00 + sP0 + sPP];
	// (1, 1) (0, 1)
	dEnergy += threeSiteInteraction[s00 + sPP + s0P];
	// (0, 1) (-1, 0)
	dEnergy += threeSiteInteraction[s00 + s0P + sM0];
	// (-1, 0) (-1, -1)
	dEnergy += threeSiteInteraction[s00 + sM0 + sMM];
	
        energy += (1.0/3.0)*dEnergy;
    }  
  }
  
  return energy;
}


ObservableType BinaryAlloyHexagonal2D::getDifferenceInNNEnergy()
{
  ObservableType energyChange {0.0};
  energyChange = getNNEnergyDelta(CurX[0], CurY[0], siteOccupation[CurX[0]][CurY[0]])
    + getNNEnergyDelta(CurX[1], CurY[1], siteOccupation[CurX[1]][CurY[1]])
    - getNNEnergyDelta(CurX[0], CurY[0], CurType[0])
    - getNNEnergyDelta(CurX[1], CurY[1], CurType[1]);
  return energyChange;
}

ObservableType BinaryAlloyHexagonal2D::getDifferenceInNNNEnergy()
{
  ObservableType energyChange {0.0};
  energyChange = getNNNEnergyDelta(CurX[0], CurY[0], siteOccupation[CurX[0]][CurY[0]])
    + getNNNEnergyDelta(CurX[1], CurY[1], siteOccupation[CurX[1]][CurY[1]])
    - getNNNEnergyDelta(CurX[0], CurY[0], CurType[0])
    - getNNNEnergyDelta(CurX[1], CurY[1], CurType[1]);
  return energyChange;
}

ObservableType BinaryAlloyHexagonal2D::getDifferenceInTriangleEnergy()
{
  ObservableType energyChange {0.0};
  energyChange = getTriangleEnergyDelta(CurX[0], CurY[0], siteOccupation[CurX[0]][CurY[0]])
    + getTriangleEnergyDelta(CurX[1], CurY[1], siteOccupation[CurX[1]][CurY[1]])
    - getTriangleEnergyDelta(CurX[0], CurY[0], CurType[0])
    - getTriangleEnergyDelta(CurX[1], CurY[1], CurType[1]);
  return energyChange;
}

void BinaryAlloyHexagonal2D::doMCMove()
{
  // Need this here since resetObservables() is not called if firstTimeGetMeasures = false
  //for (int i = 0; i < numObservables; i++)
  //  oldObservables[i] = observables[i];

  CurX[0] = unsigned(getIntRandomNumber()) % Size;
  CurY[0] = unsigned(getIntRandomNumber()) % Size;
  CurX[1] = unsigned(getIntRandomNumber()) % Size;
  CurY[1] = unsigned(getIntRandomNumber()) % Size;
  while(CurX[0] == CurX[1] && CurY[0] == CurY[1])
    {
      CurX[1] = unsigned(getIntRandomNumber()) % Size;
      CurY[1] = unsigned(getIntRandomNumber()) % Size;
    }

  CurType[0] = siteOccupation[CurX[0]][CurY[0]];
  CurType[1] = siteOccupation[CurX[1]][CurY[1]];

  siteOccupation[CurX[0]][CurY[0]] = CurType[1];
  siteOccupation[CurX[1]][CurY[1]] = CurType[0];
}


/*
void BinaryAlloyHexagonal2D::undoMCMove()
{
  siteOccupation[CurX[0]][CurY[0]] = CurType[0];
  siteOccupation[CurX[1]][CurY[1]] = CurType[1];
  restoreObservables();
}
*/

void BinaryAlloyHexagonal2D::acceptMCMove()
{
  // update "old" observables
  for (unsigned int i=0; i < numObservables; i++)
    oldObservables[i] = observables[i];
}


void BinaryAlloyHexagonal2D::rejectMCMove()
{
  siteOccupation[CurX[0]][CurY[0]] = CurType[0];
  siteOccupation[CurX[1]][CurY[1]] = CurType[1];
  
  for (unsigned int i=0; i < numObservables; i++)
    observables[i] = oldObservables[i];
}

/*
void HeisenbergHexagonal2D::buildMPIConfigurationType()
{
}
*/


void BinaryAlloyHexagonal2D::readAlloyConfigFile(const std::filesystem::path& alloyConfigFile)
{

  std::cout << "\n   BinaryAlloyHexagonal2D class reading configuration file: " << alloyConfigFile << "\n";

  std::ifstream inputFile(alloyConfigFile);
  std::string line, key;
  unsigned int numberOfSites {0};

  if (inputFile.is_open()) {

    while (std::getline(inputFile, line)) {

      if (!line.empty()) {
        std::istringstream lineStream(line);
        lineStream >> key;

        if (key.compare(0, 1, "#") != 0) {

          if (key == "TotalNumberOfSites") {
            lineStream >> numberOfSites;
            //std::cout << "   BinaryAlloyHexagonal2D: numberOfSites = " << numberOfSites << "\n";
            continue;
          }
          else if (key == "Observables") {
            unsigned int counter = 0;
            while (lineStream && counter < numObservables) {
              lineStream >> observables[counter];
              //std::cout << "   BinaryAlloyHexagonal2D: observables[" << counter << "] = " << observables[counter] << "\n";
              counter++;
            }
            continue;
          }
          else if (key == "AlloyConfiguration") {
            //std::cout << "   BinaryAlloyHexagonal2D: Alloy Configuration read: \n";
            for (unsigned int i=0; i<Size; i++) {
              for (unsigned int j=0; j<Size; j++) {
                lineStream.clear();
                std::getline(inputFile, line);               
                if (!line.empty())  lineStream.str(line);
                lineStream >> siteOccupation[i][j];
                //printf("      %2df\n", siteOccupation[i][j]);
              }
            }
            continue;
          }
        }

      }
    }

    inputFile.close();
  }

  // Sanity checks:
  assert(numberOfSites == systemSize);
  
  printf("   Initial configuration read:\n");
  for (unsigned int i=0; i<Size; i++) {
    for (unsigned int j=0; j<Size; j++)
      printf("      %2d\n", siteOccupation[i][j]);
    printf("\n");
  }

}


void BinaryAlloyHexagonal2D::initializeAlloyConfiguration(int initial)
{

  int numA = Size * Size * concentration[0] / (concentration[0] + concentration[1]);
  int numB = Size * Size - numA;
  numPerSpecies[0] = numA;
  numPerSpecies[1] = Size * Size - numA;

  // printf("Initializing Alloy: %d %s, %d %s\n", numA, species[0].c_str(), Size * Size - numA, species[1].c_str());

  switch (initial)
    {
    case 1: // Honeycomb
      numA = numB = 0;

      for (unsigned int i = 0; i < Size; i++)
	for (unsigned int j = 0; j < Size; j++)
	  if((i + j) % 3 != 0) {
	    siteOccupation[i][j] = 0;
	    numA++;
	  } else {
	    siteOccupation[i][j] = 1;
	    numB++;
	  }
      
      numPerSpecies[0] = numA;
      numPerSpecies[1] = numB;
      concentration[0] = double(numA)/double(Size * Size);
      concentration[1] = double(numB)/double(Size * Size);
      printf("Initializing Alloy: %d %s, %d %s\n", numA, species[0].c_str(), Size * Size - numA, species[1].c_str());
      printf("  Honeycomb (Fixed concentration!)\n");
      
      break;
    case 11: // Honeycomb
      numA = numB = 0;

      for (unsigned int i = 0; i < Size; i++)
	for (unsigned int j = 0; j < Size; j++)
	  if((i + j) % 3 != 0) {
	    siteOccupation[i][j] = 1;
	    numB++;
	  } else {
	    siteOccupation[i][j] = 0;
	    numA++;
	  }
      
      numPerSpecies[0] = numA;
      numPerSpecies[1] = numB;
      concentration[0] = double(numA)/double(Size * Size);
      concentration[1] = double(numB)/double(Size * Size);
      printf("Initializing Alloy: %d %s, %d %s\n", numA, species[0].c_str(), Size * Size - numA, species[1].c_str());
      printf("  Honeycomb (Fixed concentration!)\n");
      
      break;
    case 2: // Kagome
      numA = numB = 0;
      
      for (unsigned int i = 0; i < Size; i++)
	for (unsigned int j = 0; j < Size; j++)
	  if((i % 2 == 1) || (j % 2 == 0)) {
	    siteOccupation[i][j] = 0;
	    numA++;
	  } else {
	    siteOccupation[i][j] = 1;
	    numB++;
	  }

      numPerSpecies[0] = numA;
      numPerSpecies[1] = numB;
      concentration[0] = double(numA)/double(Size * Size);
      concentration[1] = double(numB)/double(Size * Size);
      printf("Initializing Alloy: %d %s, %d %s\n", numA, species[0].c_str(), Size * Size - numA, species[1].c_str());
	    printf("  Kagome (Fixed concentration!)\n");

      break;
    case 12: // Kagome
      numA = numB = 0;
      
      for (unsigned int i = 0; i < Size; i++)
	for (unsigned int j = 0; j < Size; j++)
	  if((i % 2 == 1) || (j % 2 == 0)) {
	    siteOccupation[i][j] = 1;
	    numB++;
	  } else {
	    siteOccupation[i][j] = 0;
	    numA++;
	  }

      numPerSpecies[0] = numA;
      numPerSpecies[1] = numB;
      concentration[0] = double(numA)/double(Size * Size);
      concentration[1] = double(numB)/double(Size * Size);
      printf("Initializing Alloy: %d %s, %d %s\n", numA, species[0].c_str(), Size * Size - numA, species[1].c_str());
      printf("  Kagome (Fixed concentration!)\n");

      break;
    default:
      printf("Initializing Alloy: %d %s, %d %s\n", numA, species[0].c_str(), Size * Size - numA, species[1].c_str());
      printf("  Random initial state.\n");
      for (unsigned int i = 0; i < Size; i++)
	for (unsigned int j = 0; j < Size; j++)
	  if(numA > 0) {
	    siteOccupation[i][j] = 0;
	    numA--;
	  } else {
	    siteOccupation[i][j] = 1;
	  }
      
      for (unsigned int i = 0; i < Size; i++)
	for (unsigned int j = 0; j < Size; j++)
	  {
	    unsigned int x = unsigned(getIntRandomNumber()) % Size;
	    unsigned int y = unsigned(getIntRandomNumber()) % Size;

	    int oc = siteOccupation[i][j];
	    siteOccupation[i][j] = siteOccupation[x][y];
	    siteOccupation[x][y] = oc;
	  }
    }
}

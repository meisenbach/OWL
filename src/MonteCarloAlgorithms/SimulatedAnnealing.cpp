#include <cassert>
#include <cmath>
#include <filesystem>
#include <sstream>
#include "SimulatedAnnealing.hpp"
#include "Utilities/RandomNumberGenerator.hpp"
#include "Utilities/CheckFile.hpp"
#include "Utilities/CompareNumbers.hpp"

// Constructor
SimulatedAnnealing::SimulatedAnnealing(PhysicalSystem* ps, const char* inputFile)
{

  bestObservable0 = std::numeric_limits<double>::max();
  
  if (GlobalComm.thisMPIrank == 0)
    printf("\nSimulation method: Simulated Annealing\n");

  if (std::filesystem::exists(inputFile))
    readMCInputFile(inputFile);
  else {
    std::cout << "Error: No input file for reading Simulated Annealing info. Quiting... \n";
    exit(7);
  }

  std::filesystem::create_directory("best_configurations");
  
  physical_system = ps;

  // Allocate space to store observables and other statistics
  if (physical_system->numObservables > 0) {
    averagedObservables.assign(physical_system->numObservables, 0.0);
    averagedObservablesSquared.assign(physical_system->numObservables, 0.0);
    standardDeviations.assign(physical_system->numObservables, 0.0);
    standardErrors.assign(physical_system->numObservables, 0.0);
  }

  if (std::filesystem::exists("mc.dat")) 
    timeSeriesFile = fopen("mc.dat", "a");
  else 
    timeSeriesFile = fopen("mc.dat", "w"); 

  // fprintf(timeSeriesFile, "# Thermalization: (%lu steps) \n", numberOfThermalizationSteps);
  // fprintf(timeSeriesFile, "# Temperature %8.5f\n", temperature);
  // fprintf(timeSeriesFile, "# Sequence# Temperature  Energy  MCsteps   Observables\n");

  if (simInfo.restartFlag) {
    if (std::filesystem::exists("simulated_annealing_checkpoint.dat"))
      readCheckPointFile("simulated_annealing_checkpoint.dat");
    else {
      std::cout << "\n   WARNING! Restart file 'simulated_annealing_checkpoint.dat' not found. ";
      std::cout << "\n            Performing a fresh run instead of a restarted run. \n\n";
    }
  }

}


//Destructor
SimulatedAnnealing::~SimulatedAnnealing()
{

  fclose(timeSeriesFile);

  if (GlobalComm.thisMPIrank == 0)
    printf("Exiting Simulated Annealing class... \n");

}


void SimulatedAnnealing::run() 
{

  char fileName[51];


  currentTime = lastBackUpTime = MPI_Wtime();
  if (GlobalComm.thisMPIrank == 0)
    {
      printf("   Running Simulated Annealing...\n");
      printf("   Temperature Reduction Steps: %d\n", numberOfTemperatureSteps);
    }

  // initial state:
  physical_system -> getObservables();
  bestSequenceNumber++;
  bestObservable0 = physical_system -> observables[0];
  fprintf(timeSeriesFile, "# INITIAL configuration # %d @ %f\n", bestSequenceNumber, bestObservable0);
  fprintf(timeSeriesFile, "%5d %f  %f  ", bestSequenceNumber, temperature, bestObservable0);
  writeMCFile(thermalizationStepsPerformed);
  sprintf(fileName, "best_configurations/initital%05d.dat", bestSequenceNumber);
  physical_system -> writeConfiguration(10, fileName);
  //

  while (numberOfTemperatureSteps > 0) {
    if (GlobalComm.thisMPIrank == 0)
      printf("   Temperature: %f (%d temperature steps remaining)\n",temperature, numberOfTemperatureSteps);

    clearObservables();
    physical_system -> getObservables();
    physical_system -> oldObservables[0] = physical_system -> observables[0];
    rejectedMoves = acceptedMoves = 0;
    /*
  // Thermalization (observables are not accumulated)
    while(false) {
      //  while (thermalizationStepsPerformed < numberOfThermalizationSteps) {
  //for (unsigned long int MCSteps=0; MCSteps<numberOfThermalizationSteps; MCSteps++) {

      for (unsigned long int i=0; i<numberOfMCUpdatesPerStep; i++) {

	physical_system -> doMCMove();
	physical_system -> getObservables();
  
	// Determine acceptance
	if ( exp((physical_system -> oldObservables[0] - physical_system -> observables[0]) / temperature ) > getRandomNumber2() )
	  physical_system -> acceptMCMove(); 
	else
	  physical_system -> rejectMCMove();

	// save best configuration found
	if (physical_system -> observables[0] < bestObservable0)
	  {
	    bestSequenceNumber++;
	    bestObservable0 = physical_system -> observables[0];
	    // fprintf(timeSeriesFile, "# NEW best configuration # %d @ %f\n", bestSequenceNumber, bestObservable0);
	    fprintf(timeSeriesFile, "%5d %f  %f  ", bestSequenceNumber, temperature, bestObservable0);
	    writeMCFile(thermalizationStepsPerformed);
	    sprintf(fileName, "best_configurations/best%05d.dat", bestSequenceNumber);
	    physical_system -> writeConfiguration(10, fileName);
	  }
      }

      thermalizationStepsPerformed++;
      physical_system -> getAdditionalObservables();

      // Write observables to file
      // writeMCFile(thermalizationStepsPerformed);
    
      // Write restart files at interval
      currentTime = MPI_Wtime();
      if (GlobalComm.thisMPIrank == 0) {
	if (currentTime - lastBackUpTime > checkPointInterval) {
	  writeCheckPointFiles(checkPoint);
	  lastBackUpTime = currentTime;
	}
      }
 
    }

    // fprintf(timeSeriesFile, "# End of thermalization. \n\n");
    */
    // fprintf(timeSeriesFile, "# Accumulation: (%lu steps) \n", numberOfMCSteps);
    fprintf(timeSeriesFile, "# Temperature %8.5f\n", temperature);
    // fprintf(timeSeriesFile, "# MC steps           Observables\n");
    MCStepsPerformed = 0;
    // Observable accumulation starts here
    while (MCStepsPerformed < numberOfMCSteps) {
  //for (unsigned long int MCSteps=0; MCSteps<numberOfMCSteps; MCSteps++) {

      for (unsigned long int i=0; i<numberOfMCUpdatesPerStep; i++) {
    
	physical_system -> doMCMove();
	physical_system -> getObservables();
  
	// Determine acceptance
	if ( exp((physical_system -> oldObservables[0] - physical_system -> observables[0]) / temperature ) > getRandomNumber2() ) {
	  physical_system -> acceptMCMove();
	  acceptedMoves++;
	}
	else {
	  physical_system -> rejectMCMove();
	  rejectedMoves++;
	}

	// save best configuration found
	if (physical_system -> observables[0] < bestObservable0)
	  {
	    bestSequenceNumber++;
	    bestObservable0 = physical_system -> observables[0];
	    // fprintf(timeSeriesFile, "# NEW best configuration # %d @ %f\n", bestSequenceNumber, bestObservable0);
	    printf("NEW best configuration # %d @ %f\n", bestSequenceNumber, bestObservable0);
	    fprintf(timeSeriesFile, "%5d %f  %f  ", bestSequenceNumber, temperature, bestObservable0);
	    writeMCFile(totalMCStepsPerformed);
	    sprintf(fileName, "best_configurations/best%05d.dat", bestSequenceNumber);
	    physical_system -> writeConfiguration(10, fileName);
	  }

      }
      MCStepsPerformed++;
      totalMCStepsPerformed++;

      physical_system -> getAdditionalObservables();
      accumulateObservables(); 

      // Write observables to file
      // writeMCFile(MCStepsPerformed);

      // Write restart files at interval
      currentTime = MPI_Wtime();
      if (GlobalComm.thisMPIrank == 0) {

	if (currentTime - lastBackUpTime > checkPointInterval) {
	  writeCheckPointFiles(checkPoint);
	  lastBackUpTime = currentTime;
	}

	if ((configurationWriteInterval != 0) && (totalMCStepsPerformed % configurationWriteInterval == 0)) {
	  sprintf(fileName, "configurations/config%012lu.dat", totalMCStepsPerformed);
	  physical_system -> writeConfiguration(0, fileName);
	  //sprintf(fileName, "configurations/config%012lu.xyz", MCStepsPerformed);
	  //physical_system -> writeConfiguration(1, fileName);
	  sprintf(fileName, "configurations/all-configs.dat");
	  physical_system -> writeConfiguration(2, fileName);
	}

      }

  }
  // fprintf(timeSeriesFile, "# End of accumulation. \n\n");
  fflush(timeSeriesFile);
  writeCheckPointFiles(checkPoint);

  calculateAveragesAndVariances();
  writeCheckPointFiles(endOfSimulation);
  
  // physical_system -> calculateThermodynamics(averagedObservables, averagedObservablesSquared, temperature);

  updateTemperature();
  numberOfTemperatureSteps--;
  }

}


// This implementation is similar to the one in Histogram class
void SimulatedAnnealing::readMCInputFile(const char* fileName)
{

  if (GlobalComm.thisMPIrank == 0) 
    std::cout << "   Simulated Annealing class reading input file: " << fileName << "\n";

  std::ifstream inputFile(fileName);   // TODO: check if a file stream is initialized
  std::string line, key;

  if (inputFile.is_open()) {
    
    while (std::getline(inputFile, line)) {

      if (!line.empty()) {
        
        std::istringstream lineStream(line);
        lineStream >> key;

        if (key.compare(0, 1, "#") != 0) {

	  /*
          if (key == "numberOfThermalizationSteps") {
            lineStream >> numberOfThermalizationSteps;
            //std::cout << "Metropolis: numberOfThermalizationSteps = " << numberOfThermalizationSteps << "\n";
            continue;
          }
          else 
	  */
	  if (key == "numberOfMCSteps") {
            lineStream >> numberOfMCSteps;
            //std::cout << "Metropolis: numberOfMCSteps = " << numberOfMCSteps << "\n";
            continue;
          }
          else if (key == "numberOfMCUpdatesPerStep") {
            lineStream >> numberOfMCUpdatesPerStep;
            //std::cout << "Metropolis: numberOfMCUpdatesPerStep = " << numberOfMCUpdatesPerStep << "\n";
            continue;
          }
	  else if (key == "numberOfTemperatureSteps") {
            lineStream >> numberOfTemperatureSteps;
            //std::cout << "Simulated Annealing: numberOfTemperatureSteps = " << numberOfTemperatureSteps << "\n";
            continue;
          }
          else if (key == "temperature") {
            lineStream >> temperature;
            //std::cout << "Metropolis: temperature = " << temperature << "\n";
            continue;
          }
          else if (key == "checkPointInterval") {
            lineStream >> checkPointInterval;
            //std::cout << "Metropolis: checkPointInterval = " << checkPointInterval << " seconds \n";
            continue;
          }
          else if (key == "configurationWriteInterval") {
            lineStream >> configurationWriteInterval;
            //std::cout << "Metropolis: configurationWriteInterval = " << configurationWriteInterval << " seconds \n";
            continue;
          }
	  else if (key == "bestSequenceNumber") {
            lineStream >> configurationWriteInterval;
            //std::cout << "Simulated Annealing: bestSequenceNumber = " << bestSequenceNumber << "\n";
            continue;
          }
	  else if (key == "bestObservable0") {
            lineStream >> bestObservable0;
            //std::cout << "Simulated Annealing: bestObservable0 = " << bestObservable0 << "\n";
            continue;
          }
          

        }

      }

    }
    inputFile.close();

  }

}


void SimulatedAnnealing::readCheckPointFile(const char* fileName)
{

  if (GlobalComm.thisMPIrank == 0) 
    std::cout << "   Simulated Annealing class reading checkpoint file: " << fileName << "\n";

  std::ifstream inputFile(fileName);   // TODO: check if a file stream is initialized
  std::string line, key;

  if (inputFile.is_open()) {
    
    while (std::getline(inputFile, line)) {

      if (!line.empty()) {
        
        std::istringstream lineStream(line);
        lineStream >> key;

        if (key.compare(0, 1, "#") != 0) {

          if (key == "temperature") {
            lineStream >> temperature;
            //std::cout << "Metropolis: restartTemperature = " << restartTemperature << "\n";
            continue;
          }
	  /*
          else if (key == "thermalizationStepsPerformed") {
            lineStream >> thermalizationStepsPerformed;
            //std::cout << "Metropolis: thermalizationStepsPerformed = " << thermalizationStepsPerformed << "\n";
            continue;
          }
	  */
          else if (key == "MCStepsPerformed") {
            lineStream >> MCStepsPerformed;
            //std::cout << "Metropolis: MCStepsPerformed = " << MCStepsPerformed << "\n";
            continue;
          }
          else if (key == "acceptedMoves") {
            lineStream >> acceptedMoves;
            //std::cout << "Metropolis: acceptedMoves = " << acceptedMoves << "\n";
            continue;
          }
          else if (key == "rejectedMoves") {
            lineStream >> rejectedMoves;
            //std::cout << "Metropolis: rejectedMoves = " << rejectedMoves << "\n";
            continue;
          }
          else if (key == "averagedObservables") {
            unsigned int counter = 0;
            while (lineStream && counter < physical_system->numObservables) {
              lineStream >> averagedObservables[counter];
              //std::cout << "Metropolis: averageObservables[" << counter << "] = " << averagedObservables[counter] << "\n";
              counter++;
            }
            continue;
          }
          else if (key == "standardErrors") {
            unsigned int counter = 0;
            while (lineStream && counter < physical_system->numObservables) {
              lineStream >> standardErrors[counter];
              //std::cout << "Metropolis: standardErrors[" << counter << "] = " << standardErrors[counter] << "\n";
              counter++;
            }
            continue;
          }
	  else if (key == "numberOfTemperatureSteps") {
            lineStream >> numberOfTemperatureSteps;
            //std::cout << "Simulated Annealing: numberOfTemperatureSteps = " << numberOfTemperatureSteps << "\n";
            continue;
          }
	   else if (key == "bestSequenceNumber") {
            lineStream >> bestSequenceNumber;
            //std::cout << "Simulated Annealing: bestSequenceNumber = " << bestSequenceNumber << "\n";
            continue;
          }
	  else if (key == "bestObservable0") {
            lineStream >> bestObservable0;
            //std::cout << "Simulated Annealing: bestObservable0 = " << bestObservable0 << "\n";
            continue;
          }

        }

      }

    }
    inputFile.close();

  }
  
  
  // Check consistency: thermalizationSteps
  if (thermalizationStepsPerformed >= numberOfThermalizationSteps) {
    std::cout << "\n   CAUTION! Thermalization steps performed from previous run >= numberOfThermalizationSteps in this run.";
    if (std::filesystem::exists("configurations/config_checkpoint.dat"))
      std::cout << "\n            - Configuration checkpoint file is found, thermalization will be skipped.\n\n";
    else {          // perform thermalization
      thermalizationStepsPerformed = 0;
      std::cout << "\n            - Configuration checkpoint file is not found, thermalization will be performed.\n\n";
    }    
  }

  // Check consistency: MCSteps
  if (MCStepsPerformed >= numberOfMCSteps) {
    std::cout << "\n   CAUTION! Number of MC steps performed from previous run >= numberOfMCSteps required.";
    std::cout << "\n            No further work will be performed. Quitting OWL...\n\n";
    exit(7);
  }

  // Check consistency: acceptedMoves and rejectedMoves
  assert (acceptedMoves + rejectedMoves == MCStepsPerformed * numberOfMCUpdatesPerStep);

  // Restore averagedObservables and averagedObservablesSquared for accumulation
  for (unsigned int i=0; i<physical_system->numObservables; i++) {
    standardDeviations[i] = standardErrors[i] * sqrt(double(MCStepsPerformed));
    averagedObservablesSquared[i] = (standardDeviations[i] * standardDeviations[i] + averagedObservables[i] * averagedObservables[i]) * double(MCStepsPerformed);
    averagedObservables[i] *= double(MCStepsPerformed);
  }
    
}

void SimulatedAnnealing::clearObservables()
{

  for (unsigned int i=0; i<physical_system->numObservables; i++) {
    averagedObservables[i]        = 0.0;;
    averagedObservablesSquared[i] = 0.0;
  }

}

void SimulatedAnnealing::accumulateObservables()
{

  for (unsigned int i=0; i<physical_system->numObservables; i++) {
    averagedObservables[i]        += physical_system -> observables[i];
    averagedObservablesSquared[i] += physical_system -> observables[i] * physical_system -> observables[i];
  }

}


void SimulatedAnnealing::calculateAveragesAndVariances()
{

  for (unsigned int i=0; i<physical_system->numObservables; i++) {
    averagedObservables[i]        /= double(numberOfMCSteps);
    averagedObservablesSquared[i] /= double(numberOfMCSteps);
    standardDeviations[i]          = sqrt( averagedObservablesSquared[i] - averagedObservables[i] * averagedObservables[i] );
    standardErrors[i]              = standardDeviations[i] / sqrt(double(numberOfMCSteps));
  }

}

void SimulatedAnnealing::writeMCFile(unsigned long int MCSteps)
{

  fprintf(timeSeriesFile, "%15lu ", MCSteps);
  
  for (unsigned int i=0; i<physical_system->numObservables; i++)
    fprintf(timeSeriesFile, "%15.6f ", physical_system -> observables[i]);

  fprintf(timeSeriesFile, "\n");

}


void SimulatedAnnealing::writeStatistics(OutputMode output_mode, const char* filename) 
{
  
  FILE* checkPointFile;
  if (filename != NULL)
    checkPointFile = fopen(filename, "w");
  else checkPointFile = stdout;

  switch(output_mode) {

    case endOfSimulation :

      fprintf(checkPointFile, "\n");
      fprintf(checkPointFile, "             Statistics of Simulated Annealing sampling \n");
      fprintf(checkPointFile, "   ----------------------------------------------------- \n");
      fprintf(checkPointFile, "   Simulation temperature         : %8.5f \n", temperature);
      fprintf(checkPointFile, "   Best energy found              : %f \n", bestObservable0);
      fprintf(checkPointFile, "   Number of thermalization steps :  %lu \n",  numberOfThermalizationSteps);
      fprintf(checkPointFile, "   Total number of MC steps       :  %lu \n",  numberOfMCSteps);
      fprintf(checkPointFile, "   Number of MC updates per step  :  %lu \n",  numberOfMCUpdatesPerStep);
      fprintf(checkPointFile, "   Number of accepted MC updates  :  %lu (%5.2f %%) \n", 
              acceptedMoves, double(acceptedMoves) / double(numberOfMCSteps * numberOfMCUpdatesPerStep) * 100.0);
      fprintf(checkPointFile, "   Number of rejected MC updates  :  %lu (%5.2f %%) \n", 
              rejectedMoves, double(rejectedMoves) / double(numberOfMCSteps * numberOfMCUpdatesPerStep) * 100.0);
      
      fprintf(checkPointFile, "\n");
    
      fprintf(checkPointFile, "                        Observable                          Mean          Std. error of the mean \n");
      fprintf(checkPointFile, "   ---------------------------------------------------------------------------------------------- \n");
      for (unsigned int i=0; i<physical_system -> numObservables; i++)
        fprintf(checkPointFile, "   %45s :     %12.5f         %12.5f \n", physical_system -> observableName[i].c_str(), averagedObservables[i], standardErrors[i]);    
      fprintf(checkPointFile, "\n"); 

      break;

    case checkPoint :

      fprintf(checkPointFile, "temperature                   %8.5f\n", temperature);
      fprintf(checkPointFile, "thermalizationStepsPerformed   %lu\n",  thermalizationStepsPerformed);
      fprintf(checkPointFile, "MCStepsPerformed               %lu\n",  MCStepsPerformed);
      fprintf(checkPointFile, "acceptedMoves                  %lu\n",  acceptedMoves);
      fprintf(checkPointFile, "rejectedMoves                  %lu\n",  rejectedMoves);

      fprintf(checkPointFile, "numberOfTemperatureSteps       %u\n",  numberOfTemperatureSteps);
      fprintf(checkPointFile, "bestSequenceNumber             %u\n",  bestSequenceNumber);
      fprintf(checkPointFile, "bestObservable0                %f\n",  bestObservable0);

      fprintf(checkPointFile, "averagedObservables   ");
      for (unsigned int i=0; i<physical_system -> numObservables; i++)
        fprintf(checkPointFile, "%12.5f      ", averagedObservables[i] / double(MCStepsPerformed));
      fprintf(checkPointFile, "\n");
      
      fprintf(checkPointFile, "standardErrors   ");
      for (unsigned int i=0; i<physical_system -> numObservables; i++) {
        double temp_ave       = averagedObservables[i] / double(MCStepsPerformed);
        double temp_ave2      = averagedObservablesSquared[i] / double(MCStepsPerformed);
        standardDeviations[i] = sqrt(temp_ave2 - temp_ave*temp_ave);
        standardErrors[i]     = standardDeviations[i] / sqrt(double(MCStepsPerformed));
        fprintf(checkPointFile, "%12.5f      ", standardErrors[i]);
      }
      fprintf(checkPointFile, "\n");
      
      break;

    default :
      break;
      
  }

  if (filename != NULL) fclose(checkPointFile);

}


void SimulatedAnnealing::writeCheckPointFiles(OutputMode output_mode)
{

  char fileName[51];

  switch (output_mode) {

    case endOfIteration :
      break;

    case endOfSimulation :
      sprintf(fileName, "configurations/config_final.dat");
      physical_system -> writeConfiguration(0, fileName);
      writeStatistics(endOfSimulation);
      writeStatistics(endOfSimulation, "simulated_annealing_final.dat");
      break;

    case checkPoint :
      physical_system -> writeConfiguration(0, "configurations/config_checkpoint.dat");
      writeStatistics(checkPoint, "simulated_annealing_checkpoint.dat"); 
      break;

    default :
      break;

  };

}

void SimulatedAnnealing::updateTemperature()
{
  switch(temeratureUpdateScheme) {
  default:
    temperature = 0.95 * temperature;
  };
}

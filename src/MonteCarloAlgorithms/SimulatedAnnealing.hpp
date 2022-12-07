#ifndef SIMULATED_ANNEALING_HPP
#define SIMULATED_ANNEALING_HPP

#include <fstream>
#include <vector>
#include <limits>
#include "MCAlgorithms.hpp"


class SimulatedAnnealing : public MonteCarloAlgorithm {

public :

  SimulatedAnnealing(PhysicalSystem* ps, const char* inputFile);
  ~SimulatedAnnealing();

  void run()                  override;

private :

  PhysicalSystem* physical_system;

  ObservableType bestObservable0;
  unsigned int bestSequenceNumber {0};
  int temeratureUpdateScheme {0};
  unsigned int numberOfTemperatureSteps {1};

  std::vector<ObservableType> averagedObservables;
  std::vector<ObservableType> averagedObservablesSquared;
  std::vector<ObservableType> standardDeviations;
  std::vector<ObservableType> standardErrors;

  FILE* timeSeriesFile;

  unsigned long int numberOfThermalizationSteps;
  unsigned long int numberOfMCSteps;
  unsigned long int numberOfMCUpdatesPerStep {1};
  double temperature;

  unsigned long int thermalizationStepsPerformed {0};
  unsigned long int MCStepsPerformed             {0};
  unsigned long int totalMCStepsPerformed        {0};

  void readMCInputFile(const char* fileName);  // TODO: this should move to MCAlgorithms base class (Histogram class has the same function)
  void readCheckPointFile(const char* fileName);

  void clearObservables();
  void accumulateObservables();                // TODO: this should move to MCAlgorithms base class
  void calculateAveragesAndVariances();
  
  void writeMCFile(unsigned long int MCSteps);
  void writeStatistics(OutputMode output_mode, const char* = NULL);    // TODO: this should move to MCAlgorithms base class
  void writeCheckPointFiles(OutputMode output_mode);                   // TODO: this should move to MCAlgorithms base class

  void updateTemperature();
};

#endif

// -----------------------------------------------------------------
// Run Monte Carlo calculation using mix of Metropolis
// and Wolff algorithms and periodic boundary conditions
// Print energy, average phi, specific heat, susceptibility,
// autocorrelation times, errors, correlations, etc.
#include <math.h>                   // For floor and sqrt
#include <stdio.h>                  // For printf and fprintf
#include <gsl/gsl_sf_log.h>         // For natural log
#include <gsl/gsl_math.h>           // For power
#include <gsl/gsl_fft_real.h>       // Fast Fourier Transform
#include "Lattice.hh"               // Method and data declarations
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Declare and initialize data
unsigned int randomSite = 0;

double aveEnergy        = 0;
double avePhi           = 0;
double avePhiAbs        = 0;
double squaredEnergy    = 0;
double squaredPhi       = 0;
double quartPhi         = 0;
double cumulant         = 0;
double specificHeat     = 0;
double susceptibility   = 0;

// Bimodality binning stuff
double maxPhi                   = 0;
static const unsigned int bins  = 21;
unsigned int counts[bins];
double midCounts                = 0;            // Make these double now
double maxCounts                = 0;            // to avoid casts later
double bimodality               = 0;

double scaleFactor  = 0;
double autocorTime  = 0;
double energyStDev  = 0;
double phiStDev     = 0;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Calculate autocorrelation time and necessary points of
// autocorrelation function
unsigned int calcAutocor(unsigned int sampleSize, double* phiDataAbs,
                         double* autocorrelation) {

  unsigned int i, t = 1;
  double temp;

  // Generate Chi[0] for scaling purposes
  for (i = 0; i < sampleSize; i++)
    scaleFactor += (phiDataAbs[i] * phiDataAbs[i]);
  scaleFactor /= sampleSize;
  scaleFactor -= avePhiAbs * avePhiAbs;

  autocorrelation[0] = 1;

  // Only calculate points of autocorrelation function that
  // roughly conform to exponential approximation
  while (t < sampleSize) {
    autocorrelation[t] = 0;
    for (i = 0; i < sampleSize - t; i++)
      autocorrelation[t] += (phiDataAbs[i] * phiDataAbs[i + t]);

    autocorrelation[t] /= (sampleSize - t);
    autocorrelation[t] -= avePhiAbs * avePhiAbs;
    if (autocorrelation[t] < 0)     // Not exponential!
      break;

    autocorrelation[t] /= scaleFactor;      // Scale by Chi[0]
    if (autocorrelation[t] >= autocorrelation[t - 1])   // Not exponential!
      break;

    temp = -gsl_sf_log(autocorrelation[t]);
    temp = t / temp;
    autocorTime += temp;

    t++;
  }
  return t;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Bin values of phi and calculate bimodality
// Values of phi range from -maxPhi to +maxPhi
// Bin i of n will contain values greater than maxPhi(2*i/n - 1)
// and less than maxPhi*(2*(i+1)/n - 1)  (where i = 0...n - 1)
void calcBimodality(unsigned int sampleSize, double* phiData) {
  unsigned int i, j;
  double lowerBound, upperBound;

  for (i = 0; i < bins; i++)
    counts[i] = 0;

  for (i = 0; i < sampleSize; i++)
    for (j = 0; j < bins; j++) {
      lowerBound = maxPhi * (2 * (double)j / (double)bins - 1);
      upperBound = maxPhi * (2 * ((double)j + 1) / (double)bins - 1);
      if (phiData[i] >= lowerBound && phiData[i] <= upperBound) {
        counts[j]++;
        break;
      }
    }

  midCounts = counts[(bins - 1) / 2];
  for (i = 0; i < bins; i++)
    if (counts[i] > maxCounts)
      maxCounts = counts[i];

  bimodality = 1 - (midCounts / maxCounts);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Main method runs simulation using command line parameters
int main(int argc, char** const argv) {
  unsigned int i, j, k, t, gap = 5;

  if (argc != 6) {
    fprintf(stderr, "Usage: %s muSquared lambda length ", argv[0]);
    fprintf(stderr, "init sampleSize\n");
    fflush(stderr);
    exit(1);
  }

  // Data that depends on command-line parameters
  double muSquared = atof(argv[1]) / 10000;   // Bare mass squared
  double lambda = atof(argv[2]) / 100;        // Bare coupling
  unsigned int length = atoi(argv[3]);        // Lattice length (square)
  unsigned int latticeSize = length * length;
  unsigned int init = atoi(argv[4]);          // Iterations for equilibration
  unsigned int sampleSize = atoi(argv[5]);    // Iterations for statistics

  // Number of uncorrelated measurements -- samplesize divided by
  // autocorrelation time
  double measurements = 0;

  double energyData[sampleSize];
  double phiData[sampleSize];
  double phiDataAbs[sampleSize];
  // Autocorrelation function -- not all array elements will be used
  double autocorrelation[sampleSize];

  // Position-space and momentum-space correlation functions
  unsigned int halfL = int(length / 2);
  double posCorr[halfL];
  double momCorr[halfL];
  for (i = 0; i < halfL; i++)
    posCorr[i] = 0;

  Lattice* theLattice = new Lattice(muSquared, lambda, length);

  // Initialize/equilibrate lattice
  // Do cluster update after every _gap_ Metropolis sweeps
  gap = 5;
  for (i = 0; i < init; i++) {
    for (j = 0; j < gap; j++) {
      for (k = 0; k < latticeSize; k++)
        theLattice->metropolis(k);
    }

    randomSite = (unsigned int)floor(latticeSize *
                                gsl_rng_uniform(theLattice->generator));
    theLattice->wolff(randomSite);
  }

  // Sweeps with measurements turned on
  // Do cluster update after every _gap_ Metropolis sweeps
  for (i = 0; i < sampleSize; i++) {
    for (j = 0; j < gap; j++) {
      for (k = 0; k < latticeSize; k++)
        theLattice->metropolis(k);
    }

    randomSite = (unsigned int)floor(latticeSize *
                                gsl_rng_uniform(theLattice->generator));
    theLattice->wolff(randomSite);

    energyData[i] = theLattice->calcTotalEnergy();
    aveEnergy += energyData[i];

    phiData[i] = theLattice->calcAveragePhi();
    phiDataAbs[i] = fabs(phiData[i]);
    avePhi += phiData[i];
    avePhiAbs += phiDataAbs[i];

    // Maximum (absolute) value of phi
    if (phiDataAbs[i] > maxPhi)
      maxPhi = phiDataAbs[i];

    squaredEnergy += energyData[i] * energyData[i];
    squaredPhi += phiData[i] * phiData[i];
    quartPhi += gsl_pow_4(phiData[i]);

    // Calculate spatial correlations
    theLattice->calcCorrelations(posCorr);
  }

  // Done with this
  delete theLattice;

  // Take averages
  aveEnergy     /= sampleSize;
  avePhi        /= sampleSize;
  avePhiAbs     /= sampleSize;
  squaredEnergy /= sampleSize;
  squaredPhi    /= sampleSize;
  quartPhi      /= sampleSize;
  for (i = 0; i < halfL; i++) {
    posCorr[i] /= sampleSize;
    posCorr[i] -= avePhi * avePhi;
  }

  // Add bootstrapping here if desired
  specificHeat = squaredEnergy - (aveEnergy * aveEnergy);
  specificHeat *= latticeSize;
  susceptibility = squaredPhi - (avePhiAbs * avePhiAbs);
  susceptibility *= latticeSize;

  // Estimate autocorrelation and standard deviations
  // Use average phi for autocorrelation time calculation
  // Should be roughly the same for all variables
  t = calcAutocor(sampleSize, phiDataAbs, autocorrelation);

  // Use standard formula to generate standard deviations
  // from autocorrelation time, average, sampleSize, etc
  if (t == 1) {            // Autocorrelation time zero
    autocorTime = 0;
    energyStDev = 0;
    phiStDev    = 0;
  }
  else {
    autocorTime /= (double)(t - 1.0);
    measurements = sampleSize / autocorTime;

    energyStDev = 2.0 * autocorTime / sampleSize;
    energyStDev *= squaredEnergy - (aveEnergy * aveEnergy);
    energyStDev = sqrt(energyStDev);

    phiStDev = 2.0 * autocorTime / sampleSize;
    phiStDev *= squaredPhi - (avePhiAbs * avePhiAbs);
    phiStDev = sqrt(phiStDev);
  }

  // Calculate binder cumulant
  cumulant = 1.0 - quartPhi / (3.0 * squaredPhi * squaredPhi);

  // Calculate bimodality
  calcBimodality(sampleSize, phiData);

  // Take Fourier transform of correlation functions and
  // extrapolate to zero momentum
  for (i = 0; i < halfL; i++)
    momCorr[i] = posCorr[i];
  gsl_fft_real_radix2_transform(momCorr, 1, halfL);

  printf("%lf,%lf,", muSquared, lambda);
  printf("%lf,%lf,", autocorTime, measurements);
  printf("%lf,%lf,", aveEnergy, energyStDev);
  printf("%lf,%lf,", avePhiAbs, phiStDev);
  printf("%lf,%lf,", specificHeat, susceptibility);
  printf("%lf,%lf,", cumulant, bimodality);
  printf("%lf,%lf,", avePhi, scaleFactor);

  // Position-space functions only have real components
  for (i = 0; i < halfL; i++)
    printf("%lf,", posCorr[i]);

  // Momentum-space functions are complex in general; print real
  // and imaginary parts together (half-complex arrangement)
  printf("%lf,%lf,", momCorr[0], 0.0);
  for (i = 1; i < (unsigned int)(halfL / 2); i++)
    printf("%lf,%lf,", momCorr[i], momCorr[halfL - i]);
  printf("%lf,%lf,", momCorr[(unsigned int)(halfL / 2)], 0.0);
  for (i = (unsigned int)(halfL / 2) + 1; i < (unsigned int)(halfL); i++)
    printf("%lf,%lf,", momCorr[halfL - i], -momCorr[i]);
  printf("\n");

  return 0;
}
// -----------------------------------------------------------------
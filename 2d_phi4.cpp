// -----------------------------------------------------------------
// Run Monte Carlo calculation using mix of Metropolis
// and Wolff algorithms and periodic boundary conditions
// Print energy, average phi, specific heat, suscept,
// autocorrelation times, errors, correlations, etc.
#include <math.h>                   // For floor and sqrt
#include <stdio.h>                  // For printf and fprintf
#include <ctime>                    // For clock
#include <gsl/gsl_sf_log.h>         // For natural log
#include <gsl/gsl_math.h>           // For power
#include "Lattice.hh"               // Method and data declarations
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Declare and initialize global data
unsigned int randomSite = 0;

double aveEnergy = 0.0, squaredEnergy = 0.0, energyStDev = 0.0;
double avePhi = 0.0, avePhiAbs = 0.0, squaredPhi = 0.0, phiStDev = 0.0;
double quartPhi = 0.0, cumulant = 0.0;
double specHeat = 0.0, suscept = 0.0;

// Bimodality binning stuff
// Make (integer) counts doubles to avoid casting later on
static const unsigned int bins  = 21;
unsigned int counts[bins];
double midCounts = 0.0, maxCounts = 0.0;
double maxPhi = 0.0, bimod = 0.0;

// Autocorrelation stuff -- scaleFactor is Chi[0]
double scaleFactor = 0.0, autocorTime = 0.0;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Timing helper routines
// Print time stamp
void time_stamp(char *msg) {
  time_t time_stamp;

  time(&time_stamp);
  printf("%s: %s\n", msg, ctime(&time_stamp));
  fflush(stdout);
}

// Double precision CPU time in seconds
double dclock() {
  return (((double)clock()) / CLOCKS_PER_SEC);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Calculate autocorrelation time of |phi|
// and necessary points of autocorrelation function
unsigned int calcAutocor(unsigned int meas, double *phiDataAbs,
                         double *autocor) {

  unsigned int i, t = 1;
  double td;

  // Generate Chi[0] for scaling purposes
  for (i = 0; i < meas; i++)
    scaleFactor += (phiDataAbs[i] * phiDataAbs[i]);
  scaleFactor /= meas;
  scaleFactor -= avePhiAbs * avePhiAbs;

  autocor[0] = 1;

  // Only calculate points of autocor function that
  // roughly conform to exponential approximation
  while (t < meas) {
    autocor[t] = 0;
    for (i = 0; i < meas - t; i++)
      autocor[t] += (phiDataAbs[i] * phiDataAbs[i + t]);

    autocor[t] /= (meas - t);
    autocor[t] -= avePhiAbs * avePhiAbs;
    if (autocor[t] < 0)     // Not exponential!
      break;

    autocor[t] /= scaleFactor;      // Scale by Chi[0]
    if (autocor[t] >= autocor[t - 1])   // Not exponential!
      break;

    td = -gsl_sf_log(autocor[t]);
    td = t / td;
    autocorTime += td;

    t++;
  }
  return t;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Bin values of phi and calculate bimod
// Values of phi range from -maxPhi to +maxPhi
// Bin i of n will contain values greater than maxPhi(2*i/n - 1)
// and less than maxPhi*(2*(i+1)/n - 1)  (where i = 0...n - 1)
void calcBimodality(unsigned int meas, double *phiData) {
  unsigned int i, j;
  double lowerBound, upperBound;

  for (i = 0; i < bins; i++)
    counts[i] = 0;

  for (i = 0; i < meas; i++)
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

  bimod = 1 - (midCounts / maxCounts);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Main method runs simulation using command line parameters
int main(int argc, const char **argv) {
  unsigned int i, j, k, t, gap = 5;
  double td, td2;

  if (argc != 6) {
    fprintf(stderr, "Usage: %s muSquared lambda L init meas\n", argv[0]);
    fflush(stderr);
    exit(1);
  }

  // Load and print data from command-line parameters
  double muSquared = atof(argv[1]);           // Bare mass squared
  double lambda = atof(argv[2]);              // Bare coupling
  unsigned int L = atoi(argv[3]);             // Length of (square) lattice
  unsigned int latticeSize = L * L;
  unsigned int init = atoi(argv[4]);          // Iterations for equilibration
  unsigned int meas = atoi(argv[5]);    // Iterations for statistics

  printf("Two-dimensional phi^4 theory\n");
  printf("Metropolis with Wolff cluster flip every %d sweeps\n", gap);
  time_stamp("start");
  printf("muSquared %.g\n", muSquared);
  printf("lambda %.g\n", lambda);
  printf("L %d\n", L);
  printf("init %d\n", init);
  printf("meas %d\n\n", meas);

  // Other local data
  double energyData[meas];
  double phiData[meas], phiDataAbs[meas];
  double dtime = -dclock(), wtime;

  // Number of uncorrelated measurements -- samplesize divided by
  // autocorrelation time
  double measurements = 0.0;

  // Autocorrelation function -- not all array elements will be used
  double autocor[meas];

  // Position-space and momentum-space correlation functions
  unsigned int halfL = int(L / 2);
  double posCorr[halfL], momCorr[halfL];
  for (i = 0; i < halfL; i++) {
    posCorr[i] = 0;
    momCorr[i] = 0;
  }

  Lattice *theLattice = new Lattice(muSquared, lambda, L);
  td = theLattice->calcTotalEnergy();
  td2 = theLattice->calcAveragePhi();
  printf("START %.8g %.8g\n", td, td2);

  // Initialize/equilibrate lattice
  // Do cluster update after every _gap_ Metropolis sweeps
  wtime = -dclock();
  for (i = 0; i < init; i++) {
    for (j = 0; j < gap; j++) {
      for (k = 0; k < latticeSize; k++)
        theLattice->metropolis(k);
    }

    randomSite = (unsigned int)floor(latticeSize *
                                gsl_rng_uniform(theLattice->generator));
    theLattice->wolff(randomSite);
  }
  wtime += dclock();
  printf("%d WARMUPS COMPLETED in %.4g seconds\n", init, wtime);

  // Sweeps with measurements turned on
  // Do cluster update after every _gap_ Metropolis sweeps
  for (i = 0; i < meas; i++) {
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

    // Print energy and magnetization
    printf("MEAS %.6g %.6g\n", energyData[i], phiData[i]);

    // Calculate two-point function in both position and momentum space
    // Only print (real) position-space correlator to reduce output size
    theLattice->calcCorrelations(posCorr, momCorr, phiData[i]);
  }

  // Done with this
  delete theLattice;

  // Take averages
  aveEnergy     /= meas;
  avePhi        /= meas;
  avePhiAbs     /= meas;
  squaredEnergy /= meas;
  squaredPhi    /= meas;
  quartPhi      /= meas;
  for (i = 0; i < halfL; i++) {
    posCorr[i] /= meas;
    momCorr[i] /= meas;
  }

  // TODO: If we're going to print these here
  //       (rather than analyzing them offline with some other script)
  //       Then we should include jackknife or bootstrap uncertainties...
  specHeat = squaredEnergy - (aveEnergy * aveEnergy);
  specHeat *= latticeSize;
  suscept = squaredPhi - (avePhiAbs * avePhiAbs);
  suscept *= latticeSize;

  // Estimate autocorrelation and standard deviations
  // Use average phi for autocorrelation time calculation
  // Should be roughly the same for all variables
  t = calcAutocor(meas, phiDataAbs, autocor);

  // Use standard formula to generate standard deviations
  // from autocorrelation time, average, meas, etc
  if (t == 1) {            // Autocorrelation time zero
    autocorTime = 0.0;
    energyStDev = 0.0;
    phiStDev    = 0.0;
  }
  else {
    autocorTime /= (double)(t - 1.0);
    measurements = meas / autocorTime;

    energyStDev = 2.0 * autocorTime / meas;
    energyStDev *= squaredEnergy - (aveEnergy * aveEnergy);
    energyStDev = sqrt(energyStDev);

    phiStDev = 2.0 * autocorTime / meas;
    phiStDev *= squaredPhi - (avePhiAbs * avePhiAbs);
    phiStDev = sqrt(phiStDev);
  }

  // Calculate binder cumulant
  cumulant = 1.0 - quartPhi / (3.0 * squaredPhi * squaredPhi);

  // Calculate bimod
  calcBimodality(meas, phiData);

  // Print autocorrelation information, raw averages and derived quantities
  printf("\nAUTOCOR %.6g %.6g %.6g\n", autocorTime, measurements, scaleFactor);
  printf("AVE %.6g %.6g %.6g %.6g %.6g\n",
         aveEnergy, energyStDev, avePhiAbs, phiStDev, avePhi);
  printf("DERIVED %.6g %.6g %.6g %.6g\n", specHeat, suscept, cumulant, bimod);

  // Position-space functions only have real components
  printf("AVE_CORR");
  for (i = 0; i < halfL; i++)
    printf(" %.6g", posCorr[i]);
  printf("\n");

  // Momentum-space functions are half-complex
  // Rearrange half-complex format to print real and imaginary parts together
  printf("AVE_MOM_CORR %.6g 0.0", momCorr[0]);
  for (i = 1; i < (unsigned int)(halfL / 2); i++)
    printf(" %.6g %.6g", momCorr[i], momCorr[halfL - i]);
  printf(" %.6g 0.0", momCorr[(unsigned int)(halfL / 2)]);
  for (i = (unsigned int)(halfL / 2) + 1; i < (unsigned int)(halfL); i++)
    printf(" %.6g %.6g", momCorr[halfL - i], -momCorr[i]);
  printf("\n");

  dtime += dclock();
  printf("\nTime = %.4g seconds\n", dtime);
  fflush(stdout);
  return 0;
}
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Lattice for phi^4 simulations using mu action
#ifndef _LATTICE_HH
#define _LATTICE_HH

#include "HashTable.hh"             // Hash table for searching cluster
#include <vector>                   // Lattice is vector of vectors
#include <gsl/gsl_rng.h>            // Random number generators
#include <gsl/gsl_sf_exp.h>         // Exponential functions

// A simple struct to hold a site's neighbors
struct siteNeighbors {
  unsigned int prevX;
  unsigned int nextX;
  unsigned int prevY;
  unsigned int nextY;
};
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Class definition
class Lattice {
  public:
    // Member data
    std::vector<double> lattice;    // Continuous values
    unsigned int length;            // Lattice length (square)
    unsigned int latticeSize;       // Number of sites in lattice

    double muSquared;               // Bare mass-squared
    double lambda;                  // Bare coupling

    // Neighboring lattice sites
    std::vector<siteNeighbors*> neighbors;

    HashTable *cluster;
    gsl_rng *generator;

    // Constructors, destructor
    Lattice(double m, double l, unsigned int length);
    Lattice();
    ~Lattice();

    // Set up periodic boundary conditions
    void getNeighbors(unsigned int site, siteNeighbors* toInit);

    // Calculation methods
    double calcTotalEnergy();
    double calcAveragePhi();
    void calcCorrelations(double spatialCorr[], double momCorr[],
                          double phibar);

    // Metropolis and Wolff algorithms
    void metropolis(unsigned int site);

    bool clusterCheck(unsigned int site, unsigned int toAdd);

    // Inelegant but faster
    void growClusterPos(unsigned int site);
    void growClusterNeg(unsigned int site);
    void flipCluster();
    unsigned int wolff(unsigned int site);      // Returns cluster size
};

#endif
// -----------------------------------------------------------------

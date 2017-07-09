// -----------------------------------------------------------------
// Corr/Lattice_FFT.hh
// Lattice of spins for phi^4 simulations using mu action
// Header file contains data and method declarations
// David Schaich -- daschaich@gmail.com
// Created 12 October 2005
// Last modified 18 April 2007
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Frontmatter and include directives
// Avoid multiple inclusion
#ifndef _LATTICE_HH
#define _LATTICE_HH

#include "HashTable.hh"             // Hash table for searching cluster
#include <vector>                   // Lattice is vector of vectors
#include "include/gsl_rng.h"        // Random number generators
#include "include/gsl_sf_exp.h"     // Exponential functions
// -----------------------------------------------------------------



// -----------------------------------------------------------------
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
        // ---------------------------------------------------------
        // Data!
        // Member data
        std::vector<double> lattice;    // Continuous values
        unsigned int length;            // Lattice length (square)
        unsigned int latticeSize;       // Number of sites in lattice

        double muSquared;       // Mass of particles
        double lambda;          // Coupling strength

        // Neighboring lattice sites
        std::vector<siteNeighbors*> neighbors;

        HashTable* cluster;
        gsl_rng* generator;
        // ---------------------------------------------------------



        // ---------------------------------------------------------
        // Methods!
        // Constructors, destructor
        Lattice(double m, double l, unsigned int length);
        Lattice();
        ~Lattice();

        // Set up periodic boundary conditions
        void getNeighbors(unsigned int site, siteNeighbors* toInit);

        // Calculation methods
        double calcTotalEnergy();
        double calcAveragePhi();
        void calcCorrelations(double posCorr[], double momCorr[]);

        // Simulation methods - metropolis and wolff algorithms
        void metropolis(unsigned int site);

        bool clusterCheck(unsigned int site, unsigned int toAdd);

        // Inelegant but faster
        void growClusterPos(unsigned int site);
        void growClusterNeg(unsigned int site);
        void flipCluster();
        unsigned int wolff(unsigned int site);      // Returns cluster size
        // ---------------------------------------------------------
};
// -----------------------------------------------------------------



// -----------------------------------------------------------------
#endif  // _LATTICE_HH
// -----------------------------------------------------------------


// -----------------------------------------------------------------
// Corr/Lattice.cpp
// Lattice of spins for phi^4 simulations using periodic b.c.
// Contains implementations of standard methods
// David Schaich -- daschaich@gmail.com
// Created 12 October 2005
// Last modified 23 April 2007
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Frontmatter and include directives
#include "Lattice.hh"               // Method and variable declarations
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Constructors and destructor
Lattice::Lattice(double m, double l, unsigned int x) {
    generator = gsl_rng_alloc(gsl_rng_mt19937);     // Mersenne Twister
    gsl_rng_set (generator, (unsigned int)(100 * m * l));

    // Simplify some expressions below
    muSquared = 2 + (m / 2);
    lambda = l / 4;

    length = x;
    latticeSize = length * length;
    cluster = new HashTable(latticeSize / 4);

    // Random initial state in range [-1.5, 1.5)
    for (unsigned int i = 0; i < latticeSize; i++)
        lattice.push_back(3 * gsl_rng_uniform(generator) - 1.5);

    // Set up neighbors... calculate once for all sites
    for (unsigned int i = 0; i < latticeSize; i++) {
        siteNeighbors* temp = new siteNeighbors;
        getNeighbors(i, temp);
        neighbors.push_back(temp);
    }
}

Lattice::Lattice()  {
    Lattice(-1.25, 1, 32);
}

Lattice::~Lattice() {}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set up periodic boundary conditions - only do once
void Lattice::getNeighbors(unsigned int site, siteNeighbors* toInit) {
    if ((site + 1) % length == 0)
        toInit->nextX = site + 1 - length;
    else
        toInit->nextX = site + 1;

    if (site >= latticeSize - length)
        toInit->nextY = site + length - latticeSize;
    else
        toInit->nextY = site + length;

    if (site % length == 0)
        toInit->prevX = site + length - 1;
    else
        toInit->prevX = site - 1;

    if (site < length)
        toInit->prevY = site + latticeSize - length;
    else
        toInit->prevY = site - length;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Calculation methods
// Calculate total energy by looping through lattice
double Lattice::calcTotalEnergy() {
    double totalEnergy = 0;
    double currentPhi;
    for (unsigned int i = 0; i < latticeSize; i++) {
        currentPhi = lattice[i];
        totalEnergy -= currentPhi * (lattice[neighbors[i]->nextX]
                                  +  lattice[neighbors[i]->nextY]);

        currentPhi *= currentPhi;
        // Recall muSquared redefined
        totalEnergy += muSquared * currentPhi;

        currentPhi *= currentPhi;
        // Recall lambda redefined
        totalEnergy += lambda * currentPhi;
    }
    return totalEnergy / latticeSize;
}

// Note: does not return absolute value
double Lattice::calcAveragePhi() {
    double currentPhi = 0;
    for (unsigned int i = 0; i < latticeSize; i++)
        currentPhi += lattice[i];

    return currentPhi / latticeSize;
}

// Calculates connected two-point spatial correlation functions
void Lattice::calcCorrelations(double spatialCorr[]) {
    unsigned int halfL = int(length / 2);
    double corr[halfL];             // To be added to input array
    unsigned int cur = 0;           // Current site

    for (unsigned int i = 0; i < halfL; i++)
        corr[i] = 0;

    // Loop root of paths over lattice to get full average
    for (unsigned int i = 0; i < latticeSize; i++) {
        // Loop over various correlation functions
        // One for each path length
        for (unsigned int j = 0; j < halfL; j++) {
            cur = i;

            // Move to right-most point to be added
            for (unsigned int k = 0; k < j; k++)
                cur = neighbors[cur]->nextX;

            corr[j] += lattice[i] * lattice[cur];

            // Shift up and left
            for (unsigned int k = 0; k < j; k++) {
                cur = neighbors[cur]->prevX;
                cur = neighbors[cur]->nextY;
                corr[j] += lattice[i] * lattice[cur];
            }

            // Shift down and left
            for (unsigned int k = 0; k < j; k++) {
                cur = neighbors[cur]->prevX;
                cur = neighbors[cur]->prevY;
                corr[j] += lattice[i] * lattice[cur];
            }

            // Shift down and right
            for (unsigned int k = 0; k < j; k++) {
                cur = neighbors[cur]->nextX;
                cur = neighbors[cur]->prevY;
                corr[j] += lattice[i] * lattice[cur];
            }

            // Shift up and right -- careful not to double-count
            // Also be careful not to wrap around unsigned ints
            if (j == 0)
                continue;
            for (unsigned int k = 0; k < j - 1; k++) {
                cur = neighbors[cur]->nextX;
                cur = neighbors[cur]->nextY;
                corr[j] += lattice[i] * lattice[cur];
            }
        }   // End loop over correlation functions
    }   // End loop over root sites

    // Add results to previously-accumulated results
    spatialCorr[0] += corr[0] / latticeSize;
    for (unsigned int i = 1; i < halfL; i++)
        spatialCorr[i] += corr[i] / (4 * i * latticeSize);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Metropolis method
void Lattice::metropolis(unsigned int site) {
    double currentPhi = lattice[site];

    // Generate new value
    double newValue = currentPhi + (3 * gsl_rng_uniform(generator) - 1.5);
    double temp = newValue;

    // Calculate energy difference
    double difference = (currentPhi - newValue)
                      * (lattice[neighbors[site]->nextX]
                      +  lattice[neighbors[site]->nextY]
                      +  lattice[neighbors[site]->prevX]
                      +  lattice[neighbors[site]->prevY]);

    newValue *= newValue;
    currentPhi  *= currentPhi;
    // Recall muSquared redefined
    difference += muSquared * (newValue - currentPhi);

    newValue *= newValue;
    currentPhi *= currentPhi;
    // Recall lambda redefined
    difference += lambda * (newValue - currentPhi);

    // Flip if difference <= 0, otherwise probabilistic acceptance
    if (difference <= 0)
        lattice[site] = temp;
    else if (gsl_rng_uniform(generator) < gsl_sf_exp(-difference))
        lattice[site] = temp;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Wolff methods for growing cluster and so on
// A convenience method that keeps me from having to write these few
// lines over and over again - adds to cluster probabilistically
bool Lattice::clusterCheck(unsigned int site, unsigned int toAdd) {
    if (cluster->find(toAdd))
        return false;

    // Could try calculating these once for each cluster...
    double probability = 1 - gsl_sf_exp(-2 * lattice[site] * lattice[toAdd]);

    if (gsl_rng_uniform(generator) < probability) {
        cluster->insert(toAdd);
        return true;
    }

    return false;
}

// Grows cluster from specified site - recursive
void Lattice::growClusterPos(unsigned int site) {
    unsigned int toCheck = neighbors[site]->prevX;
    if (lattice[toCheck] > 0 && clusterCheck(site, toCheck))
        growClusterPos(toCheck);

    toCheck = neighbors[site]->nextX;
    if (lattice[toCheck] > 0 && clusterCheck(site, toCheck))
        growClusterPos(toCheck);

    toCheck = neighbors[site]->prevY;
    if (lattice[toCheck] > 0 && clusterCheck(site, toCheck))
        growClusterPos(toCheck);

    toCheck = neighbors[site]->nextY;
    if (lattice[toCheck] > 0 && clusterCheck(site, toCheck))
        growClusterPos(toCheck);
}

void Lattice::growClusterNeg(unsigned int site) {
    unsigned int toCheck = neighbors[site]->prevX;
    if (lattice[toCheck] <= 0 && clusterCheck(site, toCheck))
        growClusterNeg(toCheck);

    toCheck = neighbors[site]->nextX;
    if (lattice[toCheck] <= 0 && clusterCheck(site, toCheck))
        growClusterNeg(toCheck);

    toCheck = neighbors[site]->prevY;
    if (lattice[toCheck] <= 0 && clusterCheck(site, toCheck))
        growClusterNeg(toCheck);

    toCheck = neighbors[site]->nextY;
    if (lattice[toCheck] <= 0 && clusterCheck(site, toCheck))
        growClusterNeg(toCheck);
}

// Since this method trawls through the whole cluster (which has reached
// the end of its usefulness), let's be efficient and clear it here
void Lattice::flipCluster() {
    node* temp1;
    node* temp2;

    for (unsigned int i = 0; i < cluster->tableNumber; i++) {
        temp1 = cluster->table[i];
        temp2 = cluster->table[i];
        while (temp1 != NULL) {
            lattice[temp1->value] *= -1;
            temp2 = temp2->next;
            delete temp1;
            temp1 = temp2;
        }
        cluster->table[i] = NULL;
    }
    cluster->size = 0;
}

unsigned int Lattice::wolff(unsigned int site) {
    cluster->insert(site);

    // Inelegant, but faster than having just one method
    if (lattice[site] > 0)
        growClusterPos(site);
    else
        growClusterNeg(site);

    unsigned int toReturn = cluster->size;
    flipCluster();

    return toReturn;
}
// -----------------------------------------------------------------


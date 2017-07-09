// -----------------------------------------------------------------
// Corr/Lattice_4d.cpp
// Lattice of spins for continuous simulations
// Contains implementations of standard methods
// David Schaich -- daschaich@gmail.com
// Created 4 November 2005
// Last modified 30 April 2006
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Include directives
#include "Lattice_4d.hh"        // Method and variable declarations
#include "include/gsl_math.h"   // For power
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Constructors and destructor
Lattice::Lattice(double m, double l, unsigned int x) {
    generator = gsl_rng_alloc(gsl_rng_mt19937);     // Mersenne Twister
    gsl_rng_set (generator, (unsigned int)(100 * m * l));

    muSquared = 4 + (m / 2);        // Redefine to speed up calcs
    lambda = l / 4;                 // Redefine to speed up calcs

    length = x;
    latticeSize = int(gsl_pow_4(x));
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
    Lattice(-1.25, 1, 8);
}

Lattice::~Lattice() {}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Set up periodic boundary conditions -- only do once
// Try to avoid the modulo % operation
void Lattice::getNeighbors(unsigned int site, siteNeighbors* toInit) {
    if ((site + 1) % length == 0)
        toInit->nextX = site + 1 - length;
    else
        toInit->nextX = site + 1;

    if (site % length == 0)
        toInit->prevX = site + length - 1;
    else
        toInit->prevX = site -1;

    int length2 = length * length;
    if ((int(site / length) + 1 ) % length == 0)
        toInit->nextY = site + length - length2;
    else
        toInit->nextY = site + length;

    if (int(site / length) % length == 0)
        toInit->prevY = site + length2 - length;
    else
        toInit->prevY = site - length;

    int length3 = length2 * length;
    if ((int(site / length2) + 1) % length == 0)
        toInit->nextZ = site + length2 - length3;
    else
        toInit->nextZ = site + length2;

    if(int(site / length2) % length == 0)
        toInit->prevZ = site + length3 - length2;
    else
        toInit->prevZ = site - length2;

    if((int(site / length3) + 1) % length == 0)
        toInit->nextCT = site + length3 - latticeSize;
    else
        toInit->nextCT = site + length3;

    if (int(site / length3) % length == 0)
        toInit->prevCT = site + latticeSize - length3;
    else
        toInit->prevCT = site - length3;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Calculation methods
// Calculate total energy, looping through lattice
double Lattice::calcTotalEnergy() {
    double totalEnergy = 0;
    double currentPhi;

    for (unsigned int i = 0; i < latticeSize; i++) {
        currentPhi = lattice[i];
        totalEnergy -= currentPhi
                    *  (lattice[neighbors[i]->nextX]
                    +   lattice[neighbors[i]->nextY]
                    +   lattice[neighbors[i]->nextZ]
                    +   lattice[neighbors[i]->nextCT]);

        currentPhi *= currentPhi;
        totalEnergy += muSquared * currentPhi;

        currentPhi *= currentPhi;
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
void Lattice::calcCorrelations(double posCorr[]) {
    unsigned int halfL = int(length / 2);
    unsigned int l2 = length * length;
    unsigned int l3 = l2 * length;
    double corr[halfL];             // To be added to input array
    unsigned int count[halfL];      // For normalization
    unsigned int x, y, z, ct;       // Root site
    unsigned int dx, dy, dz, dct;   // Site separations

    for (unsigned int i = 0; i < halfL; i++) {
        corr[i] = 0;
        count[i] = 0;
    }

    for (unsigned int i = 0; i < latticeSize; i++) {
        x  = i % length;
        y  = abs(int(i / length)) % length;
        z  = abs(int(i / l2)) % length;
        ct = abs(int(i / l3));

        for (unsigned int j = 0; j < latticeSize; j++) {
            dx  = abs(x  - (j % length));
            dy  = abs(y  - (int(j / length) % length));
            dz  = abs(z  - (int(i / l2) % length));
            dct = abs(ct - int(i / l3));
            if (dx  >= length / 2)
                dx   = length - dx;
            if (dy  >= length / 2)
                dy   = length - dy;
            if (dz  >= length / 2)
                dz   = length - dz;
            if (dct >= length / 2)
                dct  = length - dct;

            if (dx + dy + dz + dct < halfL) {
                corr[dx + dy + dz + dct] += lattice[i] * lattice[j];
                count[dx + dy + dz + dct]++;
            }
        }
    }

    for (unsigned int i = 0; i < halfL; i++)
        posCorr[i] += corr[i] / count[i];
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Metropolis method
void Lattice::metropolis(unsigned int site) {
    double currentPhi = lattice[site];

    // Generate new value
    double newValue = lattice[site] + (3 * gsl_rng_uniform(generator) - 1.5);
    double temp = newValue;

    // Calculate energy difference
    double difference = (currentPhi - newValue)
                      * (lattice[neighbors[site]->nextX]
                      +  lattice[neighbors[site]->prevX]
                      +  lattice[neighbors[site]->nextY]
                      +  lattice[neighbors[site]->prevY]
                      +  lattice[neighbors[site]->nextZ]
                      +  lattice[neighbors[site]->prevZ]
                      +  lattice[neighbors[site]->nextCT]
                      +  lattice[neighbors[site]->prevCT]);

    currentPhi *= currentPhi;
    newValue *= newValue;
    difference += muSquared * (newValue - currentPhi);

    currentPhi *= currentPhi;
    newValue *= newValue;
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

    toCheck = neighbors[site]->prevZ;
    if (lattice[toCheck] > 0 && clusterCheck(site, toCheck))
        growClusterPos(toCheck);

    toCheck = neighbors[site]->nextZ;
    if (lattice[toCheck] > 0 && clusterCheck(site, toCheck))
        growClusterPos(toCheck);

    toCheck = neighbors[site]->prevCT;
    if (lattice[toCheck] > 0 && clusterCheck(site, toCheck))
        growClusterPos(toCheck);

    toCheck = neighbors[site]->nextCT;
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

    toCheck = neighbors[site]->prevZ;
    if (lattice[toCheck] <= 0 && clusterCheck(site, toCheck))
        growClusterNeg(toCheck);

    toCheck = neighbors[site]->nextZ;
    if (lattice[toCheck] <= 0 && clusterCheck(site, toCheck))
        growClusterNeg(toCheck);

    toCheck = neighbors[site]->prevCT;
    if (lattice[toCheck] <= 0 && clusterCheck(site, toCheck))
        growClusterNeg(toCheck);

    toCheck = neighbors[site]->nextCT;
    if (lattice[toCheck] <= 0 && clusterCheck(site, toCheck))
        growClusterNeg(toCheck);
}

// Looping through cluster to flip it, might as well clear it at the
// same time
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


// -----------------------------------------------------------------
// Lattice for phi^4 with periodic BC
// Contains implementations of standard methods

#include <gsl/gsl_fft_real.h>       // Fast Fourier Transform
#include "Lattice.hh"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Constructors and destructor
Lattice::Lattice(double m, double l, unsigned int x) {
  unsigned int i;

  generator = gsl_rng_alloc(gsl_rng_mt19937);     // Mersenne Twister
  gsl_rng_set (generator, (unsigned int)(100 * m * l));

  // Simplify some expressions below
  muSquared = 2 + (m / 2);
  lambda = l / 4;

  length = x;
  latticeSize = length * length;
  cluster = new HashTable(latticeSize / 4);

  // Random initial state in range [-1.5, 1.5)
  for (i = 0; i < latticeSize; i++)
    lattice.push_back(3 * gsl_rng_uniform(generator) - 1.5);

  // Set up neighbors... calculate once for all sites
  for (i = 0; i < latticeSize; i++) {
    siteNeighbors *temp = new siteNeighbors;
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
// Set up periodic boundary conditions
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
  unsigned int i;
  double currentPhi, totalEnergy = 0;

  for (i = 0; i < latticeSize; i++) {
    currentPhi = lattice[i];
    totalEnergy -= currentPhi * (lattice[neighbors[i]->nextX]
                 + lattice[neighbors[i]->nextY]);

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
  unsigned int i;
  double currentPhi = 0;
  for (i = 0; i < latticeSize; i++)
    currentPhi += lattice[i];

  return currentPhi / latticeSize;
}

// Calculate and print connected two-point spatial correlation function
// Use 'Manhattan distance' -- probably not optimal...
// Subtract the lattice average phibar
// Add to input arrays rather than overwriting them
void Lattice::calcCorrelations(double posCorr[], double momCorr[],
                               double phibar) {

  unsigned int i, j, halfL = int(length / 2);
  double corr[halfL];         // To be added to input arrays
  unsigned int count[halfL];  // For normalization
  unsigned int x, y;          // Root site
  unsigned int dx, dy;        // Site separations

  for (i = 0; i < halfL; i++) {
    corr[i] = 0;
    count[i] = 0;
  }

  for (i = 0; i < latticeSize; i++) {
    x = i % length;
    y = abs(int(i / length));
    for (j = 0; j < latticeSize; j++) {
      dx = abs(x - (j % length));
      dy = abs(y - int(j / length));
      if (dx >= halfL)
        dx = length - dx;
      if (dy >= halfL)
        dy = length - dy;

      if (dx + dy < halfL) {
        corr[dx + dy] += (lattice[i] - phibar) * (lattice[j] - phibar);
        count[dx + dy]++;
      }
    }
  }

  // Normalize, accumulate and print position-space two-point function
  printf("CORR");
  for (i = 0; i < halfL; i++) {
    corr[i] /= count[i];
    posCorr[i] += corr[i];
    printf(" %.6g", corr[i]);
  }
  printf("\n");

  // Compute and accumulate momentum-space two-point function
  // Don't print it to reduce output size
  // It can be reconstructed from printed CORR data above
  gsl_fft_real_radix2_transform(corr, 1, halfL);    // In-place FFT
  for (i = 0; i < halfL; i++)
    momCorr[i] += corr[i];

  // Rearrange half-complex format to print real and imaginary parts together
//  printf("MOM_CORR %.6g 0.0", corr[0]);
//  for (i = 1; i < (unsigned int)(halfL / 2); i++)
//    printf(" %.6g %.6g", corr[i], corr[halfL - i]);
//  printf(" %.6g 0.0", corr[(unsigned int)(halfL / 2)]);
//  for (i = (unsigned int)(halfL / 2) + 1; i < (unsigned int)(halfL); i++)
//    printf(" %.6g %.6g", corr[halfL - i], -corr[i]);
//  printf("\n");
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
                       + lattice[neighbors[site]->nextY]
                       + lattice[neighbors[site]->prevX]
                       + lattice[neighbors[site]->prevY]);

  newValue *= newValue;
  currentPhi *= currentPhi;
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
// Add to cluster probabilistically (convenience method)
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

// Grow cluster from specified site - recursive
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
  unsigned int i;
  node *temp, *temp2;

  for (i = 0; i < cluster->tableNumber; i++) {
    temp = cluster->table[i];
    temp2 = cluster->table[i];
    while (temp != NULL) {
      lattice[temp->value] *= -1;
      temp2 = temp2->next;
      delete temp;
      temp = temp2;
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

  flipCluster();
  return cluster->size;
}
// -----------------------------------------------------------------


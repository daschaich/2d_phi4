// -----------------------------------------------------------------
// Basic hash table as a chained list of nodes
#ifndef _HASHTABLE_HH
#define _HASHTABLE_HH

#include <vector>
#include <cstddef>    // Now needed for NULL in C++11...
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Struct of a node in a slot in the hashtable
// Seems to reduce object-oriented overhead compared to a class
struct node {
  unsigned int value;
  node* next;
};

// Class declaration for HashTable
class HashTable {
  public:
    // Constructor, destructor
    HashTable(unsigned int tableNumber);
    HashTable();
    ~HashTable();

    // Member functions
    void insert(unsigned int site);
    bool find(unsigned int site);

    // Member data
    unsigned int size;
    unsigned int tableNumber;
    unsigned int mod;
    std::vector<node*> table;       // Vector of lists of nodes
};

#endif
// -----------------------------------------------------------------

// -----------------------------------------------------------------
// Corr/HashTable.hh
// Basic hash table with chaining
// David Schaich -- daschaich@gmail.com
// Created 23 October 2005
// Last modified 25 January 2006
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Avoid multiple inclusion in a single executable
#ifndef _HASHTABLE_HH
#define _HASHTABLE_HH
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Include directives
#include <vector>
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Structs instead of classes - reduce OOP overhead
// Struct of a node in a slot in the hashtable
struct node {
    unsigned int value;
    node* next;
};
// -----------------------------------------------------------------



// -----------------------------------------------------------------
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
// -----------------------------------------------------------------



// -----------------------------------------------------------------
#endif // _HASHTABLE_HH
// -----------------------------------------------------------------


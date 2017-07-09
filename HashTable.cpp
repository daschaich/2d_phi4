// -----------------------------------------------------------------
// Corr/HashTable.cpp
// Hash table designed for lists of nodes
// David Schaich
// Created 23 October 2005
// Last modified 25 January 2006
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Include directives
#include "HashTable.hh"
#include <vector>
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Constructors and destructor
HashTable::HashTable(unsigned int numberOfTables) {
    size = 0;
    tableNumber = numberOfTables;
    std::vector<node*>* temp = new std::vector<node*>(tableNumber, NULL);
    table = *temp;

    mod = tableNumber - 1;
}

HashTable::HashTable() {
    HashTable(4093);
}

HashTable::~HashTable() {}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Member functions - insertion, searching
void HashTable::insert(unsigned int site) {
    size++;
    unsigned int index = (17 * site - 97) & mod;

    node* toAdd = new node;
    toAdd->value = site;

    toAdd->next = table[index];
    table[index] = toAdd;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// This returns whether or not site is in cluster
bool HashTable::find(unsigned int site) {
    unsigned int index = (17 * site - 97) & mod;
    node* temp = table[index];

    while (temp != NULL) {
        if (site == temp->value)
            return true;
        temp = temp->next;
    }
    return false;
}
// -----------------------------------------------------------------


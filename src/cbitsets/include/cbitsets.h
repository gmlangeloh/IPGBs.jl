#ifndef CBITSETS_H
#define CBITSETS_H

#include <stdlib.h>
#include <stdbool.h>

typedef unsigned long WORD;

//Bitset as an array of WORDs
typedef struct {
    WORD *words;
    size_t size;
} bitset;

bitset* bitset_create(unsigned num_bits);
void bitset_set(bitset* set, unsigned index);
void bitset_fill(bitset* set, int* indices, unsigned num_indices);
bool is_disjoint(bitset *a, bitset *b);
bool not_disjoint(bitset *a, bitset *b);

#endif
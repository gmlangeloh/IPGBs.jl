#include <math.h>
#include "cbitsets.h"

bitset* bitset_create(unsigned num_bits){
    bitset* set = malloc(sizeof(bitset));
    unsigned bits_per_word = sizeof(WORD) * 8;
    set->size = ceil((double)num_bits / bits_per_word);
    set->words = calloc(set->size, sizeof(WORD));
    return set;
}

void bitset_set(bitset* set, unsigned index){
    unsigned bits_per_word = sizeof(WORD) * 8;
    unsigned word_index = index / bits_per_word;
    unsigned bit_index = index % bits_per_word;
    set->words[word_index] |= (1UL << bit_index);
}

void bitset_fill(bitset* set, int* indices, unsigned num_indices){
    for(unsigned i = 0; i < num_indices; i++){
        bitset_set(set, indices[i]);
    }
}

bool is_disjoint(bitset *a, bitset *b){
    for(unsigned i = 0; i < a->size; i++){
        if(a->words[i] & b->words[i]){
            return false;
        }
    }
    return true;
}

bool not_disjoint(bitset *a, bitset *b){
    for(unsigned i = 0; i < a->size; i++){
        if(a->words[i] & b->words[i]){
            return true;
        }
    }
    return false;
}

//int main() {
//    bitset* set = bitset_create(100);
//    int indices[] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
//    int num_indices = 9;
//    bitset_fill(set, indices, num_indices);
//    bitset* set2 = bitset_create(100);
//    int indices2[] = {10, 11, 12};
//    int num_indices2 = 3;
//    bitset_fill(set2, indices2, num_indices2);
//    printf("%d\n", is_disjoint(set, set2));
//    return 0;
//}
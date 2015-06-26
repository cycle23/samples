#include <iostream>
#include <vector>

// Works in place
// O n log n

// input: array of n numbers, unsorted
// output: same numbers, sorted in increasing order
// assume: all array entries distinct

// extra - handle duplicate entries

// partitioning around a pivot

// simple case: just use first item as pivot
// then it is about comparing to this duude, left/right
// during partition, can be randomly ordered on left or right

// 1) O(n) can be performed in linear time, no extra memory, repeated in-place swaps
// 2) reduces problem size

using namespace std;

int choosePivot(int *a, int l, int r) {
// runtime if sorted and we always choose first item as pivot:
//  n + (n-1) + (n-2) + ... + 1
// O(n^2)

// best pivot would be ideally in the middle of our item
// the median element

// if in every call we get a perfect 50/50 split, it's theta-(n log n)
// but that is unlikely

// instead, choose a pivot randomly every time

// if we always get a 25/75 split, that's good enough for O(n log n)
// half of the elements give that 25-75, so 50% of the time this is what we expect

    return l+ (rand() % (r - l));
}
// assumes p == 0 for now
// time here is actually theta(n), but can say O(n)
// where n = r - l + 1
// invariant:
//    at every step, all a[l+1]..a[i-1] are all < p
//    all a[i]..a[j-1] are > p

int pivot(int *a,int k, int l, int r) {
    // if not 0, swap pivot with 1st element as preprocessing step
    if (k != l)
        swap(a[l],a[k]);
    int p = a[l];
    int i = l;
    for (int j = l + 1; j < r; j++) {
        if (a[j] < p) {
            swap(a[++i],a[j]);
        }
    }
    swap(a[l],a[i]);
    return i;
    // single scan through array
    // invariant: i = the pivot point so far, j = the stuff not seen yet
    // 3 8 2 5 1 4 7 6
    //   ij
    // step 1: i = j = 1

    // 3 8 2 5 1 4 7 6
    //   i j
    // step 2: j++; i = 1; j = 2;

    // 3 8 2 5 1 4 7 6
    //   i   j
    // step 3: j++; i = 1; j = 3;

    // 3 2 8 5 1 4 7 6
    //     i j
    // step 4: 2 < pivot, swapped with 8, i++; i = 2; j = 3;

    // 3 2 8 5 1 4 7 6
    //     i   j
    // step 5: j++, i = 2; j = 4;

    // 3 2 5 8 1 4 7 6
    //     i   j
    // step 6: 5 > pivot, swapped with 8 (or not), i = 2; j = 4;

    // 3 2 1 8 5 4 7 6
    //       i   j
    // step 7: 1 < pivot, swapped, i++, i = 3, j = 5 (we can stop here if we assume unique..)

    // two optional steps:
    // 3 2 1 8 5 4 7 6
    //       i     j
    // step 8: j++, i = 3, j=6 (4 > pivot)

    // 3 2 1 8 5 4 7 6
    //       i       j
    // step 9: j++, i = 3, j = 6

    // 1 2 3 8 5 4 7 6
    // step 10: swap pivot with i - 1, done

}

void quickSort(int *a, int l, int r) {
    if (r - l <= 1)
        return;
    int p = choosePivot(a,l,r);
    int k = pivot(a,p,l,r);
    quickSort(a, l, k);
    quickSort(a, k, r);
}
void qSort(int *a , int n) {
    quickSort(a, 0, n);
}

// choosing good pivot
// randomized quickSort
// Analysis:
// decomposition principle
//  - linearity of expectation (hashing)
// key insight of quicksort
//  - probability of comparison
// bound of quicksort running time

int main() {
    int a[12] = {4,3,2,1,6,7,9,10,23,12,54,11};
    qSort(a,12);
    return 0;
}


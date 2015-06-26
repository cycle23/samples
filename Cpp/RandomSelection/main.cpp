#include <iostream>
#include <iostream>
#include <vector>

// ordered statistics of an array including median
// with running time linear O(n) using randomized selection

// ith smallest entry in the array

// input: array A with n distinct numbers
// and a number i in {1..n}

// output ith order statistic
// example: median

// i = (n + 1 / 2)
// i = ( n / 2 for n even)

// naive approach
// sort O (n log n)
// return ith element

// prove selection is faster than sorting

// selecting pivot is O(n) using randomization in quicksort
// O(n) time deterministic algorithm
// - pivot deterministically using median of medians
// - not as practical (ie: larger space.. imma guess recursion)

// this one is in place, for instance

// exactly.. partition around the fucking pivot. yo.
// log n calls log n depth like n^log n = n.

// if n = 1 return a[0]
// choose pivot p from A uniformly at random
// and partition now for i / j (need median but keep choosing min or max)
// let j = new index of p      (meaning n/2 calls, and each time by a constant = n^2)
// if i = j return p
// if j > i then RSelect(a, j-1, i)
// if j < i then RSelect(a + j, n-j, i-j)

// runtime analysis is due to that same 25-75% range issue
// if we assumed median, it is just one recursive call
// T(n) <= T(n/2) + O(n)
// hope: random pivot is "pretty good" "often enough"

int choosePivot(int *a, int l, int r) {
    return l+ (rand() % (r - l));
}

void selectRand(int *a, int n, int i) {
    if (r - l <= 1)
        return;
    int p = choosePivot(a,l,r);
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
    selectRand(a, l, k);
}


void selectDet(int *a, int n, int i) {
    // break A into groups of 5, sort each group
    // C = the n/5 "middle elements"
    // p = selectDet(C, n/5, n/10) // recursively computes median of C
    // partition A around p
    // if j = i return p
    // if j < i return selectDet(a, j-1, i)
    // if j > i return selectDet(A+j, n -j, i -j)
}


using namespace std;

int main() {
    cout << "Hello, World!" << endl;
    return 0;
}
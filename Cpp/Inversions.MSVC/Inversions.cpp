// Inversions.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <iostream>
#include <vector>
#include <fstream>


// sortAndCount (array A, length n)
// if n = 1 return 0
// else
// (b,x) = sortAndCount(1st half of A, n/2)
// (c,y) = sortAndCount(2nd half of A, n/2)
// (d,z) = mergeAndCountSplitInv(A, n)
// return x + y + z

// Goal: countSplitInv: linear(O(n))
// => count = O(n log n)

// If entire left half >n/2, and entire right half <n/2, then are quadratic number
// of possible inversions

// countSplitInv - similar to mergeSort

// have recursive calls both count inversions and to sort
// sorting reveals split inversions


using namespace std;

pair<vector<long>, long> mergeAndCountSplitInv(vector<long> left, vector<long> right) {
	vector<long> merged;
	long i = 0, j = 0;
	long inverted = 0;
	while (i < left.size() || j < right.size()) {
		if (i < left.size()) {
			if (j < right.size()) {
				if (left[i] < right[j]) {
					merged.push_back(left[i++]);
				}
				else {
					inverted += left.size() - i;
					merged.push_back(right[j++]);
				}
			}
			else {
				merged.push_back(left[i++]);
			}
		}
		else {
			merged.push_back(right[j++]);
		}
	}
	return make_pair(merged, inverted);
}

pair<vector<long>, long> sortAndCount(vector<long> a) {
	if (a.size() <= 1) return make_pair(a, 0);
	else if (a.size() == 2) {
		if (a[0] < a[1]) {
			return make_pair(a, 0);
		}
		else {
			vector<long> ret({ a[1], a[0] });
			return make_pair(ret, 1);
		}
	}
	else {
		size_t half = a.size() / 2;
		pair<vector<long>, long> Bx = sortAndCount(vector<long>(a.begin(), a.begin() + half));
		pair<vector<long>, long> Cy = sortAndCount(vector<long>(a.begin() + half, a.end()));
		pair<vector<long>, long> Dz = mergeAndCountSplitInv(Bx.first, Cy.first);
		return make_pair(Dz.first, Bx.second + Cy.second + Dz.second);
	}
}


int _tmain(int argc, _TCHAR* argv[])
{
	ifstream infile("IntegerArray.txt");
	vector<long> A;
	/*
	while (!infile.bad() && !infile.eof()){
		long next;
		infile >> next;
		A.push_back(next);
	}
	pair<vector<long>, long> res = sortAndCount(A);
	*/
	pair<vector<long>, long> res = sortAndCount({ 1, 2, 3, 4, 6, 5 });
	cout << res.second << endl;
	return 0;
}


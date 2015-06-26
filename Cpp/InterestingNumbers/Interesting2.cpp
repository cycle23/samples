// Interesting2.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <vector>
using namespace std;

//class InterestingDigits {
//public:
	vector<int> digits(int base) {
		vector<int> interesting;
		// base 3 to 30
		// from 2 to base - 1
		for (int i = 2; i < base; i++) {
			// 4 + 4 = 8 = 1 3 = 4
			bool isInteresting = true;
			for (int j = 1; j < 100; j++) {
				int mult = i * j;
				int sum = 0;
				int digits = 0;
				while (mult) {
					digits++;
					sum += mult % base;
					mult /= base;
				}
				if (sum % base) {
					isInteresting = false;
					break;
				}
				// we went past 3 digits, we can stop
				if (digits > 3) {
					break;
				}
			}
			if (isInteresting) {
				interesting.push_back(i);
			}
		}
		return interesting;
	}
//};

int _tmain(int argc, _TCHAR* argv[])
{
	vector<int> inter = digits(5);
	return 0;
}


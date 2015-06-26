// MathSamples.cpp
// https://www.topcoder.com/community/data-science/data-science-tutorials/mathematics-for-topcoders/

#include "stdafx.h"
#include <math.h>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <iostream>

using namespace std;

int sqrt(int n) {
	return sqrt(static_cast<double>(n));
}

class primes {
private:
	// naive isPrime
	bool isPrimeNaive(int n) {
		for (int i = 2; i < n; i++)
			if (n % i == 0) return false;
		return true;
	}
	// some optimization
	bool isPrimeBetter(int n) {
		if (n <= 1) return false;
		if (n == 2) return true;
		if (n % 2 == 0) return false;
		int m = sqrt(n);
		for (int i = 3; i <= m; i += 2)
			if (n%i == 0)
				return false;
		return true;
	}
	// Sieve of Eratosthenes
	vector<bool> sieve(int n) {
		vector<bool> prime(n+1,true);
		prime[0] = false;
		prime[1] = false;
		int m = sqrt(n);
		for (int i = 2; i <= m; i++)
			if (prime[i])
				for (int k = i*i; k <= n; k += i)
					prime[k] = false;
		return prime;
	}
	bool isPrimeSieve(int n) {
		vector<bool> s = sieve(n);
		return s[n];
	}

public:
	bool isPrime(int n) {
		return isPrimeSieve(n);
	}
};

class GCD {
private:
	int gcdNaive(int a, int b) {
		for (int i = min(a, b); i >= 1; i--)
			if (a%i ==0 && b%i ==0)
				return i;
	}
	int gcdEuclid(int a, int b) {
		if (b == 0) return a;
		return gcdEuclid(b, a%b);
	}
	// note: similar recursive method with large number, small number, remainder, recursion = 
	// linear equation, ie: ax + by = c
public:
	int gcd(int a, int b) {
		return gcdEuclid(a, b);
	}
	int lcm(int a, int b) {
		return b*a / gcd(a, b);
	}
};

class Rectangle {
private:
	int xBottomLeft, yBottomLeft, xTopRight, yTopRight;
public:
	bool intersects(Rectangle& rect) {

		// does not intersect
		//             ----
		//             |  |
		//   ----      ----
		//   |  |       
		//   ----
		if (max(xBottomLeft, rect.xBottomLeft) > min(xTopRight,rect.xTopRight) ||
			max(yBottomLeft, rect.yBottomLeft) > min(yTopRight, rect.yTopRight))
		return false;
		else return true;
	}
};

// picks theorem
// B = # of lattice points on boundary
// I = # of lattice points on interior
// Area = B/2 + I - 1?
// huh? 

// euler's formula for polygonal nets
// 
//   _____
//  / F| F\edge
//  |__|__|vertice   Face = all internal polygons and the external polygon itself
//  |  |  |
//  \_F|_F/
// V - E + F = 2
// If V=2, then there are E faces (E-1 in the middle and 1 outside)
// V - E + F = 2 - E + E = 2
// Let V = n+1
// choose vertex w and assume it is connected by G edges
// if we remove w and those edges, then net has n vertices, E-G edges, F-G+1 faces
// so we can assume that:
// (n) - (E - G) + (F - G + 1) = 2
// (n + 1) - E + F = 2

class Bases {
public:
	int toDecimal(int n, int b)
	{
		int result = 0;
		int multiplier = 1;
		while (n > 0) {
			result += n % 10 * multiplier;
			multiplier *= b;
			n /= 10;
		}
		return result;
	}
	// 2<=b<=20
	string fromDecimal(int n, int b)
	{
		const string chars = "0123456789ABCDEFGHIJ";
		ostringstream result;
		while (n > 0) {
			result << chars[n%b];
			n /= b;
		}
		string s = result.str();
		reverse(s.begin(),s.end());
		return s;
	}

};

class Fraction {
	vector<int> f;
	GCD g;
public:
	vector<int> multiply(vector<int> &a, vector<int> &b) {
		return vector<int>(a[0] * b[0], a[1] * b[1]);
	}
	vector<int> add(vector<int> &a, vector<int> &b) {
		int denom = g.lcm(a[1], b[1]);
		return vector<int>(denom/a[1]*a[0] + denom/b[1]*b[0], denom);
	}
	vector<int> reduce(vector<int> &a) {
		int b = g.gcd(a[0], a[1]);
		return vector<int>(a[0] /= b, a[1] /= b);
	}
};

class Complex {
	// m = a+ib
	// n = c+id
	// m + n = (a + ib) + (c + id)
	//       = (a + c) + i(b + d)

	// m * n = (a + ib) * (c + id)
	//       =  ac + iad + ibc + (i^2)bd
	//       = (ac - bd) + i(ad + bc)
public:
	vector<int> multiplyComplex(vector<int> &m, vector<int> n) {
		return vector<int>(m[0] * n[0] - m[1] * n[1], m[0] * n[1] + m[1] * n[0]);
	}
};

int _tmain(int argc, _TCHAR* argv[])
{
	Bases b;
	cout << b.fromDecimal(100, 16);
	return 0;
}


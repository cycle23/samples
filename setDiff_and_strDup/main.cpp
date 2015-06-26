#include "setDiff.h"
#include "stripDup.h"
#include <iostream>
using namespace std;

class testSetDiff {
public:
    static bool pass() {
        vector<int> A{1,3,5};
        vector<int> B{5,7,9};
        vector<int> expected{1,3};
        auto diff = setDiff::setDifference(A,B);
        //cout << diff << endl;
        return (diff == expected);
    }
};

class testStripDup {
public:
    static bool pass() {
        vector<string> dups{"A","A","Abracadabra","Banana","B","A","B","Abracadabra","Z","Z","z","A","Z","..."};
        vector<string> expected{"A","Abracadabra","Banana","B","Z","z","..."};
        auto uniq = stripDup::stripDuplicates(dups);
        //cout << uniq << endl;
        return (uniq == expected);
    }
};


int main() {
    cout << "test setDifference:" << testSetDiff::pass() << endl;
    cout << "test stripDup:" << testStripDup::pass() << endl;
    return 0;
}
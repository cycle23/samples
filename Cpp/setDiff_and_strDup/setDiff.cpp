//
// Created by Cody on 6/11/2015.
//

#include "setDiff.h"
#include <unordered_map>

vector<int> setDiff::setDifference(const vector<int> &A, const vector<int> &B) {
    vector<int> diff;
    unordered_map<int, bool> bMap;
    for (auto b: B) {
        bMap[b] = true;
    }
    for (auto a: A) {
        if (bMap.find(a)==end(bMap)) {
            diff.push_back(a);
        }
    }
    return move(diff);
}

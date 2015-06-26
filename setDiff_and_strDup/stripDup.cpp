//
// Created by Cody on 6/11/2015.
//

#include "stripDup.h"
#include <unordered_map>

vector<string> stripDup::stripDuplicates(const vector<string> &A) {
    vector<string> uniq;
    unordered_map<string, bool> seen;
    for (auto s : A) {
        if (seen.find(s) == end(seen)) {
            seen[s] = true;
            uniq.push_back(s);
        }
    }
    return move(uniq);
}
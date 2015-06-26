#include <iostream>
#include <vector>
#include <limits>

// https://www.topcoder.com/community/data-science/data-science-tutorials/dynamic-programming-from-novice-to-advanced/
using namespace std;

int main(int argc, char **argv) {
    // min: will keep track of minimum coins to get to each sub-sum
    vector<int> min;
    // value of each type of coin
    vector<int> coins;
    coins.resize(3);
    coins[0] = 1;
    coins[1] = 3;
    coins[2] = 5;
    // target sum
    int sum = 11;
    int steps = 0;
    int corrections = 0;
    min.resize(sum);
    for (int i = 0; i <= sum; i++) {
        min[i] = numeric_limits<int>::max();
    }
    // assume it takes 0 coins for sum 1
    min[0] = 0;
    for (int i = 1; i <= sum; i++) {
        for (int j = 0; j < coins.size(); j++) {
            // value of next coin to consider
            int Vj = coins[j];
            // if this value fits within current sub-sum 'i'
            // and adding one coin to the value of current sum
            // minus the value of the current coin is less
            // than the current minimum coin count for this sum,
            // then count the minimum for this sum accordingly
            if (Vj <= i && min[i - Vj] + 1 < min[i]) {
                min[i] = min[i - Vj] + 1;
                corrections++;
            }
            steps++;
        }
    }
    cout << "number: " << min[sum] << endl;
    cout << "steps: " << steps << endl;
    cout << "corrections: " << corrections << endl;
    steps = 0;
    corrections = 0;
    // another method:
    for (int i = 1; i <= sum; i++) {
        min[i] = numeric_limits<int>::max();
    }
    for (int j = 0; j < coins.size(); j++) {
        int Vj = coins[j];
        for (int i = Vj; i <= sum; i ++) {
            steps++;
            if (min[i - Vj] + 1 < min[i]) {
                min[i] = min[i - Vj] + 1;
                corrections++;
            }
        }
    }
    cout << "number: " << min[sum] << endl;
    cout << "steps: " << steps << endl;
    cout << "corrections: " << corrections << endl;

    int seqInts[] = {5, 3, 4, 8, 6, 7, 2};
    unsigned int seqSize = sizeof(seqInts) / sizeof(int);
    vector<int> A(seqInts, seqInts + seqSize);
    // assume that longest consecutive increasing subseq so far to this is 1
    // (subseq can have holes, we want longest increasing sequence skipping non-increasing)
    vector<int> S(seqSize,1);
    // A   S
    // 5   1
    // 3   1
    // 4   2
    // 8   3 (from 5 to here = 2, but from 4 to here = 3)
    // 6   3 (skip 8, then from 4 to here = 3)
    // 7   4
    // 2   1
    // from the 2nd number on..
    for (int i = 1; i < seqSize; i++) {
        // scan back so far to this point
        for (int j = 0; j < i; j++) {
            // if the current entry is greater than that prior entry
            // and the count of going from that prior entry to here is greater
            // than the current assumed count to here so far
            if ((A[i] > A[j]) && (S[i] < S[j] + 1)) {
                // then set the current count to here up by one
                S[i] = S[j]+1;
            }
        }
    }
    // we'd need to rescan to get the max, or keep it in another val so far..
    cout << S[seqSize-1] << endl;

    {

        // graph = 2 dimensional array representing possibility of going from one edge to another
        // between singly linked vertices
        vector<vector<int> > graph(4, vector<int>(4));
        //                           * 0
        //               <-(1)-------/ |  \---------   (2)->
        //       /----------/ ^        |            \--------\
        //      * 1       (8)/      ^  |(5)                   * 2
        //         -----            |  | |          (1)<-   -/
        //              \ <-(3)   (10) | v         /-------/
        //               ------        |   -------/(4)->
        //                     \------ * 3 /
        graph[0][0] = -1;
        graph[0][1] = -1;
        graph[0][2] = 2;
        graph[0][3] = 5;

        graph[1][0] = 8;
        graph[1][1] = -1;
        graph[1][2] = -1;
        graph[1][3] = 1;

        graph[2][0] = -1;
        graph[2][1] = 1;
        graph[2][2] = -1;
        graph[2][3] = 4;

        graph[3][0] = 10;
        graph[3][1] = 3;
        graph[3][2] = 4;
        graph[3][3] = -1;

        // step 1, start at vertex 0
        // step 2, look for the next valid vertex with shortest path to it and not yet taken
        // step 3, mark the last vertex index taken to here there now
        // step 4, repeat for the next set of vertices, until vertex = N or all unmarked vertex combinations were searched
        vector<vector<int>> seen(4, vector<int>(4, -1));
        struct used {
            int weight;
            int from;
        };
        vector<used> weight(4,{-1,-1});
        weight[0] = {0,0};
        for (int i=0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                if (graph[i][j] == -1) continue;
                if (weight[j].weight == -1)
                    weight[j] = {graph[i][j], i};
                else
                    if (weight[i].weight >0 && weight[j].weight > (weight[i].weight + graph[i][j])) {
                        weight[j] = {weight[i].weight + graph[i][j], i};
                    }
            }
        }
        if (weight[3].weight == -1) {
            cout << "no path" << endl;
        }
        else {
            cout << "score: " << weight[3].weight << endl;
            int i = 3;
            cout << "vert:" << i << endl;
            while (i != 0) {
                i = weight[i].from;
                cout << "vert:" << i << endl;
            }
        }
    }
    return 0;
}
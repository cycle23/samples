// Graphs.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <list>
#include <memory>
#include <vector>
#include <stack>
#include <queue>
#include <functional>
#include <chrono>
#include <limits>

bool AreSame(double a, double b) {
	return std::fabs(a - b) < std::numeric_limits<double>::epsilon();
}
bool combinedToleranceCompare(double x, double y)
{
	double maxXYOne = std::max({ 1.0, std::fabs(x), std::fabs(y) });

	return std::fabs(x - y) <= std::numeric_limits<double>::epsilon()*maxXYOne;
}

using namespace std;


/*
* https://kartikkukreja.wordpress.com/2013/05/28/indexed-min-priority-queue-c-implementation/
* https://github.com/kartikkukreja/blog-codes/blob/master/src/Indexed%20Min%20Priority%20Queue.cpp
* Indexed min priority queue
* Supports insertion in O(log N), deletion of any key (regardless of whether
* the key is the minimum key or not) in O(log N) and changes to key values
* in O(log N), where N is the number of
* elements currently in the PQ
*/
template <class T>
class MinIndexedPQ {
	int NMAX, N;
	vector<int> heap; 
	vector<int> index;
	vector<T> keys;

	void swap(int i, int j) {
		auto t = heap[i]; heap[i] = heap[j]; heap[j] = t;
		index[heap[i]] = i; index[heap[j]] = j;
	}

	void bubbleUp(int k)    {
		while (k > 1 && keys[heap[k / 2]] > keys[heap[k]])   {
			swap(k, k / 2);
			k = k / 2;
		}
	}

	void bubbleDown(int k)  {
		auto j = k << 1;
		while (2 * k <= N) {
			j = k << 1;
			if (j < N && keys[heap[j]] > keys[heap[j + 1]])
				j++;
			if (keys[heap[k]] <= keys[heap[j]])
				break;
			swap(k, j);
			k = j;
		}
	}

public:
	// Create an empty MinIndexedPQ which can contain atmost NMAX elements
	MinIndexedPQ(int NMAX)  {
		this->NMAX = NMAX;
		N = 0;
		keys.resize(NMAX + 1);
		heap.resize(NMAX + 1);
		index.resize(NMAX + 1);
		for (auto i = 0; i <= NMAX; i++)
			index[i] = -1;
	}

	// due to C++11 rules, we should make sure compiler uses default items if we define
	// destructor
	~MinIndexedPQ() {
	}

	MinIndexedPQ(const MinIndexedPQ&) = default;
	MinIndexedPQ& operator=(const MinIndexedPQ&) = default;

	// moving, particularly useful for our private vectors
	MinIndexedPQ(MinIndexedPQ&&) = default;
	MinIndexedPQ& operator=(MinIndexedPQ&&) = default;

	// check if the PQ is empty
	bool isEmpty()  {
		return N == 0;
	}

	// check if i is an index on the PQ
	bool contains(int i)    {
		return index[i] != -1;
	}

	// return the number of elements in the PQ
	int size()  {
		return N;
	}

	// associate key with index i; 0 < i < NMAX
	void insert(int i, T key) {
		N++;
		index[i] = N;
		heap[N] = i;
		keys[i] = key;
		bubbleUp(N);
	}

	// returns the index associated with the minimal key
	int minIndex()  {
		return heap[1];
	}

	// returns the minimal key
	T minKey()    {
		return keys[heap[1]];
	}

	// delete the minimal key and return its associated index
	// Warning: Don't try to read from this index after calling this function
	int deleteMin() {
		auto min = heap[1];
		swap(1, N--);
		bubbleDown(1);
		index[min] = -1;
		heap[N + 1] = -1;
		return min;
	}

	// returns the key associated with index i
	T keyOf(int i)    {
		return keys[i];
	}

	// change the key associated with index i to the specified value
	void changeKey(int i, T key)  {
		keys[i] = key;
		bubbleUp(index[i]);
		bubbleDown(index[i]);
	}

	// decrease the key associated with index i to the specified value
	void decreaseKey(int i, T key)    {
		keys[i] = key;
		bubbleUp(index[i]);
	}

	// increase the key associated with index i to the specified value
	void increaseKey(int i, T key)    {
		keys[i] = key;
		bubbleDown(index[i]);
	}

	// delete the key associated with index i
	void deleteKey(int i)   {
		auto ind = index[i];
		swap(ind, N--);
		bubbleUp(ind);
		bubbleDown(ind);
		index[i] = -1;
	}
};

class Edge;
using spEdge = shared_ptr<Edge>;
class Vertex {
public:
	Vertex() { }
	list<spEdge> edges;
	bool visited;
	int value;
	Vertex(const Vertex&) = delete;
	Vertex& operator=(const Vertex&) = delete;
};
using spVertex = shared_ptr<Vertex>;
class Edge {
public:
	Edge() { head = nullptr; tail = nullptr; }
	spVertex head;
	spVertex tail;
	int weight;
	Edge(const Edge&) = delete;
	Edge& operator=(const Edge&) = delete;
};
bool operator<(const shared_ptr<Edge>&lhs, const shared_ptr<Edge>& rhs) {
	return lhs->weight < rhs->weight;
}
class AdjList {
public:
	vector<spVertex> vertices;
	vector<spEdge> edges;
	void addVertex(spVertex vertex) {
		vertices[vertex->value] = vertex;
	}
	void connect(int v1, int v2, int weight) {
		auto e = make_shared<Edge>();
		e->head = vertices[v1];
		e->tail = vertices[v2];
		e->weight = weight;
		e->head->edges.push_back(e);
		e->tail->edges.push_back(e);
		edges.push_back(e);
	}
	void resize(int v) {
		vertices.resize(v);
	}
	void resetVisited() {
		for (auto& v : vertices) {
			v->visited = false;
		}
	}
};

class UnionFind {
	// maintain a partition of a set of objects
	// group = find(x);
	// union(ci,cj);
	// - linked structure for each connected component
	// - each component has arbitrary leader
	// - each vertex points to the leader of its component directly
	// - name of a component = leader vertex
	// - so can just check if they have same leader 
	//   to see if a cycle is created
	// ie: find(u) = find(v)
	// need to rewire leader components each time we merge
	// keep leader of bigger group, rewire all the smaller group in
	// just keep size field for each group
public:
	UnionFind(int size) {
		leader.resize(size); 
		for (int i = 0; i < size; i++) {
			leader[i] = make_shared <UnionNode>(-1, i, 1);
		}
	}
	int find(int v) {
		return leader[v]->leader;
	}
	// prefer i = parent, (ie:starting point)
	int combine(int i, int j) {
		if (leader[i]->size <= leader[j]->size) {
			return combineThem(i, j);
		}
		else {
			return combineThem(j, i);
		}
	}
private:
	int combineThem(int from, int to) {
		// push all the former leader's childern to new leader
		for (auto child : leader[from]->children) {
			++(leader[to]->size);
			child->leader = to;
			leader[to]->children.push_back(child);
		}
		// now push in the the former leader
		++(leader[to]->size);
		leader[from]->leader = to;
		leader[to]->children.push_back(leader[from]);

		// we are our own leader (may be first time setting)
		leader[to]->leader = to;
		// old from has no children now, so empty it and its size.
		leader[from]->size = 0;
		leader[from]->children.clear();
		return to;
	}
	class UnionNode {
	public:
		UnionNode(int l, int v, int s) : leader(l), value(v), size(s) {}
		int leader;
		int value;
		int size;
		vector<shared_ptr<UnionNode>> children;
	};
	vector<shared_ptr<UnionNode>> leader;
};

class LazyUnionFind {
	// maintain a partition of a set of objects
	// group = find(x);
	// union(ci,cj);
	// - linked structure for each connected component
	// - each component indicates linkage via root of a tree
	// - when doing union, just assign top node to root
	//  - need to do find first to find each root
	//  - we know root by fact it points to self
	// - but the naive implementation will be potentially log(n) because union
	//   could potentially build a linear tree.
	//   - so we need to keep the depth of a tree and merge into larger one
	//   - maintain the field rank[node x] for each x in X
	// - ie: maximum number of hops from some leaf to node x
	//   initial rank = 0
	// - at this point, find = log n, combine = log n.. so not that great.

	// but Path Compression makes it better
	//   - as we traverse up, we want to rewire all the items up the chain to the same parent
	//   - shallow bushy tree instead of balanced binary

	// Rank maintenance is now a little different, it does not apply now to the 
	// actual cost of a find search, ie: rank of a leaf can be >0
	// However, analysis still stands more or less

	// hopcroft ullman says this data structure works on m log*(n) (so constant in m)
	// but "tarjan" bound says O(m * alpha(n))
	// inverse ackermann function
	//   - grows much more slowly than log*n
	//   
public:
	LazyUnionFind(int size) {
		leader.resize(size);
		rank.resize(size);
		for (int i = 0; i < size; i++) {
			leader[i] = make_shared <LazyUnionNode>(i, i, 1);
			rank[i] = 0;
		}
	}
	// traverse until hit a root
	int find(int v) {
		while (v != leader[v]->leader)
			return find(leader[v]->leader);
	}
	// prefer i = parent, (ie:starting point)
	int combine(int i, int j) {
		auto s1 = find(i);
		auto s2 = find(j);
		if (rank[s1] > rank[s2]) {
			leader[s2]->leader = s1;
		}
		else {
			if (rank[s1] == rank[s2]) {
				++rank[s2];
			}
			leader[s1]->leader = s2;
		}
	}
private:
	vector<int> rank;
	class LazyUnionNode {
	public:
		LazyUnionNode(int l, int v, int s) : leader(l), value(v), size(s) {}
		int leader;
		int value;
		int size;
		vector<shared_ptr<LazyUnionNode>> children;
	};
	vector<shared_ptr<LazyUnionNode>> leader;
};

class Graph {
public:
	Graph() {}
	Graph(int vertices) {
		adjList.resize(vertices);
		for (auto i = 0; i < vertices; ++i) {
			auto v = make_shared<Vertex>();
			v->value = i;
			v->visited = false;
			addVertex(v);
		}
	}
	void addVertex(spVertex vertex) {
		N++;
		adjList.addVertex(vertex);
	}
	void connect(int v1, int v2, int weight) {
		M++;
		adjList.connect(v1, v2, weight);
	}
	// disallow override by char, (just to see that happening)
	void connect(char v1, char v2, char weight) = delete;
	
	void resetVisited(){
		visitOrder.clear();
		adjList.resetVisited();
	}
	bool connected(int v1, int v2, int& weight, int& hops) {
		return false;
	}
	bool isSparse() {
		return N < pow(M, (1.5));
	}
	bool fullyConnected() {
		return true;
	}
	void connectedComponents(vector<Graph> &connected) {
		//TODO: need proper copy by value semantic here
		connected.push_back(*this);
		for (auto& g : connected) {
			g.resetVisited();
		}
	}
	
	void DFS(){
		resetVisited();
		// start from root, then take each first edge's tail, then backtrack as needed 
		// to layers in stack
		// - mazelike
		// - DAG topology
		// - connected components in a directed graph
		stack<spVertex> vStack;
		auto& start = adjList.vertices[0];
		start->visited = true;
		visitOrder.push_back(start);
		vStack.push(start);
		while (!vStack.empty()) {
			auto& v = vStack.top();
			vStack.pop();
			for (auto& e : v->edges) {
				auto& w = e->tail;
				if (!w->visited) {
					visitOrder.push_back(w);
					w->visited = true;
					vStack.push(w);
				}
			}
		}
	}
	void BFS(){
		resetVisited();
		// start from root, then take each neighbor one at a time until no more, 
		// then start with first set of visited neighbors waiting in queue
		// - shortest path 
		// - connected undirected components
		queue<spVertex> vQueue;
		auto& start = adjList.vertices[0];
		start->visited = true;
		visitOrder.push_back(start);
		vQueue.push(start);
		while (!vQueue.empty()) {
			auto& v = vQueue.front();
			vQueue.pop();
			for (auto& e : v->edges) {
				auto& w = e->tail;
				if (!w->visited) {
					visitOrder.push_back(w);
					w->visited = true;
					vQueue.push(w);
				}
			}
		}
	}
	unique_ptr<pair<Graph, Graph>> cut() {
		// current graph is not good for this with static vertex init.. 
		Graph A(N);
		Graph B(N);
		return make_unique<pair<Graph, Graph>>(move(A), move(B));
	}
	// shortestLength can work with simple BFS only if all weights= 1
	// so use Dijkstra instead
	// X = [s] - vertices processed so far
	// A[s].minCost = 0 (use visitOrder in the end)
	unique_ptr<pair<vector<int>, vector<int>>> shortestLength(int startIndex, int endIndex) {
		resetVisited();
		vector<int> key(N, INT_MAX);
		vector<int> parent(N, -1);
		MinIndexedPQ<int> heap(N);
		auto& startVertex = adjList.vertices[startIndex];
		key[0] = 0;
		startVertex->visited = true;
		visitOrder.push_back(startVertex);
		heap.insert(0, 0);
		for (auto v = 1; v < N; v++) {
			heap.insert(v, INT_MAX);
		}
		while (!heap.isEmpty() && parent[endIndex] == -1) {
			auto uCost = heap.minKey();
			auto uIndex = heap.deleteMin();
			auto& uVert = adjList.vertices[uIndex];
			for (auto& e : uVert->edges) {
				auto& vVert = e->tail;
				auto vIndex = vVert->value;
				if (heap.contains(vIndex) && uCost + e->weight < key[vIndex]) {
					parent[vIndex] = uIndex;
					vVert->visited = true;
					visitOrder.push_back(vVert);
					key[vIndex] = uCost + e->weight;
					heap.decreaseKey(vIndex, key[vIndex]);
				}
				if (vIndex == endIndex)
					break;
			}
		}
		return make_unique<pair<vector<int>, vector<int>>>(move(key), move(parent));
	}
	int numVertices() {
		return N;
	}
	int numEdges() {
		return M;
	}
	// undirected graphs
	unique_ptr<Graph> prims() {
		// growing one vertex at a time.
		// add one edge, span one vertex adjacent
		// try cheapest edge from visited vertices to unvisited
		resetVisited();
		Graph out(numVertices());
		// could do keys = edge cost, but better way..
		// that would require also checking on whether the edge crosses frontier
		
		// instead, storing vertices as in dijkstra solution
		// - elements in heap = vertices of V - X (not yet spanned)
		// - for v in V - X, key[v] = cheapest edge (u,v) with v in X
		MinIndexedPQ<shared_ptr<Edge>> heap(N-1);
		vector<shared_ptr<Vertex>> X;
		// X = [S] (one specific vertex)
		X.push_back(adjList.vertices[0]);
		adjList.vertices[0]->visited = true;
		// T = 0
		vector<shared_ptr<Edge>> T;
		// need to add the cheapest edge for each key[v]
		for (auto e : adjList.vertices[0]->edges) {
			auto vEdge = e->tail;
			auto vIndex = vEdge->value;
			if (e->head->visited && !vEdge->visited) {
				if (!heap.contains(vIndex)) {
					heap.insert(vIndex, e);
				}
				else if (heap.keyOf(vIndex)->weight > e->weight) {
					heap.decreaseKey(vIndex, e);
				}
			}
		}
		// while X != V (all original vertices)
		while (X.size() != adjList.vertices.size()) {
			//  - let e=(u,v) -> cheapest edge of G with u in X, v !in X
			// how do we get this edge given the vertex only? 
			// need a heap that stores edges not just the weight.
			auto uEdge = heap.minKey();
			auto vIndex = heap.deleteMin();
			//  - add e to T
			T.push_back(uEdge);
			//  - add v to X
			auto vVertex = uEdge->tail;
			vVertex->visited = true;
			X.push_back(vVertex);
			// maintain the invariant for next loop
			// for all nodes that cross now from "v" into any other vertex
			// remaining, we should compare if this is a smaller weight from this 'v'
			for (auto e : vVertex->edges) {
				auto wEdge = e->tail;
				auto wIndex = wEdge->value;
				if (e->head->visited && !wEdge->visited) {
					// could delete, recompute, and re-insert
					if (!heap.contains(wIndex)) {
						heap.insert(wIndex, e);
					}
					else if (heap.keyOf(wIndex)->weight > e->weight) {
						heap.decreaseKey(wIndex, e);
					}
				}
			}
		}
		// connect all the collected edges for our output graph
		for (auto e : T) {
			out.connect(e->head->value, e->tail->value, e->weight);
		}
		return make_unique<Graph>(move(out));
	}

	unique_ptr<Graph> primsPoly() {
		// growing one vertex at a time.
		// add one edge, span one vertex adjacent
		// try cheapest edge from visited vertices to unvisited
		resetVisited();
		Graph out(numVertices());
		vector<shared_ptr<Vertex>> X;
		// X = [S] (one specific vertex)
		X.push_back(adjList.vertices[0]);
		// T = 0
		vector<shared_ptr<Edge>> T;
		// while X != V (all original vertices)
		while (X.size() != adjList.vertices.size()) {
			//  - let e=(u,v) -> cheapest edge of G with u in X, v !in X
			auto minWeight = INT_MAX;
			shared_ptr<Edge> minEdge = nullptr;
			for (auto u : X) {
				for (auto e : u->edges) {
					auto v = e->tail;
					if (!v->visited) {
						minWeight = min(minWeight, e->weight);
						minEdge = e;
					}
				}
			}
			//  - add e to T
			T.push_back(minEdge);
			//  - add v to X
			minEdge->tail->visited = true;
			X.push_back(minEdge->tail);
		}
		// connect all the collected edges for our output graph
		for (auto e : T) {
			out.connect(e->head->value, e->tail->value, e->weight);
		}
		return make_unique<Graph>(move(out));
	}
	// uses union-find data structure
	unique_ptr<Graph> kruskal(int k = -1) {
		// grows tree by pieces and coalesces at end
		// just start grabbing the smallest edges one at a time
		// but don't violate our assumptions while doing it
		// so. could also sort these first and then do it. 
		resetVisited();
		Graph out(numVertices());
		UnionFind unionFind(numVertices());
		vector<shared_ptr<Edge>> edges = adjList.edges;
		sort(begin(edges), end(edges));
		// sort edges.
		// for i = 1 to M
		auto inserted = 0;
		for (auto e : edges) {
			if (k > 0 && inserted == (k - 1))
				break;
			inserted++;
			//   - if T + [i] has no cycles, (union find structure?)
			//    (could do DFS/BFS as well)
			auto uGroup = unionFind.find(e->head->value);
			auto vGroup = unionFind.find(e->tail->value);
			// group == same -> have not assigned to another group yet
			if (-1 == vGroup) {
				out.connect(e->head->value, e->tail->value, e->weight);
				unionFind.combine(e->head->value, e->tail->value);
			}
			// union find approach:
			//  - group = find(vertex) (already there?)
			//  - add new edge(u,v) = union(group(u),group(v))
		}
		// return T
		return make_unique<Graph>(move(out));
	}

	unique_ptr<Graph> kruskalDFS() {
		// grows tree by pieces and coalesces at end
		// just start grabbing the smallest edges one at a time
		// but don't violate our assumptions while doing it
		// so. could also sort these first and then do it. 
		resetVisited();
		Graph out(numVertices());
		vector<shared_ptr<Edge>> edges = adjList.edges;
		sort(begin(edges), end(edges));
		// sort edges.
		// for i = 1 to M
		for (auto e : edges) {
			//   - if T + [i] has no cycles, (union find structure?)
			//    (could do DFS/BFS as well)

			// union find approach:
			//  - group = find(vertex) (already there?)
			//  - add new edge(u,v) = union(group(u),group(v))
		}
		// return T
		return make_unique<Graph>(move(out));
	}
	unique_ptr<Graph> maxSpaceKClusters(int k) {
		// Max-Spacing k-clusters
		// - spacing = distance between the minimum p,q where p in set a, q in set b
		// - given distance measure and # of clusters k, create largest spacing
		// step 1: initialize where each vector is in own cluster
		// - assuming we are still > k clusters
		//   - find closest pair of separated points and fuse
		//   - this is kruskal's algorithm stopping at k-1 edges
		//    - points = vertices, distances = edge-costs, point pairs = edges
		
		// Single Link Clustering
	};


private:
	AdjList adjList;
	int N;
	int M;
	vector<spVertex> visitOrder;
};

class Cache {
	// small fast memory vs. big slow memory
	// - process sequence of page requests
	// - on a fault, what do we evict to bring in new data?
	// - the "furthest in the future" algorithm is an optimal greedy algorithm
	//    - if you know in advance that you'd evict A in 7 steps, and B in 70 steps, evict B now.
	//    - of course we don't know the future, so not implementable... 
	// - this is the most ideal guideline, and motivates Least-recently-used
	//    - so look in the past, and guess that what was requested recently will be requested again soon
	//    - so evict the oldest data
	// - and is a way to measure after the fact (or, well, in fact, to tune over time on strategy)
	//    - data locality of reference works on LRU
	//    - but if we know we are doing worst, can switch, for instance, to MRU (evict most recent)
	//    - (often MRU can be used, for instance, if it is VERY expensive to draw from far past, 
	//       or if the particular use case is legacy history.. etc. 
	//        -- for instance, 2 weeks ago a VERY popular article came out, so use MRU for now...)
	// - proof = very difficult exchange argument

};

// C++ 11 cheat until 14
/*
template <class C>
auto cbegin(const C& container)->decltype(std::begin(container))
{
	// if container is const, this actually returns a const iterator.
	return std::begin(container);
}

template <class C>
auto cend(const C& container)->decltype(std::end(container))
{
	return std::end(container);
}
*/

class Scheduler {
public:
	// one shared resource: ie - processor
	// many jobs to do (N)
	// - order of sequence of jobs
	// each job as:
	//  - weight (w_[j]) "priority"
	//  - length (l_[j]) "time"
	// Completion time C_[j] = sum(l_[0]...l_[j-1])+l[j]
	class Job {
	public:
		double weight;
		double length;
		// We are storing these directly in vector, so using default copy constructors
		//Job(const Job&) = delete;
		//Job& operator=(const Job&) = delete;
	};
	int N;
	// minimizing the weighted sum of completion times
	// - min sum(j=1:n,w_[j]*C_[j]);
	// edge cases:
	//  - if all jobs have same length, schedule higher weight earlier
	//  - if all have the same weight, schedule shorter earlier
	//  - if w_[i] > w_[j] && l_[i] > l_[j]
	//     - assign "scores" to jobs that are:
	//         - 1)increasing in weight
	//         - 2)decreasing in length
	//         - *)ties can be broken arbitrarily (ie: random flip or first in list.. etc.)
	//     - ie: weight - length or weight / length
	//         - so which of these works best?
	//         - example: {{3,5},{1,2}}
	//         - 1) score(-2,-1) = takes J2 first = minsum:23
	//         - 2) score(3/5,1/2) = takes J1 first = minsum:22 (better minsum)
	// - using #2 = n - log n runtime and is correct
	// - can prove correctness of ties with similar logic to "inverted pairs" problem
	// - (n choose 2 inversions) = n*(n-1)/2 (ie - we are using a bubble sort logic in our proof)

	unique_ptr<vector<int>> schedule(vector<Job> &jobs) {
		vector<int> completionTimes;
		// - use lambda comparison
		// - implement schedule sort solution
		// note: I had tried cbegin/cend, but.. of course... this is not constant!
		sort(begin(jobs), end(jobs),
			[](const Job &a, const Job &b)
			{return a.weight / a.length < b.weight / b.length; });
		completionTimes.resize(jobs.size());
		completionTimes[0] = jobs[0].length;
		for (auto i = 1; i < jobs.size(); i++) {
			completionTimes[i] = completionTimes[i - 1] + jobs[i].length;
		}
		return make_unique<vector<int>>(move(completionTimes));
	}
};

class Huffman {
	// builds a huffman encoding tree on inputs
	// one approach: top-down, divide and conquer, partition into two sets, with ~50% of freq each
	// but not optimal;
	// huffman: bottom-up, merge
	//   - nodes, leaves, each labelled from alphabet/universe-set
	// ie: in order of smallest items each time = the greedy aspect
	// the merged score is used to compute the next time through
	// since we are looking for "two smallest" = use heap!
	// but even faster: sort inputs, then just do o(n) work (but how to sort the merge??)
	// -- mergesort? (heap sort or tries -- wouldn't need to use heap data structure, can use queue
	//                or two queues)
};


int _tmain(int argc, _TCHAR* argv[])
{
	Graph g(32);
	srand(chrono::system_clock::now().time_since_epoch().count());
	// connect 2 level 1 children to parent (p = 0, c = 1/2)
	// connect 4 level 2 children to parents (p = 1/2, c = 3/4/5/6)
	// (p=3/4/5/6, c = 7/8/9/10/11/12/13/14/15..)
	// create a balanced / breadth-first-ordered / directed graph
	for (auto parent = 0; parent < 16; ++parent) {
		for (auto child = (parent*2)+1; child <= (parent*2)+2 && child < 32; ++child) {
			// weight = 5 for the basic tree, so path-len = 5, weight = 25.
			g.connect(parent, child, 5);
		}
	}
	// now randomly connect with about 30 more edges with random weights <= 3
	for (auto edges = 0; edges < 30; ++edges) {
		int start = rand() % 31;
		int end = rand() % 31;
		while (end == start) {
			end = rand() % 31;
		}
		g.connect(start, end, rand() % 3);
	}
	g.DFS();
	g.BFS();
	auto paths = g.shortestLength(0, 31);
	Scheduler s;
	vector<Scheduler::Job> jobs{ { 3, 5 }, { 1, 2 }, { 4, 8 }, { 10, 1 }, { 18, 1 }, { 1, 18 }, { 1, 20 }, { 1, 1 }, { 20, 20 } };
	auto completionTimes = s.schedule(jobs);
	vector<Scheduler::Job> jobs2{ { 20, 20 }, { 1, 2 }, { 4, 8 }, { 18, 1 }, { 1, 18 }, { 1, 20 }, { 1, 1 }, { 3, 5 }, { 10, 1 } };
	auto completionTimes2 = s.schedule(jobs2);
	vector<Scheduler::Job> jobs3{ { 3, 5 }, { 3, 10 }, { 3, 1 }, { 3, 20 }, { 3, 1 }, { 3, 15}, { 3, 17 }, { 3, 20 }, { 1, 20 } };
	auto completionTimes3 = s.schedule(jobs3);
	Graph primMe(4);
	primMe.connect(0, 1, 1);
	primMe.connect(1, 2, 2);
	primMe.connect(2, 3, 3);
	primMe.connect(0, 3, 1);
	auto primmed = primMe.prims();
	Graph kruskalMe(4);
	kruskalMe.connect(0, 1, 1);
	kruskalMe.connect(1, 2, 2);
	kruskalMe.connect(2, 3, 3);
	kruskalMe.connect(0, 3, 1);
	auto kruskalled = kruskalMe.kruskal();
	return 0;
}


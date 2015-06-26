/**
 * Created by Cody on 6/24/2015.
 */
public class WeightedQuickUnion implements UnionFind {
    private int[] id;
    private int[] size;

    WeightedQuickUnion(int N) {
        id = new int[N];
        size = new int[N];
        for (int i=0; i < N; i++) {
            id[i] = i;
            size[i] = 1;
        }
    }
    public void dumpIds() {
        for (int i = 0; i < id.length; i++) {
            StdOut.print(id[i] + " ");
        }
        StdOut.println("");
    }

    public boolean find(int p, int q) { return root(p) == root(q); }
    // take steps to avoid tall trees during union
    // keep track of size of each tree
    // balance by linking root of smaller tree to root of larger tree (can also use height or 'rank')
    public void union(int p, int q) {
        // set the id of p's root to the id of q's root
        int i = root(p);
        int j = root(q);
        if (i == j) return;
        if (size[i] >= size[j]) {
            id[j] = id[i];
            size[i] += size[j];
        }
        else {
            id[i] = id[j];
            size[j] += size[i];
        }
    }
    private int root(int i) {
        while (id[i] != i) i = id[i];
        return i;
    }
}

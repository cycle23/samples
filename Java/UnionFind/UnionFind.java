/**
 * Created by Cody on 6/24/2015.
 */
public interface UnionFind {
    public boolean find(int p, int q);
    public void union(int p, int q);
    public void dumpIds();

}

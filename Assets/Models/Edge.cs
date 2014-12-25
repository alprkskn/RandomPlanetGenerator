using System.Collections.Generic;

public class Edge
{
    public int[] n;
    public List<int> f = new List<int>();
    public List<int> subdivided_n = new List<int>();
    public List<int> subdivided_e = new List<int>();

    public Edge()
    {
        n = new int[2];
    }

    public Edge(int start, int end)
    {
        this.n = new int[2] { start, end };
    }

    public Edge(int[] e)
    {
        this.n = e;
    }


}
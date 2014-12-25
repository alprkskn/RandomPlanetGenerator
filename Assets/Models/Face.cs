using System.Collections.Generic;
using UnityEngine;

public class Face
{
    public List<int> n = new List<int>();
    public List<int> e = new List<int>();
    public Vector3 centroid;

    public Face()
    {

    }

    public Face(int[] n, int[] e)
    {
        this.n = new List<int>(n);
        this.e = new List<int>(e);
    }
}
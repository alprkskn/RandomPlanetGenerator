using System.Collections.Generic;
using UnityEngine;

public class Node
{
    public Vector3 p;
    public List<int> e = new List<int>();
    public List<int> f = new List<int>();

    public Node()
    {

    }

    public Node(Vector3 v)
    {
        this.p = v;
    }
}
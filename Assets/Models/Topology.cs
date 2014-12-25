using System.Collections.Generic;
using UnityEngine;

public class Topology
{
    public List<Corner> corners;
    public List<Border> borders;
    public List<Tile> tiles;

    public Topology()
    {
        corners = new List<Corner>();
        borders = new List<Border>();
        tiles = new List<Tile>();
    }

    public Topology(List<Corner> corners, List<Border> borders, List<Tile> tiles)
    {
        this.corners = corners;
        this.borders = borders;
        this.tiles = tiles;
    }

    public void DebugDraw(float scale, Color color, float duration)
    {
        foreach(var b in borders)
        {
            if(b.corners.Length != 2)
            {
                Debug.Log("LOL - " + b.corners.Length);
            }
            else
            {
                Debug.DrawLine(b.corners[0].position * scale, b.corners[1].position * scale, color, duration);
            }
        }
    }
}
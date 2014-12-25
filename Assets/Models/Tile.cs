using UnityEngine;

public class Tile
{
    public int id;
    public Vector3 position;
    public Vector3 averagePosition;
    public Vector3 normal;
    public float area;
    public float radius;
    public Corner[] corners;
    public Border[] borders;
    public Tile[] tiles;

    public Tile(int id, Vector3 position, int cornerCount, int borderCount, int tileCount)
    {
        this.id = id;
        this.position = position;
        this.corners = new Corner[cornerCount];
        this.borders = new Border[borderCount];
        this.tiles = new Tile[tileCount];
    }
}
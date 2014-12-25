using UnityEngine;

public class Corner
{
    public int id;
    public Vector3 position;
    public Corner[] corners;
    public Border[] borders;
    public Tile[] tiles;
    public float area;

    public Corner(int id, Vector3 position, int cornerCount, int borderCount, int tileCount)
    {
        this.id = id;
        this.position = position;
        this.corners = new Corner[cornerCount];
        this.borders = new Border[borderCount];
        this.tiles =   new Tile[tileCount];
    }

    public Vector3 vectorTo(Corner corner)
    {
        return corner.position - this.position;
    }
}
using UnityEngine;

public class Border
{
    int id;
    public Corner[] corners;
    public Border[] borders;
    public Tile[] tiles;
    public Vector3 midpoint;

    public Border(int id, int cornerCount, int borderCount, int tileCount)
    {
        this.id = id;
        this.corners = new Corner[cornerCount];
        this.borders = new Border[borderCount];
        this.tiles = new Tile[tileCount];
    }

    public Corner oppositeCorner(Corner corner)
    {
        return (this.corners[0] == corner) ? this.corners[1] : this.corners[0];
    }

    public Tile oppositeTile(Tile tile)
    {
        return (this.tiles[0] == tile) ? this.tiles[1] : this.tiles[0];
    }

    public float length()
    {
        return Vector3.Distance(this.corners[0].position, this.corners[1].position);
    }
}
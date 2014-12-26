using System;
using System.Collections.Generic;
using UnityEngine;

public class Polyhedra
{
    public List<Node> nodes;
    public List<Edge> edges;
    public List<Face> faces;
    public Topology topology = null;

    public Polyhedra(List<Node> n, List<Edge> e, List<Face> f)
    {
        this.nodes = n;
        this.edges = e;
        this.faces = f;
    }

    public void DebugDraw(float scale, Color color, float duration)
    {
        foreach(var e in edges)
        {
            Debug.DrawLine(nodes[e.n[0]].p * scale, nodes[e.n[1]].p * scale, color, duration);
        }
    }
    public void generatePlanetTopology()
    {
        List<Corner> corners = new List<Corner>(this.faces.Count);
        List<Border> borders = new List<Border>(this.edges.Count);
        List<Tile> tiles = new List<Tile>(this.nodes.Count);
        int i;
        for(i = 0; i < this.faces.Count; ++i)
        {
            var face = this.faces[i];
            corners.Add(new Corner(i, face.centroid * 1000, face.e.Count, face.e.Count, face.n.Count));
        }

        for(i = 0; i < this.edges.Count; ++i)
        {
            var edge = this.edges[i];
            borders.Add(new Border(i, 2, 4, 2)); //edge.f.Count, this.faces[edge.f[0]].e.Count + this.faces[edge.f[1]].e.Count - 2, edge.n.Count
        }
        
        for(i = 0; i < this.nodes.Count; ++i)
        {
            var node = this.nodes[i];
            tiles.Add(new Tile(i, node.p * 1000, node.f.Count, node.e.Count, node.e.Count));
        }

        for(i = 0; i < corners.Count; ++i)
        {
            var corner = corners[i];
            var face = this.faces[i];
            for(var j = 0; j < face.e.Count; ++j)
            {
                corner.borders[j] = (borders[face.e[j]]);
            }
            for(var j = 0; j < face.n.Count; ++j)
            {
                corner.tiles[j] = (tiles[face.n[j]]);
            }
        }

        for(i = 0; i < borders.Count; ++i)
        {
            var border = borders[i];
            var edge = this.edges[i];
            var averageCorner = new Vector3(0, 0, 0);
            var n = 0;
            for (var j = 0; j < edge.f.Count; ++j)
            {
                var corner = corners[edge.f[j]];
                averageCorner += corner.position;
                border.corners[j] = corner;
                for (var k = 0; k < corner.borders.Length; ++k)
                {
                    if (corner.borders[k] != border) border.borders[n++] = corner.borders[k];
                }
            }
            border.midpoint = averageCorner / border.corners.Length;
            for (var j = 0; j < edge.n.Length; ++j)
            {
                border.tiles[j] = tiles[edge.n[j]];
            }
        }
        
        for (i = 0; i < corners.Count; ++i)
        {
            var corner = corners[i];
            for (var j = 0; j < corner.borders.Length; ++j)
            {
                corner.corners[j] = corner.borders[j].oppositeCorner(corner);
            }
        }

        for (i = 0; i < tiles.Count; ++i)
        {
            var tile = tiles[i];
            var node = this.nodes[i];
            for (var j = 0; j < node.f.Count; ++j)
            {
                tile.corners[j] = corners[node.f[j]];
            }
            for (var j = 0; j < node.e.Count; ++j)
            {
                var border = borders[node.e[j]];
                if (border.tiles[0] == tile)
                {
                    for (var k = 0; k < tile.corners.Length; ++k)
                    {
                        var corner0 = tile.corners[k];
                        var corner1 = tile.corners[(k + 1) % tile.corners.Length];
                        if (border.corners[1] == corner0 && border.corners[0] == corner1)
                        {
                            border.corners[0] = corner0;
                            border.corners[1] = corner1;
                        }
                        else if (border.corners[0] != corner0 || border.corners[1] != corner1)
                        {
                            continue;
                        }
                        tile.borders[k] = border;
                        tile.tiles[k] = border.oppositeTile(tile);
                        break;
                    }
                }
                else
                {
                    for (var k = 0; k < tile.corners.Length; ++k)
                    {
                        var corner0 = tile.corners[k];
                        var corner1 = tile.corners[(k + 1) % tile.corners.Length];
                        if (border.corners[0] == corner0 && border.corners[1] == corner1)
                        {
                            border.corners[1] = corner0;
                            border.corners[0] = corner1;
                        }
                        else if (border.corners[1] != corner0 || border.corners[0] != corner1)
                        {
                            continue;
                        }
                        tile.borders[k] = border;
                        tile.tiles[k] = border.oppositeTile(tile);
                        break;
                    }
                }
            }

            tile.averagePosition = new Vector3(0, 0, 0);
            for (var j = 0; j < tile.corners.Length; ++j)
            {
                tile.averagePosition += tile.corners[j].position;
            }
            tile.averagePosition /= tile.corners.Length;
			
            var maxDistanceToCorner = 0f;
            for (var j = 0; j < tile.corners.Length; ++j)
            {
                maxDistanceToCorner = Math.Max(maxDistanceToCorner, Vector3.Distance(tile.corners[j].position, tile.averagePosition));
            }
			
            var area = 0f;
            for (var j = 0; j < tile.borders.Length; ++j)
            {
                area += PlanetGenerate.CalculateTriangleArea(tile.position, tile.borders[j].corners[0].position, tile.borders[j].corners[1].position);
            }
            tile.area = area;
			
            tile.normal = tile.position.normalized;
            tile.radius = maxDistanceToCorner;
        }

        for (i = 0; i < corners.Count; ++i)
        {
            var corner = corners[i];
            corner.area = 0;
            for (var j = 0; j < corner.tiles.Length; ++j)
            {
                corner.area += corner.tiles[j].area / corner.tiles[j].corners.Length;
            }
        }

        this.topology = new Topology(corners, borders, tiles);
    }
    

}
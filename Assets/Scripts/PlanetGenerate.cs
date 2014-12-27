using System.Linq;
using UnityEngine;
using System.Collections;
using System;
using System.Collections.Generic;
using Random = System.Random;

[RequireComponent(typeof(MeshFilter), typeof(MeshRenderer))]
public class PlanetGenerate : MonoBehaviour 
{
    Random globalRandom;
    public int SubdivisionLevel = 2;
    public float DistortionLevel = 1f;
    int plateCount = 36;
    float oceanicRate = 0.7f;
    float heatLevel = 1.0f;
    float moistureLevel = 1.0f;

    Polyhedra planet;
	// Use this for initialization
	void Start ()
    {

        //Polyhedra p = generateIcosahedron();
        //p.DebugDraw(5f, Color.red, 10000f);

        ////Polyhedra p2 = generateSubdividedIcosahedron(2);
        ////p2.DebugDraw(10f, Color.green, 10000f);
        globalRandom = new Random(4);
        //Polyhedra p10 = Utils.generateSubdividedIcosahedron(20);
        //p10.DebugDraw(1000f, Color.yellow, 10000f);


        planet = generatePlanetMesh(SubdivisionLevel, DistortionLevel);
        //planet.DebugDraw(20f, Color.green, 10000f);
        planet.generatePlanetTopology();
        planet.topology.generatePlanetPartition();
        //planet.DebugDraw(40f, Color.red, 10000f);
        //planet.topology.DebugDraw(1f, Color.red, 10000f);
        Debug.Log(planet.topology.tiles.Count);
	    var verts = new List<Vector3>();
        var norms = new List<Vector3>();
        var tris = new List<int>();

	    var tileStart = 0;
	    foreach (var tile in planet.topology.tiles)
	    {
            if(tileStart >= 64000)
            {
                GameObject go = new GameObject("child");
                go.transform.parent = this.gameObject.transform;

                var mesh = new Mesh { name = "World Mesh" };
                go.AddComponent<MeshFilter>().mesh = mesh;
                var mr = go.AddComponent<MeshRenderer>();
                mr.material = new Material(Shader.Find("Diffuse"));
                mesh.vertices = verts.ToArray();
                //mesh.uv = UV;
                mesh.triangles = tris.ToArray();
                mesh.normals = norms.ToArray();

                tileStart = 0;
                verts.Clear();
                norms.Clear();
                tris.Clear();
            }
            verts.Add(tile.averagePosition);
            norms.Add(tile.averagePosition.normalized);
	        tileStart = verts.Count - 1;
	        foreach (var corner in tile.corners)
	        {
	            verts.Add(corner.position);
                norms.Add(tile.averagePosition.normalized);

                ////brnDEBUG
                //var cubeMark = GameObject.CreatePrimitive(PrimitiveType.Sphere);
                //cubeMark.transform.position = corner.position;
                //cubeMark.renderer.material.color = new Color(255, 0, 0);
	        }

            for (int i = 1; i <= tile.corners.Length; i++)
            {
                var mid = (i + 1) % (tile.corners.Length + 1) != 0
                    ? i + 1
                    : 1;

                tris.Add(tileStart);
                tris.Add(tileStart + mid);
                tris.Add(tileStart + i);

            }
	    }

        GameObject ch = new GameObject("child");
        ch.transform.parent = this.gameObject.transform;
        var mesh2 = new Mesh { name = "World Mesh" };
        ch.AddComponent<MeshFilter>().mesh = mesh2;
        var mr2 = ch.AddComponent<MeshRenderer>();
        mr2.material = new Material(Shader.Find("Diffuse"));
        mesh2.vertices = verts.ToArray();
        //mesh.uv = UV;
        mesh2.triangles = tris.ToArray();
        mesh2.normals = norms.ToArray();

        //GameObject go = new GameObject("icosahedron");
        //MeshFilter mf = go.AddComponent<MeshFilter>();
        //MeshRenderer mr = go.AddComponent<MeshRenderer>();


        //mf.mesh = generateIcosahedron();
        //mr.material = new Material(Shader.Find("Diffuse"));

        //go.transform.position = Vector3.zero;
	}
	
	// Update is called once per frame
	void Update () {
	    if(Input.GetMouseButtonDown(0))
        {
            int index = globalRandom.Next(0, planet.topology.partition.partitions.Count);
            Color c = new Color((float)globalRandom.NextDouble(), (float)globalRandom.NextDouble(), (float)globalRandom.NextDouble());
            foreach(var t in planet.topology.partition.tiles)
            {
                t.DebugDraw(c, 500f, 50f);
            }
            foreach(var p in planet.topology.partition.partitions)
            {
                c = new Color((float)globalRandom.NextDouble(), (float)globalRandom.NextDouble(), (float)globalRandom.NextDouble());
                foreach(var t in p.tiles)
                {
                    t.DebugDraw(c, 500f, 50f);
                }
            }
        }
	}

    void generatePlanet(int icosahedronSubdivision, float topologyDistortionRate, int plateCount) 
    {
        var icosahedron = Utils.generateIcosahedron();
        for(var i = 0; i < icosahedron.faces.Count; ++i)
        {
            var face = icosahedron.faces[i];
        }

    }

    Polyhedra generatePlanetMesh(int icosahedronSubdivision, float topologyDistortionRate)
    {
        var mesh = Utils.generateSubdividedIcosahedron(icosahedronSubdivision);
        var totalDistortion = Math.Ceiling(mesh.edges.Count * topologyDistortionRate);
		var remainingIterations = 6;
        int i;
        while(remainingIterations > 0)
        {
			var iterationDistortion = (int)Math.Floor(totalDistortion / remainingIterations);
			totalDistortion -= iterationDistortion;
			distortMesh(mesh, iterationDistortion, globalRandom);
			relaxMesh(mesh, 0.5f);
			--remainingIterations;
        }
        
        //var initialIntervalIteration = action.intervalIteration;
	
		var averageNodeRadius = Math.Sqrt(4 * Math.PI / mesh.nodes.Count);
		var minShiftDelta = averageNodeRadius / 50000 * mesh.nodes.Count;
		var maxShiftDelta = averageNodeRadius / 50 * mesh.nodes.Count;

		float priorShift, shiftDelta;
		var currentShift = relaxMesh(mesh, 0.5f);

        do
        {
            priorShift = currentShift;
            currentShift = relaxMesh(mesh, 0.5f);
            shiftDelta = Math.Abs(currentShift - priorShift);
        } while (shiftDelta >= minShiftDelta /* && action.intervalIteration - initialIntervalIteration < 300*/);

        for(i = 0; i < mesh.faces.Count; ++i)
        {
            var face = mesh.faces[i];
            var p0 = mesh.nodes[face.n[0]].p;
            var p1 = mesh.nodes[face.n[1]].p;
            var p2 = mesh.nodes[face.n[2]].p;
            face.centroid = calculateFaceCentroid(p0, p1, p2).normalized;
        }

        for(i = 0; i < mesh.nodes.Count; ++i)
        {
            var node = mesh.nodes[i];
            var faceIndex = node.f[0];
            for(var j = 1; j < node.f.Count - 1; ++j)
            {
                faceIndex = findNextFaceIndex(mesh, i, faceIndex);
                var k = node.f.IndexOf(faceIndex);
                node.f[k] = node.f[j];
                node.f[j] = faceIndex;
            }
        }

        return mesh;
    }

    bool distortMesh(Polyhedra mesh, int degree, Random random)
    {
        double totalSurfaceArea = 4 * Math.PI;
        double idealFaceArea = totalSurfaceArea / mesh.faces.Count;
        double idealEdgeLength = Math.Sqrt(idealFaceArea * 4 / Math.Sqrt(3));
        double idealFaceHeight = idealEdgeLength * Math.Sqrt(3) / 2;

        var rotationPredicate = new Func<Node, Node, Node, Node, bool>((oldNode0, oldNode1, newNode0, newNode1) =>
        {
            if (newNode0.f.Count >= 7 ||
                newNode1.f.Count >= 7 ||
                oldNode0.f.Count <= 5 ||
                oldNode1.f.Count <= 5) return false;

            var oldEdgeLength = Vector3.Distance(oldNode0.p, oldNode1.p);
            var newEdgeLength = Vector3.Distance(newNode0.p, newNode1.p);
            var ratio = oldEdgeLength / newEdgeLength;
            if (ratio >= 2 || ratio <= 0.5) return false;
            var v0 = (oldNode1.p - oldNode0.p) / oldEdgeLength;
            var v1 = (newNode0.p - oldNode0.p).normalized;
            var v2 = (newNode1.p - oldNode0.p).normalized;
            if (Vector3.Dot(v0, v1) < 0.2 || Vector3.Dot(v0, v2) < 0.2) return false;
            v0 *= -1;
            var v3 = (newNode0.p - oldNode1.p).normalized;
            var v4 = (newNode1.p - oldNode1.p).normalized;
            if (Vector3.Dot(v0, v3) < 0.2 || Vector3.Dot(v0, v4) < 0.2) return false;

            return true;
        });

        var i = 0;

        while(i < degree)
        {
            var consecutiveFailedAttempts = 0;
            var edgeIndex = random.Next(0, mesh.edges.Count);
            while(!conditionalRotateEdge(mesh, edgeIndex, rotationPredicate))
            {
                if (++consecutiveFailedAttempts >= mesh.edges.Count) return false;
                edgeIndex = (edgeIndex + 1) % mesh.edges.Count;
            }
            ++i;
        }
        return true;
    }

    float relaxMesh(Polyhedra mesh, float multiplier)
    {
        var totalSurfaceArea = 4 * Math.PI;
        var idealFaceArea = totalSurfaceArea / mesh.faces.Count;
        var idealEdgeLength = Math.Sqrt(idealFaceArea * 4 / Math.Sqrt(3));
        var idealDistanceToCentroid = idealEdgeLength * Math.Sqrt(3) / 3 * 0.9;

        var pointShifts = new List<Vector3>(mesh.nodes.Count);
        int i = 0;
        for (i = 0; i < mesh.nodes.Count; i++)
        {
            pointShifts.Add(new Vector3(0, 0, 0));
        }

        i = 0;
        while(i < mesh.faces.Count)
        {
            var face = mesh.faces[i];
            var n0 = mesh.nodes[face.n[0]];
            var n1 = mesh.nodes[face.n[1]];
            var n2 = mesh.nodes[face.n[2]];
            var p0 = n0.p;
            var p1 = n1.p;
            var p2 = n2.p;
            var e0 = Vector3.Distance(p1, p0) / idealEdgeLength; // p1.distanceTo(p0) / idealEdgeLength;
            var e1 = Vector3.Distance(p2, p1) / idealEdgeLength; // p2.distanceTo(p1) / idealEdgeLength;
            var e2 = Vector3.Distance(p0, p2) / idealEdgeLength; // p0.distanceTo(p2) / idealEdgeLength;
            var centroid = calculateFaceCentroid(p0, p1, p2).normalized;
            var v0 = centroid - p0; // centroid.clone().sub(p0);
            var v1 = centroid - p1;
            var v2 = centroid - p2;
            var length0 = v0.magnitude;
            var length1 = v1.magnitude;
            var length2 = v2.magnitude;
            v0 *= (float)(multiplier * (length0 - idealDistanceToCentroid) / length0);
            v1 *= (float)(multiplier * (length1 - idealDistanceToCentroid) / length1);
            v2 *= (float)(multiplier * (length2 - idealDistanceToCentroid) / length2);
            pointShifts[face.n[0]] += (v0);
            pointShifts[face.n[1]] += (v1);
            pointShifts[face.n[2]] += (v2);

            ++i;            
        }

        var origin = new Vector3(0, 0, 0);
        for (i = 0; i < mesh.nodes.Count; ++i)
        {
            pointShifts[i] = (mesh.nodes[i].p + (Utils.ProjectPointOnPlane(mesh.nodes[i].p, origin, pointShifts[i]))).normalized;
        }

        var rotationSupressions = new List<float>();
        for(i = 0; i < mesh.nodes.Count; ++i)
        {
            rotationSupressions.Add(0);
        }

        i = 0;
        while(i<mesh.edges.Count)
        {
            var edge = mesh.edges[i];
            var oldPoint0 = mesh.nodes[edge.n[0]].p;
            var oldPoint1 = mesh.nodes[edge.n[1]].p;
            var newPoint0 = pointShifts[edge.n[0]];
            var newPoint1 = pointShifts[edge.n[1]];
            var oldVector = (oldPoint1 - oldPoint0).normalized;
            var newVector = (newPoint1 - newPoint0).normalized;
            var suppression = (1 - Vector3.Dot(oldVector, newVector)) * 0.5;
            rotationSupressions[edge.n[0]] = Math.Max(rotationSupressions[edge.n[0]], (float)suppression);
            rotationSupressions[edge.n[1]] = Math.Max(rotationSupressions[edge.n[1]], (float)suppression);

            ++i;
        }

        float totalShift = 0;
        for(i = 0; i < mesh.nodes.Count; ++i)
        {
            var node = mesh.nodes[i];
            var point = node.p;
            var delta = point;
            node.p = Vector3.Lerp(point, pointShifts[i], 1f - (float)Math.Sqrt(rotationSupressions[i])).normalized;
            delta -= node.p;
            totalShift += delta.magnitude;
        }

        return totalShift;
    }

    Vector3 calculateFaceCentroid(Vector3 pa, Vector3 pb, Vector3 pc)
    {
        var vabHalf = (pb - pa) / 2;
        var pabHalf = pa + vabHalf;
        var centroid = ((pc - pabHalf) / 3) + pabHalf;
        return centroid;
    }

    bool conditionalRotateEdge(Polyhedra mesh, int edgeIndex, Func<Node, Node, Node, Node, bool> predicate)
    {
        var edge = mesh.edges[edgeIndex];
        var face0 = mesh.faces[edge.f[0]];
        var face1 = mesh.faces[edge.f[1]];
        var farNodeFaceIndex0 = getFaceOppositeNodeIndex(face0, edge);
        var farNodeFaceIndex1 = getFaceOppositeNodeIndex(face1, edge);
        var newNodeIndex0 = face0.n[farNodeFaceIndex0];
        var oldNodeIndex0 = face0.n[(farNodeFaceIndex0 + 1) % 3];
        var newNodeIndex1 = face1.n[farNodeFaceIndex1];
        var oldNodeIndex1 = face1.n[(farNodeFaceIndex1 + 1) % 3];
        var oldNode0 = mesh.nodes[oldNodeIndex0];
        var oldNode1 = mesh.nodes[oldNodeIndex1];
        var newNode0 = mesh.nodes[newNodeIndex0];
        var newNode1 = mesh.nodes[newNodeIndex1];
        var newEdgeIndex0 = face1.e[(farNodeFaceIndex1 + 2) % 3];
        var newEdgeIndex1 = face0.e[(farNodeFaceIndex0 + 2) % 3];
        var newEdge0 = mesh.edges[newEdgeIndex0];
        var newEdge1 = mesh.edges[newEdgeIndex1];

        if (!predicate(oldNode0, oldNode1, newNode0, newNode1)) return false;

        oldNode0.e.RemoveAt(oldNode0.e.IndexOf(edgeIndex));
        oldNode1.e.RemoveAt(oldNode1.e.IndexOf(edgeIndex));
        newNode0.e.Add(edgeIndex);
        newNode1.e.Add(edgeIndex);

        edge.n[0] = newNodeIndex0;
        edge.n[1] = newNodeIndex1;

        newEdge0.f.RemoveAt(newEdge0.f.IndexOf(edge.f[1]));
        newEdge1.f.RemoveAt(newEdge1.f.IndexOf(edge.f[0]));
        newEdge0.f.Add(edge.f[0]);
        newEdge1.f.Add(edge.f[1]);

        oldNode0.f.RemoveAt(oldNode0.f.IndexOf(edge.f[1]));
        oldNode1.f.RemoveAt(oldNode1.f.IndexOf(edge.f[0]));
        newNode0.f.Add(edge.f[1]);
        newNode1.f.Add(edge.f[0]);

        face0.n[(farNodeFaceIndex0 + 2) % 3] = newNodeIndex1;
        face1.n[(farNodeFaceIndex1 + 2) % 3] = newNodeIndex0;

        face0.e[(farNodeFaceIndex0 + 1) % 3] = newEdgeIndex0;
        face1.e[(farNodeFaceIndex1 + 1) % 3] = newEdgeIndex1;
        face0.e[(farNodeFaceIndex0 + 2) % 3] = edgeIndex;
        face1.e[(farNodeFaceIndex1 + 2) % 3] = edgeIndex;

        return true;
    }

    int getFaceOppositeNodeIndex(Face face, Edge edge)
    {
        if (face.n[0] != edge.n[0] && face.n[0] != edge.n[1]) return 0;
        if (face.n[1] != edge.n[0] && face.n[1] != edge.n[1]) return 1;
        if (face.n[2] != edge.n[0] && face.n[2] != edge.n[1]) return 2;
        else return -1;
    }

    int getEdgeOppositeFaceIndex(Edge edge, int faceIndex)
    {
        if (edge.f[0] == faceIndex) return edge.f[1];
        if (edge.f[1] == faceIndex) return edge.f[0];
        else return -1;
    }

    int findNextFaceIndex(Polyhedra mesh, int nodeIndex, int faceIndex)
    {
        var node = mesh.nodes[nodeIndex];
        var face = mesh.faces[faceIndex];
        var nodeFaceIndex = face.n.IndexOf(nodeIndex);
        var edge = mesh.edges[face.e[(nodeFaceIndex + 2) % 3]];
        return getEdgeOppositeFaceIndex(edge, faceIndex);
    }
}

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
public class Face
{
    public List<int> n = new List<int>();
    public List<int> e = new List<int>();
    public List<Tile> children = new List<Tile>();
    public Vector3 centroid;
    
    public Vector3 sphereCenter;
    public float sphereRadius;

    public Face()
    {

    }

    public Face(int[] n, int[] e)
    {
        this.n = new List<int>(n);
        this.e = new List<int>(e);
    }
}
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

    // TODO: Replace hardcoded 1000s with a scale variable.
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
				area += Utils.calculateTriangleArea(tile.position, tile.borders[j].corners[0].position, tile.borders[j].corners[1].position);
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

    public bool intersectRay(Ray ray)
    {
        if (!Utils.intersectRayWithSphere(ray, this.averagePosition, this.radius)) return false;

        if (Utils.SignedDistancePlanePoint(this.normal, this.averagePosition, ray.origin) <= 0) return false;

        var denominator = Vector3.Dot(this.normal, ray.direction);
        if (denominator == 0) return false;

        var t = -(Vector3.Dot(ray.origin, this.normal) + Utils.PlaneConstant(this.normal, this.averagePosition)) / denominator;
        var point = (ray.direction * t) + ray.origin;
        var origin = Vector3.zero;
        for(var i = 0; i < this.corners.Length; ++i)
        {
            var j = (i + 1) % this.corners.Length;
            var sideNormal = Vector3.Cross((this.corners[i].position), this.corners[j].position);
            var sideOrigin = origin;

            if (Utils.SignedDistancePlanePoint(sideNormal, sideOrigin, point) < 0) return false;
        }
        return true;
    }

    public void DebugDraw(Color color, float duration, float offset)
    {
        foreach(var b in this.borders)
        {
            if (b.corners.Length != 2)
                continue;

            Debug.DrawLine(b.corners[0].position + this.normal * offset, b.corners[1].position + this.normal * offset, color, duration);
        }
    }
}
public class Topology
{
    public List<Corner> corners;
    public List<Border> borders;
    public List<Tile> tiles;
    public SpatialPartition partition;
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

    public SpatialPartition generatePlanetPartition()
    {
        var icosahedron =  Utils.generateIcosahedron();

        // TODO: Change these 1000s with a global scale;
        for (var i = 0; i < icosahedron.faces.Count; ++i)
        {
            var face = icosahedron.faces[i];
            var p0 = icosahedron.nodes[face.n[0]].p * 1000;
            var p1 = icosahedron.nodes[face.n[1]].p * 1000;
            var p2 = icosahedron.nodes[face.n[2]].p * 1000;
            var center = (p0 + p1 + p2) / 3;
            var radius = Math.Max(Vector3.Distance(center, p0), Math.Max(Vector3.Distance(center, p2), Vector3.Distance(center, p2)));
            face.sphereCenter = center;
            face.sphereRadius = radius;
            face.children = new List<Tile>();
        }

        var unparentedTiles = new List<Tile>();
        var maxDistanceFromOrigin = 0f;

        for(var i = 0; i < this.tiles.Count; ++i)
        {
			var tile = this.tiles[i];
			maxDistanceFromOrigin = Math.Max(maxDistanceFromOrigin, tile.averagePosition.magnitude + tile.radius);
			
			var parentFound = false;
			for (var j = 0; j < icosahedron.faces.Count; ++j)
			{
				var face = icosahedron.faces[j];
				var distance = Vector3.Distance(tile.averagePosition, face.sphereCenter) + tile.radius;
				if (distance < face.sphereRadius)
				{
					face.children.Add(tile);
					parentFound = true;
					break;
				}
			}
			if (!parentFound)
			{
				unparentedTiles.Add(tile);
			}
        }

        SpatialPartition rootPartition = new SpatialPartition(Vector3.zero, maxDistanceFromOrigin, new List<SpatialPartition>(), unparentedTiles);
        for(var i = 0; i < icosahedron.faces.Count; ++i)
        {
            var face = icosahedron.faces[i];
            rootPartition.partitions.Add(new SpatialPartition(face.sphereCenter, face.sphereRadius, new List<SpatialPartition>(), face.children));
        }

        this.partition = rootPartition;
        return rootPartition;
    }
}
public class SpatialPartition
{
    public Vector3 sphereCenter;
    public float sphereRadius;
    public List<Tile> tiles;
    public List<SpatialPartition> partitions;

    public SpatialPartition(Vector3 center, float radius, List<SpatialPartition> partitions, List<Tile> tiles)
    {
        this.sphereCenter = center;
        this.sphereRadius = radius;
        this.partitions = partitions;
        this.tiles = tiles;
    }

    public bool intersectRay(Ray ray, out Tile tile)
    {
        tile = null;
        if(Utils.intersectRayWithSphere(ray, this.sphereCenter, this.sphereRadius))
        {
            for(var i = 0; i < this.partitions.Count; ++i)
            {
                Tile result;
                var intersection = this.partitions[i].intersectRay(ray, out result);
                if(intersection)
                {
                    tile = result;
                    return true;
                }
            }

            for(var i = 0; i < this.tiles.Count; ++i)
            {
                if (this.tiles[i].intersectRay(ray))
                {
                    tile = this.tiles[i];
                    return true;
                }
            }
        }
        return false;
    }
}
public static class Utils
{
    public static bool intersectRayWithSphere(Ray ray, Vector3 sCenter, float sRadius)
    {
        var v1 = sCenter - ray.origin;
        var v2 = Vector3.Project(v1, ray.direction);
        var d = Vector3.Distance(v1, v2);
        return (d <= sRadius);
    }
    
    public static Quaternion randomQuaternion(Random random)
    {
        var theta = random.NextDouble() * Math.PI * 2; // real(0, Math.PI * 2);
        var phi = Math.Acos(random.NextDouble() * 2 - 1); // realInclusive(-1, 1));
        var sinPhi = Math.Sin(phi);
        var gamma = random.NextDouble() * Math.PI * 2; // real(0, Math.PI * 2);
        var sinGamma = Math.Sin(gamma);
        return new Quaternion(
            (float)(Math.Cos(theta) * sinPhi * sinGamma),
            (float)(Math.Sin(theta) * sinPhi * sinGamma),
            (float)(Math.Cos(phi) * sinGamma),
            (float)Math.Cos(gamma));
    }

    public static Vector3 randomUnitVector(Random random)
    {
        var theta = random.NextDouble() * Math.PI * 2; // real(0, Math.PI * 2);
        var phi = Math.Acos(random.NextDouble() * 2 - 1); // realInclusive(-1, 1));
        var sinPhi = Math.Sin(phi);
        return new Vector3(
            (float)(Math.Cos(theta) * sinPhi),
            (float)(Math.Sin(theta) * sinPhi),
            (float)Math.Cos(phi));
    }

    //This function returns a point which is a projection from a point to a plane.
    public static Vector3 ProjectPointOnPlane(Vector3 planeNormal, Vector3 planePoint, Vector3 point)
    {

        float distance;
        Vector3 translationVector;

        //First calculate the distance from the point to the plane:
        distance = SignedDistancePlanePoint(planeNormal, planePoint, point);

        //Reverse the sign of the distance
        distance *= -1;

        //Get a translation vector
        translationVector = SetVectorLength(planeNormal, distance);

        //Translate the point to form a projection
        return point + translationVector;
    }

    //create a vector of direction "vector" with length "size"
    public static Vector3 SetVectorLength(Vector3 vector, float size)
    {

        //normalize the vector
        Vector3 vectorNormalized = Vector3.Normalize(vector);

        //scale the vector
        return vectorNormalized *= size;
	}

    public static float calculateTriangleArea(Vector3 pa, Vector3 pb, Vector3 pc)
    {
        float A = Vector3.Distance(pa, pb);
        float B = Vector3.Distance(pb, pc);
        float C = Vector3.Distance(pc, pa);

        float s = (A + B + C) / 2;
        float perimeter = A + B + C;
        float area = (float) Math.Sqrt(s * (s - A) * (s - B) * (s - C));
        return area;
    }

    //Get the shortest distance between a point and a plane. The output is signed so it holds information
	//as to which side of the plane normal the point is.
	public static float SignedDistancePlanePoint(Vector3 planeNormal, Vector3 planePoint, Vector3 point){
 
		return Vector3.Dot(planeNormal, (point - planePoint));
	}

    //The negative distance from the origin to the plane along the normal vector.
    public static float PlaneConstant(Vector3 planeNormal, Vector3 planePoint)
    {
        // Make sure the given normal is a unit vector;
        Vector3 norm = planeNormal.normalized;
        float dist = Vector3.Dot(planePoint, norm);

        return (dist < 0) ? dist : -dist;
    }

    public static Polyhedra generateSubdividedIcosahedron(int degree)
    {
        var icosahedron = generateIcosahedron();

        List<Node> nodes = new List<Node>();
        for(int i = 0; i < icosahedron.nodes.Count; ++i)
        {
            nodes.Add(new Node(icosahedron.nodes[i].p));
        }

        List<Edge> edges = new List<Edge>();
        for (var i = 0; i < icosahedron.edges.Count; ++i)
        {
            var edge = icosahedron.edges[i];
            //edge.subdivided_n = new List<int>();
            //edge.subdivided_e = [];
            var n0 = icosahedron.nodes[edge.n[0]];
            var n1 = icosahedron.nodes[edge.n[1]];
            var p0 = n0.p;
            var p1 = n1.p;
            var delta = p1 - p0; //p1.clone().sub(p0);
            nodes[edge.n[0]].e.Add(edges.Count);
            var priorNodeIndex = edge.n[0];
            for (var s = 1; s < degree; ++s)
            {
                var edgeIndex = edges.Count;
                var nodeIndex = nodes.Count;
                edge.subdivided_e.Add(edgeIndex);
                edge.subdivided_n.Add(nodeIndex);
                edges.Add(new Edge(priorNodeIndex, nodeIndex)); // { n: [ priorNodeIndex, nodeIndex ], f: [] });
                priorNodeIndex = nodeIndex;
                Node newnode = new Node(Vector3.Slerp(p0, p1, (float)s / degree));
                newnode.e = new List<int>() { edgeIndex, edgeIndex + 1 };
                nodes.Add(newnode); // { p: slerp(p0, p1, s / degree), e: [ edgeIndex, edgeIndex + 1 ], f: [] });
            }
            edge.subdivided_e.Add(edges.Count);
            nodes[edge.n[1]].e.Add(edges.Count);
            edges.Add(new Edge(priorNodeIndex, edge.n[1])); // { n: [ priorNodeIndex, edge.n[1] ], f: [] });
        }

        List<Face> faces = new List<Face>();
        for (var i = 0; i < icosahedron.faces.Count; ++i)
        {
            var face = icosahedron.faces[i];
            var edge0 = icosahedron.edges[face.e[0]];
            var edge1 = icosahedron.edges[face.e[1]];
            var edge2 = icosahedron.edges[face.e[2]];
            var point0 = icosahedron.nodes[face.n[0]].p;
            var point1 = icosahedron.nodes[face.n[1]].p;
            var point2 = icosahedron.nodes[face.n[2]].p;
            var delta = point1 - point0; // point1.clone().sub(point0);


            var getEdgeNode0 = (face.n[0] == edge0.n[0])
                ? new Func<int, int>((k) => { return edge0.subdivided_n[k]; })
                : new Func<int, int>((k) => { return edge0.subdivided_n[degree - 2 - k]; });
            var getEdgeNode1 = (face.n[1] == edge1.n[0])
                ? new Func<int, int>((k) => { return edge1.subdivided_n[k]; })
                : new Func<int, int>((k) => { return edge1.subdivided_n[degree - 2 - k]; });
            var getEdgeNode2 = (face.n[0] == edge2.n[0])
                ? new Func<int, int>((k) => { return edge2.subdivided_n[k]; })
                : new Func<int, int>((k) => { return edge2.subdivided_n[degree - 2 - k]; });

            var faceNodes = new List<int>();
            faceNodes.Add(face.n[0]);
            for (var j = 0; j < edge0.subdivided_n.Count; ++j)
                faceNodes.Add(getEdgeNode0(j));
            faceNodes.Add(face.n[1]);
            for (var s = 1; s < degree; ++s)
            {
                faceNodes.Add(getEdgeNode2(s - 1));
                var p0 = nodes[getEdgeNode2(s - 1)].p;
                var p1 = nodes[getEdgeNode1(s - 1)].p;
                for (var t = 1; t < degree - s; ++t)
                {
                    faceNodes.Add(nodes.Count);
                    nodes.Add(new Node(Vector3.Slerp(p0, p1, (float)t / (degree - s)))); // { p: slerp(p0, p1, t / (degree - s)), e: [], f: [], });
                }
                faceNodes.Add(getEdgeNode1(s - 1));
            }
            faceNodes.Add(face.n[2]);

            var getEdgeEdge0 = (face.n[0] == edge0.n[0])
                ? new Func<int, int>((k) => { return edge0.subdivided_e[k]; })
                : new Func<int, int>((k) => { return edge0.subdivided_e[degree - 1 - k]; });
            var getEdgeEdge1 = (face.n[1] == edge1.n[0])
                ? new Func<int, int>((k) => { return edge1.subdivided_e[k]; })
                : new Func<int, int>((k) => { return edge1.subdivided_e[degree - 1 - k]; });
            var getEdgeEdge2 = (face.n[0] == edge2.n[0])
                ? new Func<int, int>((k) => { return edge2.subdivided_e[k]; })
                : new Func<int, int>((k) => { return edge2.subdivided_e[degree - 1 - k]; });

            var faceEdges0 = new List<int>();
            for (var j = 0; j < degree; ++j)
                faceEdges0.Add(getEdgeEdge0(j));
            var nodeIndex = degree + 1;
            for (var s = 1; s < degree; ++s)
            {
                for (var t = 0; t < degree - s; ++t)
                {
                    faceEdges0.Add(edges.Count);
                    var edge = new Edge(faceNodes[nodeIndex], faceNodes[nodeIndex + 1]);
                    nodes[edge.n[0]].e.Add(edges.Count);
                    nodes[edge.n[1]].e.Add(edges.Count);
                    edges.Add(edge);
                    ++nodeIndex;
                }
                ++nodeIndex;
            }

            var faceEdges1 = new List<int>();

            nodeIndex = 1;
            for (var s = 0; s < degree; ++s)
            {
                for (var t = 1; t < degree - s; ++t)
                {
                    faceEdges1.Add(edges.Count);
                    var edge = new Edge(faceNodes[nodeIndex], faceNodes[nodeIndex + degree - s]);
                    nodes[edge.n[0]].e.Add(edges.Count);
                    nodes[edge.n[1]].e.Add(edges.Count);
                    edges.Add(edge);
                    ++nodeIndex;
                }
                faceEdges1.Add(getEdgeEdge1(s));
                nodeIndex += 2;
            }

            var faceEdges2 = new List<int>();
            nodeIndex = 1;
            for (var s = 0; s < degree; ++s)
            {
                faceEdges2.Add(getEdgeEdge2(s));
                for (var t = 1; t < degree - s; ++t)
                {
                    faceEdges2.Add(edges.Count);
                    var edge = new Edge(faceNodes[nodeIndex], faceNodes[nodeIndex + degree - s + 1]);
                    nodes[edge.n[0]].e.Add(edges.Count);
                    nodes[edge.n[1]].e.Add(edges.Count);
                    edges.Add(edge);
                    ++nodeIndex;
                }
                nodeIndex += 2;
            }

            nodeIndex = 0;
            var edgeIndex = 0;
            for (var s = 0; s < degree; ++s)
            {
                for (var t = 1; t < degree - s + 1; ++t)
                {
                    var subFace = new Face(new int[] { faceNodes[nodeIndex], faceNodes[nodeIndex + 1], faceNodes[nodeIndex + degree - s + 1] }
                        , new int[] { faceEdges0[edgeIndex], faceEdges1[edgeIndex], faceEdges2[edgeIndex] });
                    //{n: [ faceNodes[nodeIndex], faceNodes[nodeIndex + 1], faceNodes[nodeIndex + degree - s + 1], ],
                    //e: [ faceEdges0[edgeIndex], faceEdges1[edgeIndex], faceEdges2[edgeIndex], ], };
                    nodes[subFace.n[0]].f.Add(faces.Count);
                    nodes[subFace.n[1]].f.Add(faces.Count);
                    nodes[subFace.n[2]].f.Add(faces.Count);
                    edges[subFace.e[0]].f.Add(faces.Count);
                    edges[subFace.e[1]].f.Add(faces.Count);
                    edges[subFace.e[2]].f.Add(faces.Count);
                    faces.Add(subFace);
                    ++nodeIndex;
                    ++edgeIndex;
                }
                ++nodeIndex;
            }

            nodeIndex = 1;
            edgeIndex = 0;
            for (var s = 1; s < degree; ++s)
            {
                for (var t = 1; t < degree - s + 1; ++t)
                {
                    var subFace = new Face(new int[] { faceNodes[nodeIndex], faceNodes[nodeIndex + degree - s + 2], faceNodes[nodeIndex + degree - s + 1] }
                                          , new int[] { faceEdges2[edgeIndex + 1], faceEdges0[edgeIndex + degree - s + 1], faceEdges1[edgeIndex] });
                    //n: [ faceNodes[nodeIndex], faceNodes[nodeIndex + degree - s + 2], faceNodes[nodeIndex + degree - s + 1], ],
                    //e: [ faceEdges2[edgeIndex + 1], faceEdges0[edgeIndex + degree - s + 1], faceEdges1[edgeIndex], ], };
                    nodes[subFace.n[0]].f.Add(faces.Count);
                    nodes[subFace.n[1]].f.Add(faces.Count);
                    nodes[subFace.n[2]].f.Add(faces.Count);
                    edges[subFace.e[0]].f.Add(faces.Count);
                    edges[subFace.e[1]].f.Add(faces.Count);
                    edges[subFace.e[2]].f.Add(faces.Count);
                    faces.Add(subFace);
                    ++nodeIndex;
                    ++edgeIndex;
                }
                nodeIndex += 2;
                edgeIndex += 1;
            }
        }
        return new Polyhedra(nodes, edges, faces);
    }

    public static Polyhedra generateIcosahedron()
    {
        var phi = (float) ((1.0 + Math.Sqrt(5.0)) / 2.0);
        float du =(float) (1.0 / Math.Sqrt(phi * phi + 1.0));
        float dv =(float) (phi * du);

        //Vector3[] vertices = new Vector3[] { new Vector3(0, +dv, +du),
        //                                     new Vector3(0, +dv, -du),
        //                                     new Vector3(0, -dv, +du),
        //                                     new Vector3(0, -dv, -du),
        //                                     new Vector3(+du, 0, +dv),
        //                                     new Vector3(-du, 0, +dv),
        //                                     new Vector3(+du, 0, -dv),
        //                                     new Vector3(-du, 0, -dv),
        //                                     new Vector3(+dv, +du, 0),
        //                                     new Vector3(+dv, -du, 0),
        //                                     new Vector3(-dv, +du, 0),
        //                                     new Vector3(-dv, -du, 0)};
        //int[] tris = new int[] { 0, 1, 8,
        //                         0, 4, 5,
        //                         0, 5, 10,
        //                         0, 8, 4,
        //                         0, 10, 1,
        //                         1, 6, 8,
        //                         1, 7, 6,
        //                         1, 10, 7,
        //                         2, 3, 11,
        //                         2, 4, 9,
        //                         2, 5, 4,
        //                         2, 9, 3,
        //                         2, 11, 5,
        //                         3, 6, 7,
        //                         3, 7, 11,
        //                         3, 9, 6,
        //                         4, 8, 9,
        //                         5, 11, 10,
        //                         6, 9, 8,
        //                         7, 10, 11};
        //Mesh icosahedron = new Mesh();
        //icosahedron.vertices = vertices;
        //icosahedron.SetTriangles(tris, 0);
        //icosahedron.RecalculateNormals();
        //icosahedron.RecalculateBounds();

        var nodes = new List<Node>() {
            new Node(new Vector3(0, +dv, +du)),
            new Node(new Vector3(0, +dv, -du)),
            new Node(new Vector3(0, -dv, +du)),
            new Node(new Vector3(0, -dv, -du)),
            new Node(new Vector3(+du, 0, +dv)),
            new Node(new Vector3(-du, 0, +dv)),
            new Node(new Vector3(+du, 0, -dv)),
            new Node(new Vector3(-du, 0, -dv)),
            new Node(new Vector3(+dv, +du, 0)),
            new Node(new Vector3(+dv, -du, 0)),
            new Node(new Vector3(-dv, +du, 0)),
            new Node(new Vector3(-dv, -du, 0))};

        var edges = new List<Edge>() {
            new Edge(0,  1),
            new Edge(0,  4),
            new Edge(0,  5),
            new Edge(0,  8),
            new Edge(0, 10),
            new Edge(1,  6),
            new Edge(1,  7),
            new Edge(1,  8),
            new Edge(1, 10),
            new Edge(2,  3),
            new Edge(2,  4),
            new Edge(2,  5),
            new Edge(2,  9),
            new Edge(2, 11),
            new Edge(3,  6),
            new Edge(3,  7),
            new Edge(3,  9),
            new Edge(3, 11),
            new Edge(4,  5),
            new Edge(4,  8),
            new Edge(4,  9),
            new Edge(5, 10),
            new Edge(5, 11),
            new Edge(6,  7),
            new Edge(6,  8),
            new Edge(6,  9),
            new Edge(7, 10),
            new Edge(7, 11),
            new Edge(8,  9),
            new Edge(10, 11)};

        var faces = new List<Face>() {
            new Face(new int[] {0,  1,  8 },new int[] { 0,  7,  3 }),
            new Face(new int[] {0,  4,  5 },new int[] { 1, 18,  2 }),
            new Face(new int[] {0,  5, 10 },new int[] { 2, 21,  4 }),
            new Face(new int[] {0,  8,  4 },new int[] { 3, 19,  1 }),
            new Face(new int[] {0, 10,  1 },new int[] { 4,  8,  0 }),
            new Face(new int[] {1,  6,  8 },new int[] { 5, 24,  7 }),
            new Face(new int[] {1,  7,  6 },new int[] { 6, 23,  5 }),
            new Face(new int[] {1, 10,  7 },new int[] { 8, 26,  6 }),
            new Face(new int[] {2,  3, 11 },new int[] { 9, 17, 13 }),
            new Face(new int[] {2,  4,  9 },new int[] {10, 20, 12 }),
            new Face(new int[] {2,  5,  4 },new int[] {11, 18, 10 }),
            new Face(new int[] {2,  9,  3 },new int[] {12, 16,  9 }),
            new Face(new int[] {2, 11,  5 },new int[] {13, 22, 11 }),
            new Face(new int[] {3,  6,  7 },new int[] {14, 23, 15 }),
            new Face(new int[] {3,  7, 11 },new int[] {15, 27, 17 }),
            new Face(new int[] {3,  9,  6 },new int[] {16, 25, 14 }),
            new Face(new int[] {4,  8,  9 },new int[] {19, 28, 20 }),
            new Face(new int[] {5, 11, 10 },new int[] {22, 29, 21 }),
            new Face(new int[] {6,  9,  8 },new int[] {25, 28, 24 }),
            new Face(new int[] {7, 10, 11 },new int[] {26, 29, 27 })
        };

        for (var i = 0; i < edges.Count; ++i)
            for (var j = 0; j < edges[i].n.Length; ++j)
                nodes[j].e.Add(i);

        for (var i = 0; i < faces.Count; ++i)
            for (var j = 0; j < faces[i].n.Count; ++j)
                nodes[j].f.Add(i);

        for (var i = 0; i < faces.Count; ++i)
            for (var j = 0; j < faces[i].e.Count; ++j)
                edges[j].f.Add(i);


        return new Polyhedra(nodes, edges, faces);
    }
}

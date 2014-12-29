using UnityEngine;
using System.Collections;
using System.Collections.Generic;
using System;
using Random = System.Random;

public class Planet : MonoBehaviour {
    private Polyhedra planet;
    private Random globalRandom;

    public int randomSeed = 5;
    public int subdivisionLevel = 2;
    public float distortionLevel = 1f;
    int plateCount = 36;
    float oceanicRate = 0.7f;
    float heatLevel = 1.0f;
    float moistureLevel = 1.0f;

    void Awake()
    {
        //Generate();
    }

    public void Generate()
    {
        this.globalRandom = new Random(randomSeed);
        this.planet = GeneratePlanetPolyhedra();
        this.planet.GeneratePlanetTopology();
        this.planet.topology.GeneratePlanetPartition();
        this.GeneratePlanetMesh();
    }

    void Update()
    {
        Tile tile;
        this.IntersectRay(Camera.main.ScreenPointToRay(Input.mousePosition), out tile);
        if (tile != null)
        {
            tile.DebugDraw(Color.green, 0, this.transform.localScale.x / 100f, this.transform);
        }
    }

    private Polyhedra GeneratePlanetPolyhedra()
    {
        var mesh = Utils.GenerateSubdividedIcosahedron(this.subdivisionLevel);
        mesh.planet = this;
        var totalDistortion = Math.Ceiling(mesh.edges.Count * this.distortionLevel);
		var remainingIterations = 6;
        int i;
        while(remainingIterations > 0)
        {
			var iterationDistortion = (int)Math.Floor(totalDistortion / remainingIterations);
			totalDistortion -= iterationDistortion;
            mesh.DistortMesh(iterationDistortion, globalRandom);
			mesh.RelaxMesh(0.5f);
			--remainingIterations;
        }
        
        //var initialIntervalIteration = action.intervalIteration;
	
		var averageNodeRadius = Math.Sqrt(4 * Math.PI / mesh.nodes.Count);
		var minShiftDelta = averageNodeRadius / 50000 * mesh.nodes.Count;
		var maxShiftDelta = averageNodeRadius / 50 * mesh.nodes.Count;

		float priorShift, shiftDelta;
		var currentShift = mesh.RelaxMesh(0.5f);

        do
        {
            priorShift = currentShift;
            currentShift = mesh.RelaxMesh(0.5f);
            shiftDelta = Math.Abs(currentShift - priorShift);
        } while (shiftDelta >= minShiftDelta /* && action.intervalIteration - initialIntervalIteration < 300*/);

        for(i = 0; i < mesh.faces.Count; ++i)
        {
            var face = mesh.faces[i];
            var p0 = mesh.nodes[face.n[0]].p;
            var p1 = mesh.nodes[face.n[1]].p;
            var p2 = mesh.nodes[face.n[2]].p;
            face.centroid = Utils.CalculateFaceCentroid(p0, p1, p2).normalized;
        }

        for(i = 0; i < mesh.nodes.Count; ++i)
        {
            var node = mesh.nodes[i];
            var faceIndex = node.f[0];
            for(var j = 1; j < node.f.Count - 1; ++j)
            {
                faceIndex = Utils.FindNextFaceIndex(mesh, i, faceIndex);
                var k = node.f.IndexOf(faceIndex);
                node.f[k] = node.f[j];
                node.f[j] = faceIndex;
            }
        }
        this.planet = mesh;
        return mesh;
    }
    private void GeneratePlanetMesh()
    {
        Matrix4x4 trs = Matrix4x4.TRS(transform.position, transform.rotation, transform.localScale);
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
            verts.Add(trs.MultiplyPoint3x4(tile.averagePosition));
            norms.Add((transform.rotation * tile.averagePosition).normalized);

            //Debug.DrawLine(verts[verts.Count - 1], verts[verts.Count - 1] + norms[norms.Count - 1] * 25f, Color.red, 1000f);

	        tileStart = verts.Count - 1;
	        foreach (var corner in tile.corners)
	        {
	            verts.Add(trs.MultiplyPoint3x4(corner.position));
                norms.Add((transform.rotation * tile.averagePosition).normalized);

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
    }
    
    public bool IntersectRay(Ray ray, out Tile tile)
    {
        SpatialPartition partition;
        return this.IntersectRay(ray, out tile, out partition);
    }
    public bool IntersectRay(Ray ray, out Tile tile, out SpatialPartition partition)
    {
        tile = null;
        partition = null;

        if(planet.topology.partition.IntersectRay(ray, out tile, this.transform))
        {
            return true;
        }
        //foreach (var p in planet.topology.partition.partitions)
        //{
        //    if (p.IntersectRay(ray, out tile, this.transform))
        //    {
        //        partition = p;
        //        tile.DebugDraw(Color.red, 0.001f, 2f, this.transform);
        //        foreach (var t in p.tiles)
        //        {
        //            t.DebugDraw(Color.green, 0, 1f, this.transform);
        //        }
        //        return true;
        //    }
        //}
        return false;
    }
	// Use this for initialization
}

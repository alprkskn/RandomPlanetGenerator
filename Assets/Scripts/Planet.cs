using UnityEngine;
using System.Collections;
using System.Collections.Generic;
using System;
using Random = System.Random;

public class Planet : MonoBehaviour {

    private Polyhedra surface;
    private Random globalRandom;
    private List<Plate> tectonicPlates;

    public int randomSeed = 5;
    public int subdivisionLevel = 2;
    public float distortionLevel = 1f;
    public int plateCount = 36;
    public float oceanicRate = 0.7f;
    public float heatLevel = 1.0f;
    public float moistureLevel = 1.0f;

    void Awake()
    {
        //Generate();
    }

    public void Generate()
    {
        this.globalRandom = new Random(randomSeed);
        this.surface = GeneratePlanetPolyhedra();
        this.surface.GeneratePlanetTopology();
        this.surface.topology.GeneratePlanetPartition();
        this.GeneratePlanetTectonicPlates();

        //DrawTectonicPlates(6000f);
        this.GeneratePlanetMesh();
    }

    void Update()
    {
        Tile tile;
        this.IntersectRay(Camera.main.ScreenPointToRay(Input.mousePosition), out tile);
        if (tile != null)
        {
            tile.DebugDraw(Color.red, 0, this.transform.localScale.x / 100f, this.transform);
            foreach(var t in tile.plate.tiles)
            {
                t.DebugDraw(Color.green, 0, this.transform.localScale.x / 150f, this.transform);
            }
        }

    }

    void DrawTectonicPlates(float duration)
    {
        foreach(var p in this.tectonicPlates)
        {
            Color c = new Color((float)globalRandom.NextDouble(), (float)globalRandom.NextDouble(), (float)globalRandom.NextDouble());
            foreach(var t in p.tiles)
            {
                t.DebugDraw(c, duration, this.transform.localScale.x / 150f, this.transform);
            }
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
        this.surface = mesh;
        return mesh;
    }
    private void GeneratePlanetMesh()
    {
        Matrix4x4 trs = Matrix4x4.TRS(transform.position, transform.rotation, transform.localScale);
	    var verts = new List<Vector3>();
        var norms = new List<Vector3>();
        var tris = new List<int>();

	    var tileStart = 0;
	    foreach (var tile in surface.topology.tiles)
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

	        tileStart = verts.Count - 1;
	        foreach (var corner in tile.corners)
	        {
	            verts.Add(trs.MultiplyPoint3x4(corner.position));
                norms.Add((transform.rotation * tile.averagePosition).normalized);
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
    
    private void GeneratePlanetTerrain()
    {

    }

    #region Tectonic Plate Generation
    private List<Plate> GeneratePlanetTectonicPlates()
    {
        var plates = new List<Plate>();
        var platelessTiles = new List<Tile>();
        var platelessTilePlates = new List<Plate>();
        var topology = this.surface.topology;

        var failedCount = 0;
        while(plates.Count < plateCount && failedCount < 10000)
        {
            var corner = topology.corners[globalRandom.Next(0, topology.corners.Count)];
            var adjacentToExistingPlate = false;
            for(var i = 0; i < corner.tiles.Length; ++i)
            {
                if(corner.tiles[i].plate != null)
                {
                    adjacentToExistingPlate = true;
                    failedCount += 1;
                    break;
                }
            }
            if (adjacentToExistingPlate) continue;

            failedCount = 0;

            var oceanic = (globalRandom.NextDouble() < this.oceanicRate);
            var plate = new Plate(new Color((float)globalRandom.NextDouble(), (float)globalRandom.NextDouble(), (float)globalRandom.NextDouble()),
                                    Utils.RandomUnitVector(globalRandom),
                                    (float) (globalRandom.NextDouble() * (Math.PI / 15f) - (Math.PI / 30f)),
                                    (float) (globalRandom.NextDouble() * (Math.PI / 15f) - (Math.PI / 30f)),
                                    (float) ((oceanic) ? (globalRandom.NextDouble() * 0.5 - 0.8) : (globalRandom.NextDouble() * 0.4 + 0.1)),
                                    oceanic,
                                    corner);
            plates.Add(plate);

            for(var i = 0; i < corner.tiles.Length; ++i)
            {
                corner.tiles[i].plate = plate;
                plate.tiles.Add(corner.tiles[i]);
            }
            for(var i = 0; i < corner.tiles.Length; ++i)
            {
                var tile = corner.tiles[i];
                for(var j = 0; j < tile.tiles.Length; ++j)
                {
                    var adjacentTile = tile.tiles[j];
                    if(adjacentTile.plate == null)
                    {
                        platelessTiles.Add(adjacentTile);
                        platelessTilePlates.Add(plate);
                    }
                }
            }
        }

        while(platelessTiles.Count > 0)
        {
            // XXX: Try something else instead of globalRandom.Next() if this fails.
            var tileIndex = (int)Math.Floor(Math.Pow(globalRandom.NextDouble(), 2) * platelessTiles.Count);
			var tile = platelessTiles[tileIndex];
			var plate = platelessTilePlates[tileIndex];
            platelessTiles.RemoveAt(tileIndex); // splice(tileIndex, 1);
            platelessTilePlates.RemoveAt(tileIndex); // splice(tileIndex, 1);
			if (tile.plate == null)
			{
				tile.plate = plate;
				plate.tiles.Add(tile);
				for (var j = 0; j < tile.tiles.Length; ++j)
				{
					if (tile.tiles[j].plate == null)
					{
						platelessTiles.Add(tile.tiles[j]);
						platelessTilePlates.Add(plate);
					}
				}
			}
        }
        calculateCornerDistancesToPlateRoot(plates);
        this.tectonicPlates = plates;
        return plates;
    }

    struct DistanceCorner
    {
        public Corner corner;
        public float distanceToPlateRoot;
        public DistanceCorner(Corner c, float d)
        {
            this.corner = c;
            this.distanceToPlateRoot = d;
        }
    }

    void calculateCornerDistancesToPlateRoot(List<Plate> plates)
    {
        var distanceCornerQueue = new List<DistanceCorner>();
        for(var i = 0; i < plates.Count; ++i)
        {
            var corner = plates[i].root;
            corner.distanceToPlateRoot = 0;
            for(var j = 0; j < corner.corners.Length; ++j)
            {
                distanceCornerQueue.Add(new DistanceCorner(corner.corners[j], corner.borders[j].Length()));
            }
        }

        while(true)
        {
            if (distanceCornerQueue.Count == 0) return;

            var iEnd = distanceCornerQueue.Count;
            for(var i = 0; i < iEnd; ++i)
            {
                var front = distanceCornerQueue[i];
                var corner = front.corner;
                var distanceToPlateRoot = front.distanceToPlateRoot;
                if(corner.distanceToPlateRoot == 0 || corner.distanceToPlateRoot > distanceToPlateRoot)
                {
                    corner.distanceToPlateRoot = distanceToPlateRoot;
                    for(var j = 0; j < corner.corners.Length; ++j)
                    {
                        distanceCornerQueue.Add(new DistanceCorner(corner.corners[j], distanceToPlateRoot + corner.borders[j].Length()));
                    }
                }
            }

            distanceCornerQueue.RemoveRange(0, iEnd);
            distanceCornerQueue.Sort(delegate(DistanceCorner A, DistanceCorner B)
            {
                return A.distanceToPlateRoot.CompareTo(B.distanceToPlateRoot);
            });
        }
    }
    #endregion

    #region Intersection/Tile Selection
    public bool IntersectRay(Ray ray, out Tile tile)
    {
        List<Tile> partition;
        return this.IntersectRay(ray, out tile, out partition);
    }
    public bool IntersectRay(Ray ray, out Tile tile, out List<Tile> partition)
    {
        tile = null;
        partition = null;
        if(surface.topology.partition.IntersectRay(ray, out tile, out partition, this.transform))
        {
            return true;
        }
        return false;
    }
    #endregion
}

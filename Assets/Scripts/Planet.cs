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
        this.GeneratePlanetElevation();
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
    #region PlanetElevation
    void GeneratePlanetElevation()
    {
        IdentifyBoundaryBorders(this.surface.topology.borders);
        var boundaryCorners = CollectBoundaryCorners(this.surface.topology.corners);
        var boundaryCornerInnerBorderIndices = CalculatePlateBoundaryStress(boundaryCorners);
        BlurPlateBoundaryStress(boundaryCorners, 3, 0.4f);
        var elevationBorderQueue = PopulateElevationBorderQueue(boundaryCorners, boundaryCornerInnerBorderIndices);
        ProcessElevationBorderQueue(elevationBorderQueue);
        CalculateTileAverageElevations(this.surface.topology.tiles);
    }
    void IdentifyBoundaryBorders(List<Border> borders)
    {
        for (var i = 0; i < borders.Count; ++i)
        {
            var border = borders[i];
            if (border.tiles[0].plate != border.tiles[1].plate)
            {
                border.betweenPlates = true;
                border.corners[0].betweenPlates = true;
                border.corners[1].betweenPlates = true;
                border.tiles[0].plate.boundaryBorders.Add(border);
                border.tiles[1].plate.boundaryBorders.Add(border);
            }
        }
    }
    List<Corner> CollectBoundaryCorners(List<Corner> corners)
    {
        var boundaryCorners = new List<Corner>();
        for (var j = 0; j < corners.Count; ++j)
        {
            var corner = corners[j];
            if (corner.betweenPlates)
            {
                boundaryCorners.Add(corner);
                corner.tiles[0].plate.boundaryCorners.Add(corner);
                if (corner.tiles[1].plate != corner.tiles[0].plate) corner.tiles[1].plate.boundaryCorners.Add(corner);
                if (corner.tiles[2].plate != corner.tiles[0].plate && corner.tiles[2].plate != corner.tiles[1].plate) corner.tiles[2].plate.boundaryCorners.Add(corner);
            }
        }
        return boundaryCorners;
    }
    int[] CalculatePlateBoundaryStress(List<Corner> boundaryCorners)
    {
        var boundaryCornerInnerBorderIndices = new int[boundaryCorners.Count];
        for(var i = 0; i < boundaryCorners.Count; ++i)
        {
            var corner = boundaryCorners[i];
            corner.distanceToPlateRoot = 0;

            Border innerBorder = null;
            int innerBorderIndex = 0;
            for(var j = 0; j < corner.borders.Length; ++j)
            {
                var border = corner.borders[j];
                if(!border.betweenPlates)
                {
                    innerBorder = border;
                    innerBorderIndex = j;
                    break;
                }
            }

            if(innerBorder != null)
            {
                boundaryCornerInnerBorderIndices[i] = innerBorderIndex;
                var outerBorder0 = corner.borders[(innerBorderIndex + 1) % corner.borders.Length];
                var outerBorder1 = corner.borders[(innerBorderIndex + 2) % corner.borders.Length];
                var farCorner0 = outerBorder0.OppositeCorner(corner);
                var farCorner1 = outerBorder1.OppositeCorner(corner);
                var plate0 = innerBorder.tiles[0].plate;
                var plate1 = outerBorder0.tiles[0].plate != plate0 ? outerBorder0.tiles[0].plate : outerBorder0.tiles[1].plate;
                var boundaryVector = farCorner0.vectorTo(farCorner1);
                var boundaryNormal = Vector3.Cross(boundaryVector, corner.position);
                var stress = CalculateStress(plate0.CalculateMovement(corner.position), plate1.CalculateMovement(corner.position), boundaryVector, boundaryNormal);
                corner.pressure = stress.pressure;
                corner.shear = stress.shear;
            }
            else
            {
                boundaryCornerInnerBorderIndices[i] = -1;
                var plate0 = corner.tiles[0].plate;
                var plate1 = corner.tiles[1].plate;
                var plate2 = corner.tiles[2].plate;
                var boundaryVector0 = corner.corners[0].vectorTo(corner);
                var boundaryVector1 = corner.corners[1].vectorTo(corner);
                var boundaryVector2 = corner.corners[2].vectorTo(corner);
                var boundaryNormal0 = Vector3.Cross(boundaryVector0, corner.position);
                var boundaryNormal1 = Vector3.Cross(boundaryVector1, corner.position);
                var boundaryNormal2 = Vector3.Cross(boundaryVector2, corner.position);
                var stress0 = CalculateStress(plate0.CalculateMovement(corner.position), plate1.CalculateMovement(corner.position), boundaryVector0, boundaryNormal0);
                var stress1 = CalculateStress(plate1.CalculateMovement(corner.position), plate2.CalculateMovement(corner.position), boundaryVector1, boundaryNormal1);
                var stress2 = CalculateStress(plate2.CalculateMovement(corner.position), plate0.CalculateMovement(corner.position), boundaryVector2, boundaryNormal2);

                corner.pressure = (stress0.pressure + stress1.pressure + stress2.pressure) / 3f;
                corner.shear = (stress0.shear + stress1.shear + stress2.shear) / 3f;
            }
        }

        return boundaryCornerInnerBorderIndices;
    }
    struct Stress
    {
        public float pressure;
        public float shear;
        public Stress(float p, float s)
        {
            this.pressure = p;
            this.shear = s;
        }
    }
    Stress CalculateStress(Vector3 movement0, Vector3 movement1, Vector3 boundaryVector, Vector3 boundaryNormal)
    {
        var relativeMovement = movement0 - movement1;
        var pressureVector = Vector3.Project(relativeMovement, boundaryNormal);
        var pressure = pressureVector.magnitude;
        if (Vector3.Dot(pressureVector, boundaryNormal) > 0) pressure = -pressure;
        var shear = Vector3.Project(relativeMovement, boundaryVector).magnitude;
        return new Stress(2f / (1f + (float)Math.Exp(-pressure / 30f)) - 1f, 2f / (1f + (float)Math.Exp(-shear / 30f)) - 1f);
    }
    void BlurPlateBoundaryStress(List<Corner> boundaryCorners, int stressBlurIterations, float stressBlurCenterWeighting)
    {
        var newCornerPressure = new float[boundaryCorners.Count];
        var newCornerShear = new float[boundaryCorners.Count];

        for (var i = 0; i < stressBlurIterations; ++i)
        {
            for (var j = 0; j < boundaryCorners.Count; ++j)
            {
                var corner = boundaryCorners[j];
                var averagePressure = 0f;
                var averageShear = 0f;
                var neighborCount = 0;
                for (var k = 0; k < corner.corners.Length; ++k)
                {
                    var neighbor = corner.corners[k];
                    if (neighbor.betweenPlates)
                    {
                        averagePressure += neighbor.pressure;
                        averageShear += neighbor.shear;
                        ++neighborCount;
                    }
                }
                newCornerPressure[j] = corner.pressure * stressBlurCenterWeighting + (averagePressure / neighborCount) * (1 - stressBlurCenterWeighting);
                newCornerShear[j] = corner.shear * stressBlurCenterWeighting + (averageShear / neighborCount) * (1 - stressBlurCenterWeighting);
            }

            for (var j = 0; j < boundaryCorners.Count; ++j)
            {
                var corner = boundaryCorners[j];
                if (corner.betweenPlates)
                {
                    corner.pressure = newCornerPressure[j];
                    corner.shear = newCornerShear[j];
                }
            }
        }
    }
    struct ElevationBorder
    {
        public struct Origin
        {
            public Corner corner;
            public float pressure;
            public float shear;
            public Plate plate;
            public Func<float[], float> calculateElevation;

            public Origin(Corner c, Plate p, Func<float[], float> f)
            {
                this.corner = c;
                this.pressure = c.pressure;
                this.shear = c.shear;
                this.plate = p;
                this.calculateElevation = f;
            }
        }

        public Origin origin;
        public Border border;
        public Corner corner;
        public Corner nextCorner;
        public float distanceToPlateBoundary;

        public ElevationBorder(Origin o, Border b, Corner c, Corner n, float d)
        {
            this.origin = o;
            this.border = b;
            this.corner = c;
            this.nextCorner = n;
            this.distanceToPlateBoundary = d;
        }
    }
    List<ElevationBorder> PopulateElevationBorderQueue(List<Corner> boundaryCorners, int[] boundaryCornerInnerBorderIndices)
    {
        var elevationBorderQueue = new List<ElevationBorder>();
        for (var i = 0; i < boundaryCorners.Count; ++i)
        {
            var corner = boundaryCorners[i];

            var innerBorderIndex = boundaryCornerInnerBorderIndices[i];
            if (innerBorderIndex != -1)
            {
                var innerBorder = corner.borders[innerBorderIndex];
                var outerBorder0 = corner.borders[(innerBorderIndex + 1) % corner.borders.Length];
                var plate0 = innerBorder.tiles[0].plate;
                var plate1 = outerBorder0.tiles[0].plate != plate0 ? outerBorder0.tiles[0].plate : outerBorder0.tiles[1].plate;

                Func<float[], float> calculateElevation;

                if (corner.pressure > 0.3)
                {
                    corner.elevation = Math.Max(plate0.elevation, plate1.elevation) + corner.pressure;
                    if (plate0.oceanic == plate1.oceanic)
                        calculateElevation = CalculateCollidingElevation;
                    else if (plate0.oceanic)
                        calculateElevation = CalculateSubductingElevation;
                    else
                        calculateElevation = CalculateSuperductingElevation;
                }
                else if (corner.pressure < -0.3)
                {
                    corner.elevation = Math.Max(plate0.elevation, plate1.elevation) - corner.pressure / 4;
                    calculateElevation = CalculateDivergingElevation;
                }
                else if (corner.shear > 0.3)
                {
                    corner.elevation = Math.Max(plate0.elevation, plate1.elevation) + corner.shear / 8;
                    calculateElevation = CalculateShearingElevation;
                }
                else
                {
                    corner.elevation = (plate0.elevation + plate1.elevation) / 2;
                    calculateElevation = CalculateDormantElevation;
                }

                var nextCorner = innerBorder.OppositeCorner(corner);
                if (!nextCorner.betweenPlates)
                {
                    elevationBorderQueue.Add(new ElevationBorder(new ElevationBorder.Origin(corner, plate0, calculateElevation),
                                                                    innerBorder, corner, nextCorner, innerBorder.Length()));
                }
            }
            else
            {
                var plate0 = corner.tiles[0].plate;
                var plate1 = corner.tiles[1].plate;
                var plate2 = corner.tiles[2].plate;

                //elevation = 0;

                if (corner.pressure > 0.3)
                {
                    corner.elevation = Math.Max(plate0.elevation, Math.Max(plate1.elevation, plate2.elevation)) + corner.pressure;
                }
                else if (corner.pressure < -0.3)
                {
                    corner.elevation = Math.Max(plate0.elevation, Math.Max(plate1.elevation, plate2.elevation)) + corner.pressure / 4f;
                }
                else if (corner.shear > 0.3)
                {
                    corner.elevation = Math.Max(plate0.elevation, Math.Max(plate1.elevation, plate2.elevation)) + corner.shear / 8f;
                }
                else
                {
                    corner.elevation = (plate0.elevation + plate1.elevation + plate2.elevation) / 3;
                }
            }
        }
        return elevationBorderQueue;
    }
    // For the following elevation methods i had to use an array
    // as argument instead of 6 floats for some implementation issues.
    // Func<> implementation did not seem to support functions with more
    // than 4 arguments. So the array contains the following in order:
    // distanceToPlateBoundary, distanceToPlateRoot, boundaryElevation, plateElevation, pressure, shear
    float CalculateCollidingElevation(float[] args)
    {
        var t = args[0] / (args[0] + args[1]);
        if (t < 0.5f)
        {
            t = t / 0.5f;
            return args[3] + (float)Math.Pow(t - 1, 2) * (args[2] - args[3]);
        }
        else
        {
            return args[3];
        }
    }
    float CalculateSuperductingElevation(float[] args)
    {
        var t = args[0] / (args[0] + args[1]);
        if (t < 0.2f)
        {
            t = t / 0.2f;
            return args[2] + t * (args[3] - args[2] + args[4] / 2f);
        }
        else if (t < 0.5f)
        {
            t = (t - 0.2f) / 0.3f;
            return args[3] + (float)Math.Pow(t - 1, 2) * args[4] / 2f;
        }
        else
        {
            return args[3];
        }
    }
    float CalculateSubductingElevation(float[] args)
    {
        var t = args[0] / (args[0] + args[1]);
        return args[3] + (float)Math.Pow(t - 1, 2) * (args[2] - args[3]);
    }
    float CalculateDivergingElevation(float[] args)
    {
        var t = args[0] / (args[0] + args[1]);
        if (t < 0.3f)
        {
            t = t / 0.3f;
            return args[3] + (float)Math.Pow(t - 1, 2) * (args[2] - args[3]);
        }
        else
        {
            return args[3];
        }
    }
    float CalculateShearingElevation(float[] args)
    {
        var t = args[0] / (args[0] + args[1]);
        if (t < 0.2f)
        {
            t = t / 0.2f;
            return args[3] + (float)Math.Pow(t - 1, 2) * (args[2] - args[3]);
        }
        else
        {
            return args[3];
        }
    }
    float CalculateDormantElevation(float[] args)
    {
        var t = args[0] / (args[0] + args[1]);
        var elevationDifference = args[2] - args[3];
        var a = 2 * elevationDifference;
        var b = -3 * elevationDifference;
        return t * t * elevationDifference * (2 * t - 3) + args[2];
    }
    void ProcessElevationBorderQueue(List<ElevationBorder> elevationBorderQueue)
    {
        while (elevationBorderQueue.Count > 0)
        {
            var iEnd = elevationBorderQueue.Count;
            for (var i = 0; i < iEnd; ++i)
            {
                var front = elevationBorderQueue[i];
                var corner = front.nextCorner;
                if (corner.elevation == 0.0f)
                {
                    corner.distanceToPlateBoundary = front.distanceToPlateBoundary;
                    corner.elevation = front.origin.calculateElevation(new float[6] {
				corner.distanceToPlateBoundary,
				corner.distanceToPlateRoot,
				front.origin.corner.elevation,
				front.origin.plate.elevation,
				front.origin.pressure,
				front.origin.shear});

                    for (var j = 0; j < corner.borders.Length; ++j)
                    {
                        var border = corner.borders[j];
                        if (!border.betweenPlates)
                        {
                            var nextCorner = corner.corners[j];
                            var distanceToPlateBoundary = corner.distanceToPlateBoundary + border.Length();
                            if (nextCorner.distanceToPlateBoundary == 0.0f || nextCorner.distanceToPlateBoundary > distanceToPlateBoundary)
                            {
                                elevationBorderQueue.Add(new ElevationBorder(front.origin, border, corner, nextCorner, distanceToPlateBoundary));
                            }
                        }
                    }
                }
            }
            elevationBorderQueue.RemoveRange(0, iEnd);
            elevationBorderQueue.Sort(delegate(ElevationBorder A, ElevationBorder B) {
                return A.distanceToPlateBoundary.CompareTo(B.distanceToPlateBoundary);
            });
        }
    }
    void CalculateTileAverageElevations(List<Tile> tiles)
    {
        for (var i = 0; i < tiles.Count; ++i)
        {
            var tile = tiles[i];
            var elevation = 0f;
            for (var j = 0; j < tile.corners.Length; ++j)
            {
                elevation += tile.corners[j].elevation;
            }
            tile.elevation = elevation / tile.corners.Length;
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
        if (surface.topology.partition.IntersectRay(ray, out tile, out partition, this.transform))
        {
            return true;
        }
        return false;
    }
    #endregion
}

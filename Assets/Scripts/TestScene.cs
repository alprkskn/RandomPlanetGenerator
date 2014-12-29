using UnityEngine;
using System.Collections;
using System.Collections.Generic;

public class TestScene : MonoBehaviour {

    public int planetCount;
    public float minScale;
    public float maxScale;
    public Vector3 extents;

    private List<GameObject> planets;

	// Use this for initialization
	void Start () {
        this.planets = new List<GameObject>();
        for(int i = 0; i <planetCount; i++)
        {
            GameObject go = new GameObject("Planet " + i);
            planets.Add(go);
            Planet p = go.AddComponent<Planet>();
            p.subdivisionLevel = Random.Range(1, 40);
            go.transform.position = new Vector3(Random.Range(-extents.x, extents.x), Random.Range(-extents.y, extents.y), Random.Range(-extents.z, -extents.z));
            float scale = Random.Range(minScale, maxScale);
            go.transform.localScale = Vector3.one * scale;
            p.Generate();
        }

        
        SmoothFollow sf = Camera.main.GetComponent<SmoothFollow>();
        int index = Random.Range(0, planets.Count - 1);
        sf.target = planets[index].transform;
        sf.distance = 2 * sf.target.localScale.x;
        sf.height = 0;
	}
	
	// Update is called once per frame
	void Update () {

        SmoothFollow sf = Camera.main.GetComponent<SmoothFollow>();
        MouseOrbit mo = Camera.main.GetComponent<MouseOrbit>();

        if(sf.enabled)
        {
            if(Vector3.Distance(sf.target.position, Camera.main.transform.position) - sf.distance < sf.distance / 100f)
            {
                sf.enabled = false;
                mo.enabled = true;
                mo.target = sf.target;
                mo.distance = sf.distance;
            }
        }

        if(Input.GetKeyDown(KeyCode.Space))
        {
            mo.enabled = false;
            int index = Random.Range(0, planets.Count - 1);
            sf.target = planets[index].transform;
            sf.distance = 2 * sf.target.localScale.x;
            sf.height = 0;
            sf.enabled = true;
        }
	}
}

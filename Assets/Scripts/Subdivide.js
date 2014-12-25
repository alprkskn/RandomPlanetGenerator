//var material : Material;
var buildmesh : boolean = false;
private var verts : Vector3[];
private var norms : Vector3[];
private var uvs : Vector2[];
private var trigs : int[];
private var mesh : Mesh;
private var originalMesh : Mesh;
//var subdivision : boolean = false;
//var useobjectsmaterial : boolean = false;
//var forceaddmeshcollider : boolean = false;
//var sides : boolean = false;
//var middle : boolean = false;
 
function Update () {
	if(Input.GetKeyDown("1"))
		subdivide(false);
	if(Input.GetKeyDown("2"))
		subdivide(true);
	if(Input.GetKeyDown("x")){
		CopyMesh(originalMesh, mesh);
	}
}
 
function Awake() {
 
	///////////////////////////////////////////////////////////
	//// initialize
	///////////////////////////////////////////////////////////
	//var buildmesh=false;
	if(GetComponent(MeshFilter) == null){
		buildmesh=true;
		gameObject.AddComponent(MeshFilter);
	}
	if(GetComponent(MeshRenderer) == null)
		gameObject.AddComponent(MeshRenderer);
	//if(forceaddmeshcollider) {
	//	if(GetComponent(MeshCollider) == null)
	//		gameObject.AddComponent(MeshCollider);
	//}
 
	updatemesh();
 
	//if(!useobjectsmaterial)
	//	GetComponent(MeshRenderer).material = material;
 
	///////////////////////////////////////////////////////////
	//// build mesh
	///////////////////////////////////////////////////////////
 
	if(buildmesh) {
 
		//// vertices
 
		verts = [Vector3(0,-1,0),Vector3(1,1,0),Vector3(-1,1,0)
		,Vector3(0,1,-1)];
 
		//// uvs
 
		uvs = [Vector2(0,0),Vector2(0,1),Vector2(1,0),Vector2(0,0)];
 
		//// triangles
 
		trigs = [0,1,2,1,3,2];
 
		applymesh();
		mesh.RecalculateNormals();
 
		Debug.Log(verts.length);
	}
 
	originalMesh=new Mesh();
	CopyMesh(mesh, originalMesh);
}
 
function subdivide(center : boolean) {
 
	verts = mesh.vertices;
	trigs = mesh.triangles;
	uvs = mesh.uv;
	norms = mesh.normals;
 
	Debug.Log("enter subdividing: "+verts.length);
 
	var nv = new Array(verts);
	var nt = new Array(trigs);
	var nu = new Array(uvs);
	var nn = new Array(norms);
 
	if(!center) {
 
	for(var i : int = 2;i<trigs.length;i+=3) {
 
			var p0trigwho : int = trigs[i-2];
			var p1trigwho : int = trigs[i-1];
			var p2trigwho : int = trigs[i];
 
			var p0trigwhere : int = i-2;
			var p1trigwhere : int = i-1;
			var p2trigwhere : int = i;
 
			var p0 : Vector3 = verts[p0trigwho];
			var p1 : Vector3 = verts[p1trigwho];
			var p2 : Vector3 = verts[p2trigwho];
 
			var pn0 : Vector3 = norms[p0trigwho];
			var pn1 : Vector3 = norms[p1trigwho];
			var pn2 : Vector3 = norms[p2trigwho];
 
			var pu0 : Vector2 = uvs[p0trigwho];
			var pu1 : Vector2 = uvs[p1trigwho];
			var pu2 : Vector2 = uvs[p2trigwho];
 
			var p0mod : Vector3 = (p0+p1)/2;	
			var p1mod : Vector3 = (p1+p2)/2;
			var p2mod : Vector3 = (p0+p2)/2;
 
			var pn0mod : Vector3 = ((pn0+pn1)/2).normalized;	
			var pn1mod : Vector3 = ((pn1+pn2)/2).normalized;
			var pn2mod : Vector3 = ((pn0+pn2)/2).normalized;
 
			var pu0mod : Vector2 = (pu0+pu1)/2;	
			var pu1mod : Vector2 = (pu1+pu2)/2;
			var pu2mod : Vector2 = (pu0+pu2)/2;
 
			var p0modtrigwho = nv.length;
			var p1modtrigwho = nv.length+1;
			var p2modtrigwho = nv.length+2;
 
			nv.push(p0mod);
			nv.push(p1mod);
			nv.push(p2mod);
 
			nn.push(pn0mod);
			nn.push(pn1mod);
			nn.push(pn2mod);
 
			nu.push(pu0mod);
			nu.push(pu1mod);
			nu.push(pu2mod);
 
			nt[p0trigwhere] = p0trigwho;
			nt[p1trigwhere] = p0modtrigwho;
			nt[p2trigwhere] = p2modtrigwho;
 
			nt.push(p0modtrigwho);
			nt.push(p1modtrigwho);
			nt.push(p2modtrigwho);
 
			nt.push(p0modtrigwho);
			nt.push(p1trigwho);
			nt.push(p1modtrigwho);
 
			nt.push(p2modtrigwho);
			nt.push(p1modtrigwho);
			nt.push(p2trigwho);
		}
 
	}else{
 
		for(var ii : int = 2;ii<trigs.length;ii+=3) {
 
			var p0trigwhomi : int = trigs[ii-2];
			var p1trigwhomi : int = trigs[ii-1];
			var p2trigwhomi : int = trigs[ii];
 
			var p0trigwheremi : int = ii-2;
			var p1trigwheremi : int = ii-1;
			var p2trigwheremi : int = ii;
 
			var p0mi : Vector3 = verts[p0trigwhomi];
			var p1mi : Vector3 = verts[p1trigwhomi];
			var p2mi : Vector3 = verts[p2trigwhomi];
 
			var p0mn : Vector3 = norms[p0trigwhomi];
			var p1mn : Vector3 = norms[p1trigwhomi];
			var p2mn : Vector3 = norms[p2trigwhomi];
 
			var p0mu : Vector2 = uvs[p0trigwhomi];
			var p1mu : Vector2 = uvs[p1trigwhomi];
			var p2mu : Vector2 = uvs[p2trigwhomi];
 
			var p0modmi : Vector3 = (p0mi+p1mi+p2mi)/3;
			var p0modmn : Vector3 = ((p0mn+p1mn+p2mn)/3).normalized;
			var p0modmu : Vector2 = (p0mu+p1mu+p2mu)/3;	
 
			var p0modtrigwhomi = nv.length;
 
			nv.push(p0modmi);
			nn.push(p0modmn);
			nu.push(p0modmu);
 
			nt[p0trigwheremi] = p0trigwhomi;
			nt[p1trigwheremi] = p1trigwhomi;
			nt[p2trigwheremi] = p0modtrigwhomi;
 
			nt.push(p0modtrigwhomi);
			nt.push(p1trigwhomi);
			nt.push(p2trigwhomi);
 
			nt.push(p0trigwhomi);
			nt.push(p0modtrigwhomi);
			nt.push(p2trigwhomi);
 
		}
 
	}
 
	verts = nv.ToBuiltin(Vector3);
	norms = nn.ToBuiltin(Vector3);
	uvs = nu.ToBuiltin(Vector2);
	trigs = nt.ToBuiltin(int);
 
	//applyuvs();
	applymesh();
	//mesh.RecalculateNormals();
 
	Debug.Log("exit subdividing: "+verts.length);
}
 
function applyuvs() {
	uvs = new Vector2[verts.length];
	for(var i : int = 0;i<verts.length;i++)
		uvs[i] = Vector2(verts[i].x,verts[i].y);
}
 
function updatemesh() {
	//mesh = new Mesh();
	mesh = GetComponent(MeshFilter).mesh;
}
 
function applymesh() {
	print(verts.length);
	if(verts.length > 65000){
		Debug.Log("Exiting... Too many vertices");
		return;
	}
	mesh.Clear();
	mesh.vertices = verts;
	mesh.uv = uvs;
	if(mesh.uv2 != null)
		mesh.uv2 = uvs;
	mesh.normals = norms;
	mesh.triangles = trigs;
	mesh.RecalculateBounds();
	if(GetComponent(MeshCollider) != null)
		GetComponent(MeshCollider).sharedMesh = mesh;
	updatemesh();
}
 
function CopyMesh(fromMesh : Mesh, toMesh : Mesh){
	toMesh.Clear();
	toMesh.vertices=fromMesh.vertices;
	toMesh.normals=fromMesh.normals;
	toMesh.uv=fromMesh.uv;
	toMesh.triangles=fromMesh.triangles;
}
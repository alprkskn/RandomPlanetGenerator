var target : Transform;
var distance = 10.0;
var cameraSpeed = 5;

var xSpeed = 175.0;
var ySpeed = 75.0;

var yMinLimit = 20; //Lowest vertical angle in respect with the target.
var yMaxLimit = 80;

var minDistance = 5; //Min distance of the camera from the target
var maxDistance = 20;

private var x = 0.0;
private var y = 0.0;

function Start () {
    var angles = transform.eulerAngles;
    x = angles.y;
    y = angles.x;

    // Make the rigid body not change rotation
    if (rigidbody)
        rigidbody.freezeRotation = true;
}

function LateUpdate () {
    if (target && camera) {

        //Zooming with mouse
        distance += Input.GetAxis("Mouse ScrollWheel")*distance;
        distance = Mathf.Clamp(distance, minDistance, maxDistance);

        //Detect mouse drag;
        if(Input.GetMouseButton(0))   {
   
            x += Input.GetAxis("Mouse X") * xSpeed * 0.02;
            y -= Input.GetAxis("Mouse Y") * ySpeed * 0.02;       
        }
        y = ClampAngle(y, yMinLimit, yMaxLimit);
             
        var rotation = Quaternion.Euler(y, x, 0);
        var position = rotation * Vector3(0.0, 0.0, -distance) + target.position;
         
        transform.position = Vector3.Lerp (transform.position, position, cameraSpeed*Time.deltaTime);
        transform.rotation = rotation;     
    }
}

static function ClampAngle (angle : float, min : float, max : float) {
    if (angle < -360)
        angle += 360;
    if (angle > 360)
        angle -= 360;
    return Mathf.Clamp (angle, min, max);
} 
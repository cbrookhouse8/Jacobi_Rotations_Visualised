// https://stemkoski.github.io/Three.js/Shapes.html

var scene;
var camera;
var renderer;
var light;
var objects;
var dim;
var sep;
var pyrHeight;
var orbitRadius = 500;
var orbitTheta = 0;
var Jacobi; 
var maxLogVal;
var interval = 0;
var currentMatrix;

function init() {
	scene = new THREE.Scene();
	camera = new THREE.PerspectiveCamera(75, 1, 1, 10000);
    camera.position.set(orbitRadius*Math.sin(orbitTheta),
                        0,
                        orbitRadius*Math.cos(orbitTheta));

	renderer = new THREE.WebGLRenderer( {antialias: true} );
	var canvas = renderer.domElement;
	document.body.appendChild(canvas);

	objects = [];
    dim = 10;
    sep = 30;
    pyrHeight = 4*sep;

    Jacobi = new EigenDecomposer(dim);

    maxLogVal = Math.log(Jacobi.getMaxAbsElement()+1);
    
    currentMatrix = Jacobi.getM();

	for (var i = 0; i < dim*dim ; i++) {
		addObject();
	}

	light = addLight();

	window.addEventListener("resize", window_onResize);
	window_onResize();
	renderLoop();
}

function addObject() {
	var geometry;
    
    // Using wireframe materials to illustrate shape details.

	var darkMaterial = new THREE.MeshBasicMaterial({ color: /* 0xffffcc */ 0x0099FF});
	var wireframeMaterial = new THREE.MeshBasicMaterial( { color: 0x000000, wireframe: true, transparent: true } ); 
	var multiMaterial = [darkMaterial, wireframeMaterial]; 

   
      var rad = sep / Math.pow(2,0.5);      
          geometry = new THREE.CylinderGeometry(0,   // radius at the top
                                                rad, // radius at the bottom
                                                pyrHeight,  // height
                                                4,   // segments around the radius
                                                4);  // segments along height

    var pyramid = THREE.SceneUtils.createMultiMaterialObject(geometry,multiMaterial);
	    //pyramid.position.set(0, 50, -100);

    /*
	var material = new THREE.MeshPhongMaterial({
		color: 0xff0000
	});
    */

//	var wireframeMaterial = new THREE.MeshBasicMaterial( { color: 0x000000, wireframe: true, transparent: true } ); 

    // remember for hex colors beginning with # to replace '#' with '0x'
    /*
    var material = new THREE.MeshBasicMaterial({
                                       color: 0x0066ff,
                                       wireframe:true
                                      });
    */
    
//	var mesh = new THREE.Mesh(geometry,wireframeMaterial);
        
       //             resetConfig(mesh);
    var mesh = pyramid;
                    scene.add(mesh);
                    objects.push(mesh);
        mesh.rotation.y = Math.PI/4;
    	return mesh;
}

function resetConfig(mesh) {
	mesh.position.y = -100;
	mesh.rotationSpeed = 0.001 + Math.random() * Math.PI / 60;
	mesh.orbitT = Math.random() * 2 * Math.PI;
	mesh.orbitRadius = 10 + Math.random() * 300;
	mesh.orbitSpeed = 0.01 + Math.random() * Math.PI / 100;
	mesh.verticalSpeed = 1 + Math.random();

	mesh.material.color.setRGB(Math.random(), Math.random(), Math.random());
	mesh.blendSpeed = 0.1;
	mesh.blendSign = 1;

	var r = Math.random();
	if (r < 1 / 3) mesh.blendChannel = "r";
	else if (r < 2 / 3) mesh.blendChannel = "g";
	else mesh.blendChannel = "b";

}

function addLight() {
	var light = new THREE.PointLight(0xffffff, 1, 0);
	light.position.set(50, 0, 0);
	scene.add(light);
	return light;
}

function window_onResize() {
	renderer.setSize(window.innerWidth, window.innerHeight);
	camera.aspect = window.innerWidth / window.innerHeight;
	camera.updateProjectionMatrix();
	renderer.render(scene, camera);
}

function renderLoop() {
    var n = dim;
    if (interval === 0) {
       var stat = Jacobi.nextElement();
       currentMatrix = Jacobi.getM();
    } 

	for (var i = 0; i < n*n; i++) {
        var x_id = i % n;
        var z_id = ~~(i / n);
		var object = objects[i];

        object.position.y = -200;
        object.position.x = -sep*(dim-1)/2 + x_id*sep;
        object.position.z = -sep*(dim-1)/2 + z_id*sep;
        
        var absolute_element = Math.log(Math.abs(currentMatrix[z_id*n+x_id])+1);
        if (maxLogVal < absolute_element)  {
            maxLogVal = absolute_element;
            i = 0; 
        }
        var sf = absolute_element/maxLogVal; 
        // if scale.y == 0, then threeJS throws errors!
        object.scale.y = constrain(sf,0.00001,1); 
        object.position.y += (pyrHeight*sf)/2; 

        // see Junior's code for example use of blenders
	}

    orbitTheta += 0.005;
    camera.position.x = orbitRadius*Math.sin(orbitTheta);
    camera.position.z = orbitRadius*Math.cos(orbitTheta);                   0,
    camera.rotation.y = orbitTheta;

	renderer.render(scene,camera);
	requestAnimationFrame(renderLoop);
    interval = interval < 20 ? interval + 1 : 0;
}

function constrain(x,lower,upper) {
    if (x < lower) {
        return lower;
    } else if (x > upper) {
        return upper;
    } else {
        return x;
    }
}

init();

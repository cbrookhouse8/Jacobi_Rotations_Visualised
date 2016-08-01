// looks really grainy 
// https://stemkoski.github.io/Three.js/Shapes.html

var scene;
var camera;
var renderer;
var light;
var objects;
var dim;
var sep;

function init() {
	scene = new THREE.Scene();
	camera = new THREE.PerspectiveCamera(75, 1, 1, 10000);

	renderer = new THREE.WebGLRenderer();
	var canvas = renderer.domElement;
	document.body.appendChild(canvas);

	objects = [];
    dim = 10;
    sep = 30;
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

	var darkMaterial = new THREE.MeshBasicMaterial({ color: 0xffffcc });
	var wireframeMaterial = new THREE.MeshBasicMaterial( { color: 0x000000, wireframe: true, transparent: true } ); 
	var multiMaterial = [darkMaterial, wireframeMaterial]; 

   
      var rad = sep / Math.pow(2,0.5); 
          geometry = new THREE.CylinderGeometry(0,   // radius at the top
                                                rad, // radius at the bottom
                                                sep*4,  // height
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
	for (var i = 0; i < objects.length; i++) {
        var x_id = i % dim;
        var z_id = ~~(i / dim);
		var object = objects[i];

		object.rotation.y = Math.PI/4;

        object.position.y = -100;
        object.position.x = -sep*(dim-1)/2 + x_id*sep;
        object.position.z = -400 - sep*(dim-1)/2 + z_id*sep;

        object.scale.y = 1 + 0.1*Math.pow(dim-1-z_id,1.2); 
        
        // see Junior's code for example use of blenders
	}

	renderer.render(scene,camera);
	requestAnimationFrame(renderLoop);
}

init();

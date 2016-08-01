var sqWidth = 30;

var n = 20;           // dimension of symmetric matrix
var M;         // symmetric matrix to be diagonalized
var V;         // becomes the matrix of eigenvectors (i.e. R_i * ... * R_1) 
var D;           // stores the diagonal entries (eigenvalues of M)

var sweep;           // used to count how many sweeps have been done
var q,p;             // used to set and read M(p,q) 
var tresh;         // the threshold
var theta;
var tau;           // convenience variable: sin(theta) / (1 + cos(theta)   

var s;             // sin(theta)
var c;             // cos(theta)
var t;             // tan(theta)

var h;             // convenience variable: Use h as tM(p,q) or as M(q,q)-M(p,p)  

// helper arrays

var b;           // helps with calculation of D (eigenvalues)

var z;           // accumulates terms of the form tan(theta)*a[p][q] (eq. 11.1.14)
                 // helps with calculation of D

var maxSweeps;       // maximum number of sweeps we'll allow

function setup() {
	var container = document.getElementById("myContainer");
	var myCanvas = createCanvas(600,600);
			myCanvas.parent(container);
		
	    M = new Array(n*n);
        V = new Array(n*n);
        D = new Array(n); // to be initialized to the diagonal of M[][]

        b = new Array(n); // to be initialized to the diagonal of M[][]
        z = new Array(n); // will be initialized to the zero vector 

       
       // generate an arbitrary symmetric matrix of dimension n
       for (var i = 0 ; i < n ; i++) {
         for (var j = i ; j < n ; j++) {
           var r1 = random(0,1);
               r1 = r1 < 0.5 ? r1+0.3 : r1;
           var val = Math.round(20*r1)+0.0;
             M[j*n+i] = val;
             if (i != j) M[i*n+j] = val;
         }
       }
       
				// {{19,15,16},{15,16,12},{16,12,12}}
       /* 
        M[0] = 19;
        M[1] = 15;
        M[2] = 16;
        M[3] = 15;
        M[4] = 16;
        M[5] = 12;
        M[6] = 16;
        M[7] = 12;
        M[8] = 12;
      */ 
       
        console.log("");
        console.log("Initial matrix ");
        console.log("");
        printMat(M,n,"M");
        
        // Initialize the eigenvector matrix to the identity matrix
        for (var i = 0 ; i < n ; i++) {
            for (var j = 0 ; j < n ; j++) {
               V[j*n+i] = 0.0;
            }
            V[i*n+i] = 1.0;
        }
       
        console.log(""); 
        console.log("Initial V ");
        console.log("");
        printMat(V,n,"V");
        
        // Initialize b and d to the diagonal of M[][]
        for (var i = 0 ; i < n ; i++) {
           b[i] = D[i] = M[i*n+i]+0;
           
           // accumulates terms of the form tan(theta)*a[p][q] (eq. 11.1.14)
           z[i] = 0.0;
        }

            maxSweeps = 6;
            sweep = 1;
            p = 0;
            q = p+1;
   
		
}

function draw() {
  background('white');
    fill(100);
  	textAlign(LEFT,CENTER);
  
  if (frameCount%1 == 0) {
     if (sweep != maxSweeps) {
        nextElement();
     }
  }

  drawSymmetricMatrix(  width/2,      // x
                        height/2,     // y
                        sqWidth,                     // square width
                        5,                           // padding
                        width,                       // screen width
                        height,                      // screen height
                        20,                          // maximum value  in the matrix
                        M,D,n);
    	
}

function nextElement() {        
        var epsilon = 0.00004539992; 
        
            // Check if M has converged to a diagonal matrix
             
            // Check the normal return: relies on quadratic convergence to machine underflow.
    
            var sum = 0.0;
    
            // sum (upper) off-diagonal elements...
            for (var i = 0 ; i < n-1 ; i++) {
                for (var j = i+1 ; j < n ; j++) {
                    sum += Math.abs(M[j*n+i]);
                }
            } 
            
            // Maybe move this to the bottom so that
            // one extra sweep isn't done
            if(sum === 0.0){
                console.log("");
                console.log("Decomposition Complete in "+sweep+" sweeps");
                sweep = maxSweeps;

                // sort the eigenvectors by the magnitude of the Eigenvalues

                // each index in D corresponds to a column in V 
                
                // bubble sort 
                var sorted = false;
                var indexMaps = new Array(n);
                for (var j = 0 ; j < indexMaps.length ; j++) indexMaps[j] = j; 
                
                    while (!sorted) {
                        sorted = true;
                        for (var ct = 1 ; ct < D.length ; ct++) {
                              if (Math.abs(D[ct]) > Math.abs(D[ct-1])) {
                                sorted = false;
                                var temp = D[ct];
                                D[ct] = D[ct-1];
                                D[ct-1] = temp;

                                var tempInt = indexMaps[ct];
                                indexMaps[ct] = indexMaps[ct-1];
                                indexMaps[ct-1] = tempInt; 
                              }
                        }
                    }
                    
                    console.log("");

                    for (var i = 0 ; i < D.length ; i++) {
                        console.log("Lambda (check) "+i+" : "+D[i]);
                    }

                    console.log("");

                    var newV = new Array(n*n);

                    // now use indexMaps to swap the vectors
                    for (var j = 0 ; j < n ; j++) {
                        for (var i = 0 ; i < n ; i++) {
                            // indexMaps[j] says that at column J of V' there should be V[indexMaps[j]]
                            newV[j*n+i] = V[indexMaps[j]*n+i];
                        }
                    }
                   
                    for (var i = 0 ; i < n*n ; i++) {
                        var row = Math.floor(i/n);
                        var col = (i+n)%n;    
                        V[col*n+row] = newV[col*n+row];
                    }
                
                console.log("");
                for (var i = 0 ; i < D.length ; i++) {
                    console.log("Lambda "+i+" : "+D[i]);
                }
                console.log("");
                printMat(V,n,"V");
            }

            // THE ALGORITHM PROPER
           
            // for now, this block isn't of use- but perhaps use tresh later 
            if (sweep < 4) {
                tresh = 0.2*sum/(n*n);
            } else {
                tresh = 0.0;
            }
    
            // iterate over M[][] above the main diagonal
            // if p < n-1...
            //   if q < n...
                  console.log(""); 
                  console.log("In sweep "+sweep+" eliminate M("+p+","+q+") = "+M[q*n+p]);
                  
                  var pq_sizeCheck = 100.0*Math.abs(M[q*n+p]);

                    if (Math.abs(M[q*n+p]) < epsilon) {
                          M[q*n+p] = 0;
                          console.log("");
                          console.log("M["+p+"]["+q+"] is zeroed");
                        }
                    
                    // After four sweeps, skip the rotation if 
                    // the off-diagonal element is small.
                  
                  /*
                    if (i > 4 && (float)(Math.abs(d[p])+pq_sizeCheck) == (float)Math.abs(d[p])
                      && (float)(Math.abs(d[q])+pq_sizeCheck) == (float)Math.abs(d[q])) { 
                        // skip the rotation
                        // zero M(p,q)
                        M[p][q] = 0.0;
                    } else if (Math.abs(M[p][q]) > tresh) {
                       // do the rotation
                  */
    
                      if (true) {           // remove this line at some point
                                            // and uncomment the above  
    
                        // Use h here to be: h = M(q,q)-M(p,p)
                        h = D[q]-D[p];    // see 11.1.7, 11.1.8 which will give tan(theta)
    
                        // set tan(theta) or its proxy
                        if ((Math.abs(h)+pq_sizeCheck) === Math.abs(h)) {
                             t = M[q*n+p]/h;
                        } else {
                             theta = 0.5*h/M[q*n+p]; // Equation (11.1.10). 
                             t = 1.0/(Math.abs(theta)+Math.sqrt(1.0+theta*theta));
                             if (theta < 0.0) t = -t;
                        }
                    
                        // using tan(theta), find cos(theta), sin(theta) and tau   
                        c = 1.0/(Math.sqrt(1+t*t)); // 11.1.11
                        s = t*c;             // 11.1.12
                        tau = s/(1.0+c);     // 11.1.18
    
                        // Repurpose h here for use as tM(p,q).
                        // This term features in 11.1.14-15
                        h = t*M[q*n+p];
    
                        // z accumulates these terms of form
                        // tan(theta)*M[p][q] (eq. 11.1.14)
                        
                        // remember, z[] is reinitialized to 0,
                        // after every sweep
    
                        z[p] -=  h; 
                        z[q] +=  h; 
    
                        // Set M'(p,p), M'(q,q), M'(p,q)
    
                        // M'(p,p) = M(p,p)-tM(p,q) (11.1.14)
                        // M'(q,q) (11.1.15)
    
                        D[p] -=  h; // == M'(p,p) 
                        D[q] +=  h; // == M'(q,q)
    
                        // Set M'(p,q) = 0; (11.1.13)
                        M[q*n+p] = 0.0;
    
                        // Now Update:
                        // M'(j,p), M'(j,q), M'(p,j), M'(q,j)  
    
                        // remember, p < q from how the loops 
                        // are configured
    
                        // Case of rotations 1 ≤ j < p.
                        // M’(j,p) and M’(j,q) for 1 <= j < p  
                        for (var j = 0; j < p; j++) {
                            ROTATE(M,j,p,j,q,s,tau);
                        }
                      
                        // Case of rotations p < j < q.
                        // M’(p,j) and M’(j,q) for p <= j < q  
                        for (var j = p+1 ; j < q ; j++) {
                            ROTATE(M,p,j,j,q,s,tau);
                        }
                        
                        // Case of rotations q < j ≤ n.
                        // M’(p,j) and M’(q,j) for q < j <= n 
                        for (var j = q+1 ; j < n ; j++) {
                            ROTATE(M,p,j,q,j,s,tau);
                        }
                        
                        // update the rotation "matrix"
                        // we only need to modify the cols, p and q
                        for (var j = 0 ; j < n ; j++) {
                            ROTATE(V,j,p,j,q,s,tau); 
                        }
                        
                    } // if true
                    
                    if (q < n-1) {
                     q++; 
                    } else {
                     p++; 
                     q = p+1;
                    }
                    
                    if (p == n-1) {
                         console.log("");
                         console.log("End of sweep "+sweep);
                         sweep++; // on to next sweep
                         p = 0;
                         q = p+1; 
        
                        // Update d with the sum of tapq, and reinitialize z.
                        for (var ip = 0; ip < n; ip++) {
                            b[ip] += z[ip];
                            D[ip] = b[ip];
                            z[ip] = 0.0;
                        }
                    }
    } // end of nextElement()

function printMat(mat,n,name) {
		for (var i = 0 ; i < n ; i++) {
			for (var j = 0 ; j < n ; j++) {
				console.log(name+"("+i+","+j+") = "+mat[j*n+i]);
			} 
		}
}

function ROTATE(a,i,j,k,l,s,tau) {
                 var temp_1 = a[j*n+i];
                 var temp_2 = a[l*n+k];
                 a[j*n+i] = temp_1-s*(temp_2+temp_1*tau);
                 a[l*n+k] = temp_2+s*(temp_1-temp_2*tau);
}

function drawSymmetricMatrix(x,y,w,pad,wid,hei,maxVal,mat,diag,n) {
      rectMode(CENTER);
      textAlign(CENTER,CENTER);
      var xp;
      var yp;

      for (var i = 0 ; i < n ; i++) {
        for (var j = 0 ; j < n ; j++) {
          
          // we're only interested in showing the upper triangle
          // Nb. we need 0 <= log(val), therefore input must be shifted 
          // by 1 since our input with go down to 0; log(0) = -infinity
    
          // there is a lot more work here than is necessary
          // need to use the symmetry of this particular matrix
          
          if (i < j) { 
           fill(heatMap(Math.log(Math.abs(mat[j*n+i])+1),Math.log(maxVal+1)));
           xp = x-n*w/2+w/2+j*w;
           yp = y-n*w/2+w/2+i*w;
           rect(xp,yp,w,w);
          } else if (i == j) { // show the diagonal entries
           fill(heatMap(Math.log(Math.abs(diag[i])+1),Math.log(maxVal+1)));
           xp = x-n*w/2+w/2+j*w;
           yp = y-n*w/2+w/2+i*w;
           rect(xp,yp,w,w); 
          } else {
           fill(heatMap(Math.log(Math.abs(mat[i*n+j])+1),Math.log(maxVal+1)));
           xp = x-n*w/2+w/2+j*w;
           yp = y-n*w/2+w/2+i*w;
           rect(xp,yp,w,w); 
          }
          
           fill(0);
        }
      }  
}

 // http://www.andrewnoske.com/wiki/Code_-_heatmaps_and_color_gradients
 function heatMap(val,maxVal) {
      var nc = 4;  // i.e. rows, i.e. number of colors
      var value = val > maxVal ? 1 : map(val,0,maxVal,0,1);
    
      // row major format
      var colors = [0,0,255,    // blue
                    0,255,0,    // green
                    255,255,0,  // yellow
                    255,0,0];   // red
                                   
      // A static array of 4 colors:  (blue,   green,  yellow,  red) using {r,g,b} for each.
     
      // desired color will be between these two indexes in "color"
      var idx1;
      var idx2;
      var fractBetween = 0;  // Fraction between "idx1" and "idx2" is where our value is.
     
        value = value * 3;          // map value from 0 - 1 to 0 - 4
        idx1  = Math.floor(value);                  // Desired color will be after this index
        idx2  = idx1 >= 3 ? 3 : idx1+1;     // ... and up to this index.
        fractBetween = value - idx1; // Distance between the two indexes (0-1).
        
        var cols = 3;
    
      // this looks like a linear transformation. With a translation operation
      // check this
      var red   = (colors[idx2*cols+0] - colors[idx1*cols+0])*fractBetween + colors[idx1*cols+0];
      var green = (colors[idx2*cols+1] - colors[idx1*cols+1])*fractBetween + colors[idx1*cols+1];
      var blue  = (colors[idx2*cols+2] - colors[idx1*cols+2])*fractBetween + colors[idx1*cols+2];
      
      return color(red,green,blue);
}

// --- DOM utility functions ---

function appendChildren(dad) {
	for (var i=1 ; i < arguments.length ; i++) {
			dad.appendChild(arguments[i]);
	}
}

function $(s) {
	return document.getElementById(s);
}

function c$(s) {
	return document.createElement(s);
}

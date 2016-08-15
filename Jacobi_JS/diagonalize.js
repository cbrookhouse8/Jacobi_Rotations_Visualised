/*
Though comments may have helped to clarify the Numerical Recipes algorithm
There's still so much work to do to clean this up


Have yet to finish making sure that the init and nextElement functions
are accessing class variables using the 'this' keyword.
*/


function EigenDecomposer(rows) {
        console.log("Eigendecomposition");
        this.n = rows;         // dimension of the matrix
        this.m = rows;         // keep it a square matrix

        this.M;                // symmetric matrix to be diagonalized
        this.V;                // becomes the matrix of eigenvectors (i.e. R_i * ... * R_1) 
        this.D;                // stores the diagonal entries (eigenvalues of M)

        this.sweep;            // used to count how many sweeps have been done
        this.p;
        this.q;                // used to set and read M(p,q) 
        this.tresh;            // the threshold
        this.theta;
        this.tau;              // convenience variable: sin(theta) / (1 + cos(theta)   

        this.s;                // sin(theta)
        this.c;                // cos(theta)
        this.t;                // tan(theta)

        this.h;                // convenience variable: Use h as tM(p,q) or as M(q,q)-M(p,p)  

        // helper arrays

        this.b;           // helps with calculation of D (eigenvalues)

        this.z;           // accumulates terms of the form tan(theta)*a[p][q] (eq. 11.1.14)
                          // helps with calculation of D

        this.maxSweeps;       // maximum number of sweeps we'll allow
        this.maxValue;        // maximum absolute value in the matrix
        this.init();

        this.maxValue;
}

EigenDecomposer.prototype.init = function() {
       var n = this.n;

       this.M = new Array(n*n);
       this.V = new Array(n*n);
       this.D = new Array(n); // to be initialized to the diagonal of M[][]

       this.b = new Array(n); // to be initialized to the diagonal of M[][]
       this.z = new Array(n); // will be initialized to the zero vector 
       
       // generate an arbitrary symmetric matrix of dimension n
       this.M = this.RealSymmetricMatrix(n);
       this.maxValue = this.findMaxAbsoluteValue(this.M);
       console.log(this.maxValue);       
       // Initialize the eigenvector matrix to the identity matrix
       this.V = this.IdentityMatrix(n); 
       
        // Initialize b and d to the diagonal of M[][]
       for (var i = 0 ; i < n ; i++) {
           this.b[i] = this.D[i] = this.M[i*n+i]+0;
           
           // accumulates terms of the form tan(theta)*a[p][q] (eq. 11.1.14)
           this.z[i] = 0.0;
       }

       this.maxSweeps = 6;
       this.sweep = 1;
       this.p = 0;
       this.q = this.p+1;
}

EigenDecomposer.prototype.IdentityMatrix = function(n) {
       var M = new Array(n*n);
       for (var i = 0 ; i < n ; i++) {
            for (var j = 0 ; j < n ; j++) M[j*n+i] = 0.0;
            M[i*n+i] = 1.0;
       }
       return M;
}

EigenDecomposer.prototype.RealSymmetricMatrix = function(n) {
       var M = new Array(n*n);
       for (var i = 0 ; i < n ; i++) {
         for (var j = i ; j < n ; j++) {
           var r1 = Math.random();
               r1 = r1 < 0.5 ? r1+0.3 : r1;
           var val = Math.round(20*r1)+0.0;
               M[j*n+i] = val;
               if (i != j) M[i*n+j] = val;
         }
       }
       return M;
}

// sum upper off diagonal elements of a square matrix M
EigenDecomposer.prototype.sumUpperOffDiagonals = function(M,n) {
        var sum = 0.0;
        for (var i = 0 ; i < n-1 ; i++) {
            for (var j = i+1 ; j < n ; j++) {
                sum += Math.abs(this.M[j*n+i]);
            }
        }
        return sum;
}

EigenDecomposer.prototype.nextElement = function() {
        var epsilon = 0.00004539992; 
        
            // Check if M has converged to a diagonal matrix
             
            // Check the normal return: relies on quadratic convergence to machine underflow.
            var n = this.n; 
            var sum = this.sumUpperOffDiagonals(this.M,n);    

            // move this to the bottom so that
            // one extra sweep isn't done?
            if(sum === 0.0){
                console.log("");
                console.log("Decomposition Complete in "+this.sweep+" sweeps");
                this.sweep = this.maxSweeps;
                this.V = this.sortEigenvectors(this.V,this.D,n);
            }

            // THE ALGORITHM PROPER
           
            // for now, this block isn't of use- but perhaps use tresh later 
            if (this.sweep < 4) {
                this.tresh = 0.2*sum/(n*n);
            } else {
                this.tresh = 0.0;
            }
    
            // iterate over M[][] above the main diagonal
            // if this.p < n-1...
            //   if this.q < n...
                  console.log(""); 
                  console.log("In sweep "+this.sweep+" eliminate M("+this.p+","+this.q+") = "+this.M[this.q*n+this.p]);
                  
                  var pq_sizeCheck = 100.0*Math.abs(this.M[this.q*n+this.p]);

                    if (Math.abs(this.M[this.q*n+this.p]) < epsilon) {
                          this.M[this.q*n+this.p] = 0;
                          console.log("");
                          console.log("M["+this.p+"]["+this.q+"] is zeroed");
                    }
                    
                    // After four sweeps, skip the rotation if 
                    // the off-diagonal element is small.
                  
                  /*
                    if (i > 4 && (float)(Math.abs(d[this.p])+pq_sizeCheck) == (float)Math.abs(d[this.p])
                      && (float)(Math.abs(d[this.q])+pq_sizeCheck) == (float)Math.abs(d[this.q])) { 
                        // skip the rotation
                        // zero M(this.p,this.q)
                        M[this.p][this.q] = 0.0;
                    } else if (Math.abs(M[this.p][this.q]) > tresh) {
                       // do the rotation
                  */
    
                      if (true) {           // remove this line at some point
                                            // and uncomment the above  
    
                        // Use h here to be: h = M(this.q,this.q)-M(this.p,this.p)
                        this.h = this.D[this.q]-this.D[this.p];    // see 11.1.7, 11.1.8 which will give tan(theta)
    
                        // set tan(theta) or its proxy
                        if ((Math.abs(this.h)+pq_sizeCheck) === Math.abs(this.h)) {
                             this.t = this.M[this.q*n+this.p]/this.h;
                        } else {
                             this.theta = 0.5*this.h/this.M[this.q*n+this.p]; // Equation (11.1.10). 
                             this.t = 1.0/(Math.abs(this.theta)+Math.sqrt(1.0+this.theta*this.theta));
                             if (this.theta < 0.0) this.t = -this.t;
                        }
                    
                        // using tan(this.theta), find cos(theta), sin(theta) and this.tau   
                        this.c = 1.0/(Math.sqrt(1+this.t*this.t)); // 11.1.11
                        this.s = this.t*this.c;             // 11.1.12
                        this.tau = this.s/(1.0+this.c);     // 11.1.18
    
                        // Repurpose h here for use as tM(p,q).
                        // This term features in 11.1.14-15
                        this.h = this.t*this.M[this.q*n+this.p];
    
                        // z accumulates these terms of form
                        // tan(theta)*M[p][q] (eq. 11.1.14)
                        
                        // remember, z[] is reinitialized to 0,
                        // after every sweep
    
                        this.z[this.p] -=  this.h; 
                        this.z[this.q] +=  this.h; 
    
                        // Set M'(this.p,this.p), M'(this.q,this.q), M'(this.p,this.q)
    
                        // M'(this.p,this.p) = M(this.p,this.p)-tM(this.p,this.q) (11.1.14)
                        // M'(this.q,this.q) (11.1.15)
    
                        this.D[this.p] -=  this.h; // == M'(this.p,this.p) 
                        this.D[this.q] +=  this.h; // == M'(this.q,this.q)
    
                        // Set M'(this.p,this.q) = 0; (11.1.13)
                        this.M[this.q*n+this.p] = 0.0;
    
                        // Now Update:
                        // M'(j,this.p), M'(j,this.q), M'(this.p,j), M'(this.q,j)  
    
                        // remember, this.p < this.q from how the loops 
                        // are configured
    
                        // Case of rotations 1 ≤ j < this.p.
                        // M’(j,this.p) and M’(j,this.q) for 1 <= j < this.p  
                        for (var j = 0; j < this.p; j++) {
                            this.applyRotations(this.M,j,this.p,j,this.q,this.s,this.tau);
                        }
                      
                        // Case of rotations this.p < j < this.q.
                        // M’(this.p,j) and M’(j,this.q) for this.p <= j < this.q  
                        for (var j = this.p+1 ; j < this.q ; j++) {
                            this.applyRotations(this.M,this.p,j,j,this.q,this.s,this.tau);
                        }
                        
                        // Case of rotations this.q < j ≤ n.
                        // M’(this.p,j) and M’(this.q,j) for this.q < j <= n 
                        for (var j = this.q+1 ; j < n ; j++) {
                            this.applyRotations(this.M,this.p,j,this.q,j,this.s,this.tau);
                        }
                        
                        // update the rotation "matrix"
                        // we only need to modify the cols, this.p and this.q
                        for (var j = 0 ; j < n ; j++) {
                            this.applyRotations(this.V,j,this.p,j,this.q,this.s,this.tau); 
                        }
                        
                    } // if true
                    
                    if (this.q < n-1) {
                     this.q++; 
                    } else {
                     this.p++; 
                     this.q = this.p+1;
                    }
                    
                    if (this.p == n-1) {
                         console.log("");
                         console.log("End of sweep "+this.sweep);
                         this.sweep++; // on to next sweep
                         this.p = 0;
                         this.q = this.p+1; 
        
                        // Update d with the sum of tapq, and reinitialize z.
                        for (var ip = 0; ip < n; ip++) {
                            this.b[ip] += this.z[ip];
                            this.D[ip] = this.b[ip];
                            this.z[ip] = 0.0;
                        }
                    }
} // end of nextElement()

EigenDecomposer.prototype.getMaxAbsElement = function() {
    return this.findMaxAbsoluteValue(this.M);
}

EigenDecomposer.prototype.getM = function() {
     var n = this.n;
     var curM = new Array(n*n);
     for (var i = 0 ; i < n ; i++) {
        for (var j = 0 ; j < n ; j++) {
          if (i < j) { 
           // copy from the upper diagonal
           curM[j*n+i] = this.M[j*n+i];
          } else if (i == j) { 
           // copy from the array of diagonal entries
           curM[j*n+i] = this.D[i];
          } else {
           // using the symmetry of M, 
           // copy from the upper triangle
           curM[j*n+i] = this.M[i*n+j];
          }
        }
     }
     return curM;
}

EigenDecomposer.prototype.findMaxAbsoluteValue = function(M) {
   var curMax = 0; 
   for (var i in M) {
       var temp = Math.abs(M[i]);
       if (temp > curMax) curMax = temp;
   }
   return curMax;
}

EigenDecomposer.prototype.sortEigenvectors = function(V,D,n) {
                // sort the eigenvectors by the magnitude of the Eigenvalues
                // (add sorting of eigenvalues to visualisation?) 

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

                    return newV;
}

EigenDecomposer.prototype.applyRotations = function(a,i,j,k,l,s,tau) {
                 var temp_1 = a[j*this.n+i];
                 var temp_2 = a[l*this.n+k];
                 a[j*this.n+i] = temp_1-s*(temp_2+temp_1*tau);
                 a[l*this.n+k] = temp_2+s*(temp_1-temp_2*tau);
}

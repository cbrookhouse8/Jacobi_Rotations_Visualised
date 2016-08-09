function EigenDecomposer(rows,cols) {
        console.log("Eigendecomposition");
        this.n = rows;
        this.m = cols;

        this.n = 20;           // dimension of symmetric matrix
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

        this.init();
}

EigenDecomposer.prototype.init = function() {
        var n = this.n;

        this.M = new Array(n*n);
        this.V = new Array(n*n);
        this.D = new Array(n); // to be initialized to the diagonal of M[][]

        this.b = new Array(n); // to be initialized to the diagonal of M[][]
        this.z = new Array(n); // will be initialized to the zero vector 
       
       // generate an arbitrary symmetric matrix of dimension n
       for (var i = 0 ; i < n ; i++) {
         for (var j = i ; j < n ; j++) {
           var r1 = random(0,1);
               r1 = r1 < 0.5 ? r1+0.3 : r1;
           var val = Math.round(20*r1)+0.0;
               this.M[j*n+i] = val;
               if (i != j) this.M[i*n+j] = val;
         }
       }

       // Initialize the eigenvector matrix to the identity matrix
       for (var i = 0 ; i < n ; i++) {
            for (var j = 0 ; j < n ; j++) V[j*n+i] = 0.0;
            this.V[i*n+i] = 1.0;
       }
       
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

EigenDecomposer.prototype.nextElement = function() {
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


EigenDecomposer.prototype.ROTATE = function(a,i,j,k,l,s,tau) {
                 var temp_1 = this.a[j*this.n+i];
                 var temp_2 = this.a[l*this.n+k];
                 a[j*this.n+i] = temp_1-s*(temp_2+temp_1*tau);
                 a[l*this.n+k] = temp_2+s*(temp_1-temp_2*tau);
}

// before becoming a repository
// this file was stored as 
// numerical_recipes_4

// now trying to animate this a bit

// add a log message for when a val is close to 0 (by Lengyel's criterion)

// a[][] is our symmetric matrix
// a[][] is transformed such that the superdiagonals are 0, the subdiagonals/diagonal
// are the original entries.
// d[] is the eigenvalues
// v[1..n][1..n] is the matrix of eigenvectors
// nrot is the number of Jacobi rotations required for convergence

// Function computes all eigenvalues and eigenvectors of a real symmetric matrix a[1..n][1..n].
// On output, elements of a above the diagonal are destroyed.
// d[1..n] returns the eigenvalues of a. v[1..n][1..n] is a matrix whose columns contain,
// on output, the normalized eigenvectors of a. nrot returns the number of Jacobi rotations
// that were required.

int n = 60;
float sqWidth = 10;
int refresh = 1;

float[][] a;
float[][] v;

int j;               // used as counter in the final rotation steps
int i;               // used to count how many sweeps have been done
int iq,ip;           // used to set and read from a[][]
float tresh;         // the threshold
float theta;
float tau;           // sin(theta) / (1 + cos(theta)   
float t;             // tan(theta)
float sm;            // used as a sum variable e.g. when summing off diag to check state of a[][]
float s;             // sin(theta)
float c;             // cos(theta)
float h;             // Use h as tM(p,q) or as M(q,q)-M(p,p), as well as a 
                     // value holder in ROTATE()

float[] b;             // helps with calculation of d (eigenvalues)
float[] d;             // stores the diagonal entries (eigenvalues)
float g;             // used in the Rotate function
                     // and also for the condition that tests which 
                     // optimization steps to take

float[] z;           // accumulates terms of the form tan(theta)*a[p][q] (eq. 11.1.14)
                     // helps with calculation of d

int maxSweeps;       // maximum number of sweeps we'll allow
int nrot = 0;        // is the number of Jacobi rotations required for convergence

void setup() {
  size(600,600);
  background(255);
  noStroke();
  
    a = new float[n][n];
    v = new float[n][n];
    b = new float[n]; // initialized to the diagonal of a[][]
    d = new float[n]; // initialized to the diagonal of a[][]
    z = new float[n]; // will be initialized to the zero vector 
  
   // generate an arbitrary symmetric matrix of dimension n
   for (i = 0 ; i < n ; i++) {
     for (j = i ; j < n ; j++) {
       float r1 = random(0,1);
             r1 = r1 < 0.5 ? r1+0.3 : r1;
       float val = round(20*r1)+0.0;
         a[i][j] = val;
         if (i != j) a[j][i] = val;
     }
   }
    
    // Initialize the eigenvector matrix to the identity matrix
    for (ip = 0 ; ip < n ; ip++) {
        for (iq = 0; iq < n; iq++) {
           v[ip][iq] = 0.0;
        }
        v[ip][ip] = 1.0;
    }
    
    // Initialize b and d to the diagonal of a[][]
    for (ip = 0 ; ip < n ; ip++) {
       b[ip] = d[ip] = a[ip][ip];
       
       // accumulates terms of the form tan(theta)*a[p][q] (eq. 11.1.14)
       z[ip] = 0.0;
    }

        maxSweeps = 6;
        i = 1;
        ip = 0;
        iq = ip+1;
} // end of setup

void draw() {
    background(255);
    
    drawMatrix(width/2, // x
               height/2,  // y
               sqWidth,                     // square width
               5,                           // padding
               width,                       // screen width
               height,                      // screen height
               20,                          // maximum value  in the matrix
               a,d,n);
    
    float epsilon = 0.00004539992; 
    
    if (i <= maxSweeps && frameCount%refresh == 0) {
      /*
        // The normal return, which relies on quadratic convergence to machine underflow.
        // ...on the first three sweeps. ...thereafter.
        sm = 0.0; // sum variable

        // sum (upper) off-diagonal elements...
        for (ip = 0; ip < n-1; ip++) {
            for (iq = ip+1; iq < n; iq++) {
                sm += Math.abs(a[ip][iq]);
            }
        } 
        
        if(sm == 0.0){
            free_vector(z,1,n);
            free_vector(b,1,n); 
            return;
        }
        
        if (i < 4) {
            tresh = 0.2*sm/(n*n);
        } else {
            tresh = 0.0;
        }
    */

        // iterate over a[][] above the main diagonal
        //if (ip < n-1) {
        //if(iq < n) {
             
              println("In sweep "+i+" eliminate M("+ip+","+iq+") = "+a[ip][iq]);
              
                g = 100.0*Math.abs(a[ip][iq]);
                
                if (Float.isNaN(a[ip][iq])) {
                   a[ip][iq] = 0;
                   println(true); 
                }
                
                if (Math.abs(a[ip][iq]) < epsilon) {
                  a[ip][iq] = 0;
                  println("a["+ip+"]["+iq+"] is zeroed");
                }
                
                // After four sweeps, skip the rotation if 
                // the off-diagonal element is small.

        //        if (i > 4 && (float)(Math.abs(d[ip])+g) == (float)Math.abs(d[ip])
        //          && (float)(Math.abs(d[iq])+g) == (float)Math.abs(d[iq])) { 
        //            // skip the rotation
        //            a[ip][iq] = 0.0;
        //        } else if (Math.abs(a[ip][iq]) > tresh) {
                    // do the rotation

                  if (true) {           // remove this line at some point
                                        // and uncomment the above  

                    // Use h here to be: h = M(q,q)-M(p,p)
                    h = d[iq]-d[ip];    // see 11.1.7, 11.1.8 which will give tan(theta)

                    // set tan(theta) or its proxy
                    if ((float)(Math.abs(h)+g) == (float)Math.abs(h)) {
                         t = (a[ip][iq])/h;
                    } else {
                         theta = 0.5*h/(a[ip][iq]); // Equation (11.1.10). 
                         t = 1.0/((float)Math.abs(theta)+(float)Math.sqrt(1.0+theta*theta));
                         if (theta < 0.0) t = -t;
                    }
                
                    // using tan(theta), find cos(theta), sin(theta) and tau   
                    c = 1.0/((float)Math.sqrt(1+t*t)); // 11.1.11
                    s = t*c;             // 11.1.12
                    tau = s/(1.0+c);     // 11.1.18

                    // Repurpose h here for use as tM(p,q).
                    // This term features in 11.1.14-15
                    h = t*a[ip][iq];

                    // z accumulates these terms of form
                    // tan(theta)*a[p][q] (eq. 11.1.14)
                    
                    // remember, z[] is reinitialized to 0,
                    // after every sweep

                    z[ip] -=  h; 
                    z[iq] +=  h; 

                    // Set M'(p,p), M'(q,q), M'(p,q)

                    // M'(p,p) = M(p,p)-tM(p,q) (11.1.14)
                    // M'(q,q) (11.1.15)

                    d[ip] -=  h; // == M'(p,p) 
                    d[iq] +=  h; // == M'(q,q)

                    // Set M'(p,q) = 0; (11.1.13)
                    a[ip][iq] = 0.0;

            /* 
                    define ROTATE(a,i,j,k,l) { 
                        g = a[i][j];
                        h = a[k][l];
                        a[i][j] = g-s*(h+g*tau);
                        a[k][l] = h+s*(g-h*tau);
                    } 
            */
                    // Now Update:
                    // M'(j,p), M'(j,q), M'(p,j), M'(q,j)  

                    // h,g and repurposed here

                    // remember, p < q from how the loops 
                    // are configured

                    // Case of rotations 1 ≤ j < p.
                    // M’(j,p) and M’(j,q) for 1 <= j < p  
                    for (j = 0; j < ip; j++) {
                        ROTATE(a,j,ip,j,iq,g,h,s,tau);
                    }
                  
                    // Case of rotations p < j < q.
                    // M’(p,j) and M’(j,q) for p <= j < q  
                    for (j = ip+1; j < iq; j++) {
                        ROTATE(a,ip,j,j,iq,g,h,s,tau);
                    }
                    
                    // Case of rotations q < j ≤ n.
                    // M’(p,j) and M’(q,j) for q < j <= n 
                    for (j = iq+1; j < n; j++) {
                        ROTATE(a,ip,j,iq,j,g,h,s,tau);
                    }
                    
                    // update the rotation "matrix"
                    // we only need to modify the cols, p and q
                    for (j = 0; j < n; j++) {
                        ROTATE(v,j,ip,j,iq,g,h,s,tau); 
                    }
                    
                   // nrot++;
                } // if true
                
                if (iq < n-1) {
                 iq++; 
                } else {
                 ip++; 
                 iq = ip+1;
                }
                
                if (ip == n-1) {
                 i++; // onto next sweep
                 ip = 0;
                 iq = ip+1; 


                    // Update d with the sum of tapq, and reinitialize z.
                    for (int iip = 0; iip < n; iip++) {
                        b[iip] += z[iip];
                        d[iip] = b[iip];
                        z[iip] = 0.0;
                    }
                    
                }
        } // end of sweep condition
}

void drawMatrix(float x, float y, float w, float pad, float wid, float hei, float maxVal, float[][] a, float[] diag, int n) {
//void drawMatrix(float x, float y, float w, int n) {
  rectMode(CENTER);
  textAlign(CENTER,CENTER);
  float xp;
  float yp;
  for (int i = 0 ; i < n ; i++) {
    for (int j = 0 ; j < n ; j++) {
      
      // we're only interested in showing the upper triangle
      // Nb. we need 0 <= log(val), therefore input must be shifted 
      // by 1 since our input with go down to 0; log(0) = -infinity

      // there is a lot more work here than is necessary
      // need to use the symmetry of this particular matrix
      
      if (i < j) { 
       //fill(0,0,255,map(log(Math.abs(a[i][j])+1),0,log(maxVal+1),0,255));
       fill(heatMap(log(Math.abs(a[i][j])+1),log(maxVal+1)));
       xp = x-n*w/2+w/2+j*w;
       yp = y-n*w/2+w/2+i*w;
       rect(xp,yp,w,w);
      } else if (i == j) { // show the diagonal entries
      //fill(0,0,255,map(log(Math.abs(diag[i])+1),0,log(maxVal+1),0,255));
         fill(heatMap(log(Math.abs(diag[i])+1),log(maxVal+1)));
       xp = x-n*w/2+w/2+j*w;
       yp = y-n*w/2+w/2+i*w;
       rect(xp,yp,w,w); 
      } else {
       //fill(0,0,255,map(log(Math.abs(a[j][i])+1),0,log(maxVal+1),0,255));
       fill(heatMap(log(Math.abs(a[j][i])+1),log(maxVal+1)));
       xp = x-n*w/2+w/2+j*w;
       yp = y-n*w/2+w/2+i*w;
       rect(xp,yp,w,w); 
      }
      
       fill(0);
       //text(a[i][j],x-n*w/2+w/2+j*w,y-n*w/2+w/2+i*w);
    }
  }  
}

void ROTATE(float[][] a,int i,int j,int k,int l, /* added */ float g, float h, float s, float tau) { 
             g = a[i][j];
             h = a[k][l];
             a[i][j] = g-s*(h+g*tau);
             a[k][l] = h+s*(g-h*tau);
}

void printMatrix(float[][] mat, int n) {
    for (int i = 0 ; i < n ; i++) {
     for (int j = 0 ; j < n ; j++) {
        println("M("+i+","+j+") = "+mat[i][j]);
     } 
    }
}

// http://www.andrewnoske.com/wiki/Code_-_heatmaps_and_color_gradients
color heatMap(float val, float maxVal) {
  int nc = 4;  // i.e. rows, i.e. number of colors
  float value = val > maxVal ? 1 : map(val,0,maxVal,0,1);

  // row major format
  int[] colors = new int[] {0,0,255,    // blue
                            0,255,0,    // green
                            255,255,0,    // yellow
                            255,0,0};   // red
                                           
  // A static array of 4 colors:  (blue,   green,  yellow,  red) using {r,g,b} for each.
 
  // desired color will be between these two indexes in "color"
  int idx1;
  int idx2;
  float fractBetween = 0;  // Fraction between "idx1" and "idx2" is where our value is.
 
    value = value * 3;          // map value from 0 - 1 to 0 - 4
    idx1  = floor(value);                  // Desired color will be after this index
    idx2  = idx1 >= 3 ? 3 : idx1+1;     // ... and up to this index.
    fractBetween = value - ((float) idx1); // Distance between the two indexes (0-1).
    
    int cols = 3;

  // this looks like a linear transformation. With a translation operation
  // check this
  float red   = (colors[idx2*cols+0] - colors[idx1*cols+0])*fractBetween + colors[idx1*cols+0];
  float green = (colors[idx2*cols+1] - colors[idx1*cols+1])*fractBetween + colors[idx1*cols+1];
  float blue  = (colors[idx2*cols+2] - colors[idx1*cols+2])*fractBetween + colors[idx1*cols+2];
  
  return color(red,green,blue);
}
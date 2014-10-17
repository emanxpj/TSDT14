Representation
===

2D
==

* Hom. rep of point and line (proj. space, proj elements, equivalence)

* Point and line at inf

* Transformations (3x3 matrices)
  -Affine (parallell lines stay parallell)
    -Rigid transformations ( rot + translation)
    -Similarity transformations( Rigid + scaling)
    -Reflections
    -Skewing

  -Homographies. (always reversible and general 3x3)

3D
==

* Hom rep of points and planes.

* Plucker coord for lines.

* Affine transformations and Homographies

Distances
==

* Normalization needed (P and D type)

* Point to Point distance (signed)

* Point to line distance (signed)

Pinhole camera
==

* Definition
  - image plane (virtual image plane)
  - camera center, n
  - projection line
  - principle point and axis
  - Camera matrix
     - Cn = 0
  -Internal/external parameters;
	K = (f y u		(R t)
	       f v	
		 1

  -Equivalent cameras
      	-2 images of a single plane
	   	-points in a one image could be transformed 
		 by a simple homography, due to constricted to a single plane. 
       	- Rotation of the image around the camera center and take the next image.

3D rotations SO(3)
==
  
   - Matrix rep: R € SO(3) => R’R = I and det(R) = ±1
   - Angle rep: (a,n), 
	*intuitive but contains a singularity at R = I.
	in that case n is not well defined. (uniquely). 
	*(-a,-n) is equivalent to (a,n), multiple values
	* a and a + 2*pi*k equivalent. 
	
   - Rodrigues formula, (n,a) => R 
   - Unit quaternions q€H ; |q| = 1; q = (cos(a/2), sin(a/2)*n) (4-dim). 
	p_rot =-q*p_x*q^(-1);  where p_x is a pure quaternion. 
   	- Going from q -> R is ’simple’ (you don’t need to know this by heart).
   	- Double Embedding. For each rotation there are exactly two quaternion
		representations and they only differ in sign.

   - Matrix exponential representations.


Estimation
===

* What is estimation? Fitting a model to data, with the help of an error function.
	The error function is dependent of the model parameters and then try
	to minimize the error function. 
	- Different options for the error function exist.

* Geometric error. 
	The most common error function is based around the idea of minimizing
	  distances.
	-This can lead to nonlinear equations and typically these are hard to
	  minimize. Complex optimizations and good first guess needed. i.e. much more 		  advanced methods needed.

* Algebraic error. ( e = |Az|^2 ) ; A = data, z = model parameters. 
	Much less intuitive than the geometric error but much easier to minimize. 
	 -x*p = 0, y*l=0, L*x = 0, if the point is on the plane, line etc.
	 -Each point of data can give one row in A. With the addition of the noise
	  we can no longer get the error to be zero, but we can now minimize the norm
	  or Az. 

	-InHomogeneous method
	-Homogeneous method
		-SVD
		-SVD-profile
		-Solution space (preferably 1 dimensional, but sometimes larger
				 i.e. non-well defined solution.)
		-Data normalisation, Hartley.
		 	-Why it is done.
			 -By changing the coordinate system before optimizing 
			  the profile gets more distinct.
			-How it is done. 
			 -A scaling and a translation (specific)

Matrix exp 
==

*OPP (SOPP). 
	Used for estimation of rigid transformation.

  

	  

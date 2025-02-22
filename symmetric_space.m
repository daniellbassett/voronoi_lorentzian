//voronoi_lorentzian - Daniel Bassett, 2025
//--------------------------------------------------symmetric_space.m--------------------------------------------------
/*	
	Sets up the Lorentzian cone model of hyperbolic n-space using the user-defined parameters from init.

	Provides functions for working with the cone in the context of the Voronoi algorithm.	
*/

import "init.m" : n, standard_form, q, multithread, num_threads;


//----------Lorentzian cone construction----------
//Ambient inner product space
if standard_form then //use [1, ..., 1] bilinear form
	B := MatrixRing(Rationals(), n+1) ! 1; //Identity matrix
else
	B := DiagonalMatrix(Rational(), n+1, q cat [1]);
end if;

V := VectorSpace(Rationals(), n+1, B); //Equipped with bilinear form B

//Minkowski form
signatureMatrix := DiagonalMatrix(Rationals(), [-1 : i in [1..n]] cat [1]); //signature n,1

function minkowskiForm(v,w)
	return InnerProduct(v, w*signatureMatrix);
end function;

function minkowskiNorm(v) //minkowskiForm(v,v)
	return InnerProduct(v, v*signatureMatrix);
end function;

//Lorentzian cone
function isInteriorPoint(v)
	if v[n+1] lt 0 then
		return false;
	else
		return minkowskiNorm(v) gt 0;
	end if;
end function;

//voronoi_lorentzian - Daniel Bassett, 2025
//--------------------------------------------------symmetric_space.m--------------------------------------------------
/*	
	Sets up the Lorentzian cone model of hyperbolic n-space using the user-defined parameters from init.

	Provides general functions for working with the cone:
InnerProduct computes the Euclidean inner product of two vectors
minkowskiForm computes the signature n,1 inner product of two vectors
minkowskiNorm computes the corresponding quadratic form on a vector
isInteriorPoint checks if a point lies in the interior of the Lorentzian cone

	Provides special functions for working with the cone in the context of the Voronoi algorithm:
boundaryPoints computes integral points on the boundary of the cone, which are vertices of the tessellation
*/

import "init.m" : n, standard_form, q, multithread, num_threads;


//----------Lorentzian cone construction----------
//Ambient inner product space
if standard_form then //use [1, ..., 1] bilinear form
	B := MatrixRing(Rationals(), n+1) ! 1; //Identity matrix
else
	B := DiagonalMatrix(Rationals(), n+1, q cat [1]);
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


//----------Integral boundary points----------
/*We iterate over the first n-1 coordinates, and check if this gives a possible solution for the nth coordinate,
given a fixed last coordinate ('height')*/

if q[n] eq 1 then //can avoid division logic, ~10% faster
	function boundaryPoints(height)
		Vlow := VectorSpace(Rationals(), n, DiagonalMatrix(Rationals(), n, [-q[i] : i in [1..n-1]] cat [1]));
		
		points := [];
		
		v := Vlow ! 0;
		v[n] := height;
		
		heightSquared := height^2;
		lastCoordBound := q[n-1] * v[n-1]^2; //to check when the last coordinate being varied is too big
		
		while lastCoordBound le heightSquared do
			norm := InnerProduct(v,v);
			square, coord := IsSquare(norm); //~10% faster version when q[n] = 1
			
			if square then
				Append(~points, V ! ([v[i] : i in [1..n-1]] cat [coord, height]));
			end if;
			
			v[1] +:= 1;
			for i in [1..n-3] do
				if q[i] * v[i]^2 gt heightSquared then
					v[i+1] +:= 1;
					v[i] := 0; //if this is too big then it just got increased, so the previous one has been set to 0; induction.
				end if;
			end for;
			
			if q[n-2] * v[n-2]^2 gt heightSquared then //increasing last varied coord; remember to update lastCoordBound
				v[n-1] +:= 1;
				v[n-2] := 0;
				
				lastCoordBound := q[n]*v[n-1]^2;
			end if;
		end while;
		
		return points;
	end function;
else
	function boundaryPoints(height)
		Vlow := VectorSpace(Rationals(), n, DiagonalMatrix(Rationals(), n, [-q[i] : i in [1..n-1]] cat [1]));
		
		points := [];
		
		v := Vlow ! 0;
		v[n] := height;
		
		heightSquared := height^2;
		lastCoordBound := q[n-1] * v[n-1]^2; //to check when the last coordinate being varied is too big
		
		while lastCoordBound le heightSquared do
			norm := InnerProduct(v,v);
			if Integers() ! norm mod Integers() ! q[n] eq 0 then
				target := ExactQuotient(norm, q[n]);
				square, coord := IsSquare(target);
				
				if square then
					Append(~points, V ! ([v[i] : i in [1..n-1]] cat [coord, height]));
				end if;
			end if;

			v[1] +:= 1;
			for i in [1..n-3] do
				if q[i] * v[i]^2 gt heightSquared then
					v[i+1] +:= 1;
					v[i] := 0; //if this is too big then it just got increased, so the previous one has been set to 0; induction.
				end if;
			end for;
			
			if q[n-2] * v[n-2]^2 gt heightSquared then //increasing last varied coord; remember to update lastCoordBound
				v[n-1] +:= 1;
				v[n-2] := 0;
				
				lastCoordBound := q[n]*v[n-1]^2;
			end if;
		end while;
		
		return points;
	end function;
end if;

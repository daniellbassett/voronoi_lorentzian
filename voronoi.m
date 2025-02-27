//voronoi_lorentzian - Daniel Bassett, 2025
//--------------------------------------------------voronoi.m--------------------------------------------------
/*
TODO: write description
*/

import "init.m" : n;
import "symmetric_space.m" : boundaryPoints, minkowskiForm, minkowskiNorm, cellNormal;


//----------Minimal vector calculation----------
function minimalVectors(p, boundary_point_database, min, minimal_vectors) //calculates the minimum, and minimal vectors, of interior point p. allows input of some known vectors 
	if #known_vecs eq 0 then
		min := Infinity();
		minimal_vectors := [];
		bound := Infinity();
	else
		bound := 2 * p[n+1] * min / minkowskiNorm(p);
	end if;
	
	p_signs := [Sign(p[i]) : i in [1..n+1]]; //minimal vectors will have entries whose signs match those of p
	for i in [1..n+1] do
		if p_signs[i] eq 0 then
			p_signs[i] := 1; //allows defining w[i] := v[i] * p_signs[i] later without having to check for p_signs[i] = 0 at that point.
		end if;
	end for;
		
	//Iterate over boundary vectors by height
	bound := Infinity(); //initial bound, will be decreased as upper bounds to min are found.
	height := 1;
	while height lt bound do
		if height gt #boundary_point_database then
			Append(~boundary_point_database, boundaryPoints(height));
		end if;
		
		for v in boundary_points_database[height] do
			w := V ! [v[i] * p_signs[i] : i in [1..n+1]];
			
			p_component := minkowskiForm(p, w);
			if p_component lt min then //new minimum found
				min := p_component;
				minimal_vectors := [w];
				
				bound := 2 * p[n+1] * min / minkowskiNorm(p); //new upper bound on the height based on new minimum, derived using cauchy-schwarz
			elseif p_component eq min then //same minimum, new minimal vector
				Append(~minimal_vectors, w);
			end if;
		end for;
	end while;
	
	return min, minimal_vectors, boundary_point_database;
end function;

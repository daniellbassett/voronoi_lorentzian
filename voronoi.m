//voronoi_lorentzian - Daniel Bassett, 2025
//--------------------------------------------------voronoi.m--------------------------------------------------
/*
TODO: write description
*/

import "init.m" : n, standard_form;
import "symmetric_space.m" : V, boundaryPoints, minkowskiForm, minkowskiNorm, signOrbit, cellNormal, isInteriorPoint;


//----------Minimal vector calculation----------
function minimalVectors(p, boundary_point_database, min, minimal_vectors) //calculates the minimum, and minimal vectors, of interior point p. allows input of some known vectors 
	if #minimal_vectors eq 0 then
		min := Infinity();
		height_bound := Infinity();
	else
		height_bound := 2 * p[n+1] * min / minkowskiNorm(p);
	end if;
	
	print minkowskiNorm(p);
	
	p_signs := [Sign(p[i]) : i in [1..n+1]]; //minimal vectors will have entries whose signs match those of p
	p_zero_indices := [];
	for i in [1..n+1] do
		if p_signs[i] eq 0 then
			p_signs[i] := 1; //allows defining w[i] := v[i] * p_signs[i] later
			Append(~p_zero_indices, i); //will be used to swap the signs of minimal vectors where applicable
		end if;
	end for;
		
	//Iterate over boundary vectors by height
	height := 1;
	while height lt height_bound do
		if height gt #boundary_point_database then
			Append(~boundary_point_database, boundaryPoints(height));
		end if;
		
		for v in boundary_point_database[height] do
			w := V ! [v[i] * p_signs[i] : i in [1..n+1]];
			
			p_component := minkowskiForm(p, w);
			print p_component, min;
			if p_component lt min then //new minimum found
				min := p_component;
				
				sign_flip_indices := []; //both signs possible when p[i] = 0 and v[i] =/= 0
				for i in p_zero_indices do
					if v[i] ne 0 then
						Append(~sign_flip_indices, i);
					end if;
				end for;
				
				minimal_vectors := signOrbit(w, sign_flip_indices);
				
				print min;
				height_bound := 2 * p[n+1] * min / minkowskiNorm(p); //new upper bound on the height based on new minimum, derived using cauchy-schwarz
			elif p_component eq min then //same minimum, new minimal vector
				sign_flip_indices := []; //both signs possible when p[i] = 0 and v[i] =/= 0
				for i in p_zero_indices do
					if v[i] ne 0 then
						Append(~sign_flip_indices, i);
					end if;
				end for;
				
				for vec in signOrbit(w, sign_flip_indices) do
					Append(~minimal_vectors, vec);
				end for;
			end if;
		end for;
		
		height +:= 1;
	end while;
	
	return min, minimal_vectors, boundary_point_database;
end function;

//----------Perfect vector calculations---------
if standard_form then //[0,...,0,1] is a perfect vector
	function buildInitialPoint(boundary_point_database)
		return V ! ([0 : i in [1..n]] cat [1]), boundary_point_database;
	end function;
else
	function buildInitialPoint(boundary_point_database)
		p := V ! ([0 : i in [1..n]] cat [1]); //Start with some interior point
		
		min, min_vecs, boundary_point_database := minimalVectors(p, boundary_point_database, Infinity(), []);
		p := p / min; //rescale to have minimum 1
		
		rank := Rank(VerticalJoin(min_vecs));
		if rank eq n+1 then //perfect vector found
			return p, min_vecs, boundary_point_database;
		end if;
		
		height := 1;
		min_rho_plus := Infinity();
		min_rho_minus := Infinity();
		while true do //try to add vectors to increase rank, until reach a perfect vector
			normal := cellNormal(min_vecs);
			
			if height gt #boundary_point_database then
				Append(~boundary_point_database, boundaryPoints(height));
			end if;
			
			for v in boundary_point_database[height] do //try to add vectors from signOrbit(v) as minimal vectors
				for w in signOrbit(v, [1..n]) do //may be greater restriction possible, but this is not at all a critical point of the code for performance
					normal_component := minkowskiForm(w, normal);
					
					if normal_component ne 0 then
						p_component := minkowskiForm(w, p);
						rho := (1-p_component) / normal_component;
						
						if (rho gt 0 and rho lt min_rho_plus) or (rho lt 0 and rho gt min_rho_minus) then //if minvec, then we move a minimal amount in the normal (signed) direction
							candidate := p + normal * (1 - p_component)/normal_component; //vector with minkowskiForm 1 with min_vecs and w. if minimum 1 then rank improvement
							
							if isInteriorPoint(candidate) then
								candidate_min, candidate_min_vecs, boundary_point_database := minimalVectors(candidate, boundary_point_database, 1, min_vecs cat [w]);
								
								if candidate_min eq 1 then
									min_vecs := candidate_min_vecs;
									p := candidate;
									
									if Rank(VerticalJoin(min_vecs)) eq n+1 then
										return p, min_vecs, boundary_point_database;
									end if;
									
									min_rho_plus := Infinity(); //travelling in a new  normal direction
									min_rho_minus := Infinity();
								end if;
							end if;
						end if;
					end if;
				end for;
			end for;
			
			height +:= 1; //try again at next height
		end while;
	end function;
end if;

buildInitialPoint([]);

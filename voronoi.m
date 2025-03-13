//voronoi_lorentzian - Daniel Bassett, 2025
//--------------------------------------------------voronoi.m--------------------------------------------------
/*
TODO: write description
*/

import "init.m" : n, standard_form;
import "symmetric_space.m" : V, boundaryPoints, minkowskiForm, minkowskiNorm, signOrbit, cellNormal, isInteriorPoint, facets, equivalent, J;


//----------Minimal vector calculation----------
function minimalVectors(p, boundary_point_database, min, minimal_vectors) //calculates the minimum, and minimal vectors, of interior point p. allows input of some known vectors 
	if #minimal_vectors eq 0 then
		min := Infinity();
		height_bound := Infinity();
	else
		height_bound := 2 * p[n+1] * min / minkowskiNorm(p);
	end if;
	
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
			if p_component lt min then //new minimum found
				min := p_component;
				
				sign_flip_indices := []; //both signs possible when p[i] = 0 and v[i] =/= 0
				for i in p_zero_indices do
					if v[i] ne 0 then
						Append(~sign_flip_indices, i);
					end if;
				end for;
				
				minimal_vectors := signOrbit(w, sign_flip_indices);
				
				height_bound := 2 * p[n+1] * min / minkowskiNorm(p); //new upper bound on the height based on new minimum, derived using cauchy-schwarz
			elif p_component eq min then //same minimum, new minimal vector
				sign_flip_indices := []; //both signs possible when p[i] = 0 and v[i] =/= 0
				for i in p_zero_indices do
					if v[i] ne 0 then
						Append(~sign_flip_indices, i);
					end if;
				end for;
				
				for vec in signOrbit(w, sign_flip_indices) do
					if not vec in minimal_vectors then
						Append(~minimal_vectors, vec);
					end if;
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
		p := V ! ([0 : i in [1..n]] cat [1]);
		_, minimal_vectors := minimalVectors(p, [], Infinity(), []);
		return p, minimal_vectors, boundary_point_database;
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

function neighbour(p, min_vecs, facet_vecs, boundary_point_database)
	normal := cellNormal(facet_vecs);
	
	for v in min_vecs do
		if minkowskiForm(v, normal) lt 0 then
			normal *:= -1; //choose sign so that all minvecs have minkowskiForm(v,n) >= 0
			break;
		end if;
	end for;
	
	//figure out which signs are fixed by p and facet_vecs
	p_signs := [Sign(p[i]) : i in [1..n+1]];
	
	sign_flip_indices := [];
	for i in [1..n] do
		if p_signs[i] ne 0 then			
			if Sign(normal[i]) eq -p_signs[i] then //can't tell for now, so have to keep checking
				Append(~sign_flip_indices, i);
			end if;
		else
			if Sign(normal[i]) ne 0 then
				p_signs[i] := Sign(normal[i]); //constructed perfect vector will have the same sign as the normal vector
			else
				p_signs[i] := 1;
				Append(~sign_flip_indices, i); //constructed perfect vector will be zero in the ith component
			end if;
		end if;
	end for;
	
	//now try to attach a new minimal vector
	bound := Infinity();
	height := 1;
	min_rho := Infinity();
	
	while height lt bound do		
		if height gt #boundary_point_database then
			Append(~boundary_point_database, boundaryPoints(height));
		end if;
		
		for v in boundary_point_database[height] do
			w := V ! [v[i] * p_signs[i] : i in [1..n+1]];
			
			for vec in signOrbit(w, sign_flip_indices) do
				if not minkowskiNorm(vec) eq 0 then
					print v;
				end if;
				normal_component := minkowskiForm(vec, normal);
				
				if normal_component lt 0 then //vec lies on the correct side of the facet
					rho := (1 - minkowskiForm(vec, p)) / normal_component;
					
					if rho lt min_rho then //possible minimal vector of neighbour
						candidate := p + rho * normal;
						
						if isInteriorPoint(candidate) then
							candidate_min, candidate_min_vecs, boundary_point_database := minimalVectors(candidate, boundary_point_database, 1, facet_vecs cat [vec]);
							
							if candidate_min eq 1 then //neighbour found
								return candidate, candidate_min_vecs, boundary_point_database;
							end if;
							
							new_bound := 2 * candidate[n+1] / minkowskiNorm(candidate);//not clear if the bound is increasing with rho, so only updating bound when definitely valid to do so
							if new_bound lt bound then
								bound := new_bound;
							end if;
						end if;
						
						min_rho := rho;
						
						for i in [1..#sign_flip_indices] do
							if Sign(candidate[sign_flip_indices[i]]) eq p_signs[sign_flip_indices[i]] then //sign of neighbour agrees with p in coordinate sign_flip_indices[i]; no longer need to flip
								Remove(~sign_flip_indices, i);
							end if;
						end for;
					end if;
				end if;
			end for;
		end for;
		
		height +:= 1;
	end while;
end function;


//----------Voronoi algorithm----------
function voronoi(boundary_point_database)
	initial_point, initial_min_vecs, boundary_point_database := buildInitialPoint(boundary_point_database);
	
	print "Class 1 with ", #initial_min_vecs, "minimal vectors";
	print initial_point;
	
	perfect_point_list := [initial_point];
	minimal_vectors_list := [initial_min_vecs];
	
	next_untested := 1; //index in perfect_point_list of the next class of perfect vectors that have not had their neighbours checked
	while next_untested le #perfect_point_list do
		facet_list := facets(minimal_vectors_list[next_untested]);
		
		print #facet_list, "facets of", next_untested;
		count := 1;
		for facet in facet_list do
			print "facet", count;
			new_neighbour, neighbour_min_vecs, boundary_point_database := neighbour(perfect_point_list[next_untested], minimal_vectors_list[next_untested], facet, boundary_point_database);
			
			equiv := false;
			for i in [1..#perfect_point_list] do				
				if equivalent(new_neighbour, neighbour_min_vecs, perfect_point_list[i], minimal_vectors_list[i]) then
					equiv := true;
					break;
				end if;
			end for;
			
			if not equiv then //new class found
				Append(~perfect_point_list, new_neighbour);
				Append(~minimal_vectors_list, neighbour_min_vecs);
				print "Class", #perfect_point_list, "with", #neighbour_min_vecs, "minimal vectors";
				print new_neighbour;
			end if;
			count +:= 1;
		end for;
		
		next_untested +:= 1;
	end while;
	
	return perfect_point_list, minimal_vectors_list;
end function;

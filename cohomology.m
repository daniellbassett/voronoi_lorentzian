//voronoi_lorentzian - Daniel Bassett, 2025
//--------------------------------------------------cohomology.m--------------------------------------------------
/*
TODO: write description
*/

import "congruence_modules.m" : coinduced_action_matrix, invariant_module, module_action;

//----------Group operations---------
function right_cosets(H, G) //finds a transversal for the right cosets of H in G
	transversal := [];
	coset_union := [];
	
	for g in G do
		if g notin coset_union then
			Append(~transversal, g);
			
			for h in H do
				Append(~coset_union, h*g);
			end for;
		end if;
	end for;
	
	return transversal;
end function;


function intersection(H, K) //finds the intersection of two subgroups
	cap := [];
	
	for h in H do
		if h in K then
			Append(~cap, h);
		end if;
	end for;
	
	return cap;
end function;

//----------Cochains----------
function cochain_modules(level, module, orientable_cell_representatives, representative_stabilisers)
	cochains_per_cell := [];
	//cochains := [];
	
	for dim in [1..#orientable_cell_representatives] do
		Append(~cochains_per_cell, []);
		
		for generators in representative_stabilisers[dim] do
			Append(~cochains_per_cell[dim], invariant_module(level, module, generators));
		end for;
		
		//Append(~cochains, DirectSum(cochains_per_cell[dim]));
	end for;
	
	return /*cochains,*/ cochains_per_cell; ///*the full cochain in each dimension, and*/ the summands from the cells in that dimension
end function;

//----------Coboundary maps----------
function coboundary_map(level, module, orientable_reps_high, stabilisers_high, orientations_high, facets, reps_low, orientable_low_indices, facet_rep_indices, facet_rep_witnesses, facet_orientations, facet_stabilisers, facet_rep_compatibility, high_low_compatibility, /*stabiliser_cosets,*/ low_cochains, high_cochains)
	/*
		orientable_reps_high is a list of the orientable representatives of higher-dimensional cells
		stabilisers_high is a list of the stabilisers of the above representatives
		orientations_high are a basis for the vector spaces spanned by the above representatives
		
		facets is a list of lists of *all* facets of each of the above representatives
		reps_low is a list of representatives for the lower-dimensional cells (not necessarily orientable?)
		facet_rep_indices is a list of which representative the given facets are equivalent to
		facet_rep_witnesses is a list of witnesses to the above equivalences
		
		orientable_low_indices is a list of which low-dimensional representatives are orientable
		facet_orientations are a basis for the vector spaces spanned by each facet
		facet_stabilisers are the stabiliser of each facet
		facet_rep_compatibility is whether the orientation of the facet is compatible with its equivalent representative
		high_low_compatibility is whether the orientation of the facet is compatible with its higher-dimensional cell
		stabiliser_cosets is a transversal for (Gamma_sigma intersect Gamma_{tau_prime}) \ Gamma_sigma
		
		facet data is stored first by (co)dimension, and then by which higher-dimensional representative it is a facet of
	*/	
	
	total_matrix := <>;
	
	for i in [1..#orientable_reps_high] do //sigma = orientable_reps_high[i] is the highdim cell
		sigma_matrix := <RMatrixSpace(Rationals(), Dimension(low_cochains[j]), Dimension(high_cochains[i])) ! 0 : j in [1..#orientable_low_indices]>;
		for j in [1..#facets[i]] do
			tau_prime := facets[i][j]; //facet of sigma
			
			if facet_rep_indices[i][j] in orientable_low_indices then //tau (and hence tau_prime) orientable
				tau := reps_low[facet_rep_indices[i][j]]; //tau_prime equivalent to representative tau
				gamma := facet_rep_witnesses[i][j]; //tau gamma = tau_prime
				
				//calculate (Gamma_sigma intersect Gamma_{tau_prime}) \ Gamma_sigma 
				/*
				stabiliser_intersection := intersection(stabilisers_high[i], facet_stabilisers[i][j]);
				stabiliser_cosets := right_cosets(stabiliser_intersection, stabilisers_high[i]);
				
				print stabiliser_cosets;
				*/
				
				low_cochain_summand := Index(orientable_low_indices, facet_rep_indices[i][j]); //index of the summand in low_cochains that tau corresponds to
				coord_list := [];
				
				W := VectorSpace(Rationals(), Dimension(high_cochains[i]));
				
				for v in Basis(low_cochains[low_cochain_summand]) do
					image := high_cochains[i] ! 0; //summand of high_cochains corresponding to sigma
					
					/*
					if #stabiliser_cosets gt 1 then
						print "cosets:",#stabiliser_cosets;
					end if;
					
					for coset in stabiliser_cosets do
						image +:= v * facet_rep_compatibility[i][j] * high_low_compatibility[i][j] * coinduced_action_matrix(level, module, gamma * coset);
					end for;
					*/
					
					image +:= v * facet_rep_compatibility[i][j] * high_low_compatibility[i][j] * coinduced_action_matrix(level, module, gamma);
					
					//add into the matrix using the chosen bases of the cochains
					Append(~coord_list, W ! Coordinates(high_cochains[i], image));
				end for;
				
				sigma_matrix[facet_rep_indices[i][j]] := VerticalJoin(coord_list);
			end if;
		end for;
		
		Append(~total_matrix, VerticalJoin(sigma_matrix));
	end for;
	
	return HorizontalJoin(total_matrix);
end function;

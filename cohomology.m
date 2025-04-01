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
			Append(~traversal, g);
			
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
function cochain_modules(level, module, orientable_cell_representatives, representative_stabilisers, invariant_module)
	cochains_per_cell := [];
	//cochains := [];
	
	for dim in [1..#orientable_cell_representatives] do
		Append(~cochains_per_cell, []);
		
		for generators in representative_stabilisers do
			Append(~cochains_per_cell[dim], invariant_module(level, module, generators));
		end for;
		
		//Append(~cochains, DirectSum(cochains_per_cell[dim]));
	end for;
	
	return /*cochains,*/ cochains_per_cell; ///*the full cochain in each dimension, and*/ the summands from the cells in that dimension
end function;

//----------Coboundary maps----------
function coboundary_map(level, module, orientable_reps_high, stabilisers_high, orientations_high, facets, reps_low, orientable_low_indices, orientability_low, stabilisers_low, facet_rep_indices, facet_rep_witnesses, facet_orientations, low_cochains, high_cochains)
	total_matrix := [];
	
	for i in [1..#orientable_reps_high] do //sigma = orientable_reps_high[i] is the highdim cell
		sigma_matrix := [MatrixSpace(Rationals(), Dimension(low_cochains[j]), Dimension(high_cochains[i])) ! 0 : j in [1..#orientable_low_indices]];
		for j in [1..#facets[i]] do
			tau_prime := facets[i][j]; //facet of sigma
			
			if orientability_low[facet_rep_indices[i][j]] then //tau_prime orientable
				tau := reps_low[facet_rep_indices[i][j]]; //tau_prime equivalent to representative tau
				gamma := facet_rep_witnesses[i][j]; //tau gamma = tau_prime
				
				//calculate (Gamma_sigma intersect Gamma_{tau_prime}) \ Gamma_sigma 
				stabiliser_intersection := intersection(stabilisers_high[i], stabilisers_low[facet_rep_indices[i][j]]);
				stabiliser_cosets := right_cosets(stabiliser_intersection, stabilisers_high[i]);
				
				low_cochain_summand := orientable_low_indices[facet_rep_indices[i][j]]; //index of the summand in low_cochains that tau corresponds to
				coord_list := [];
				for v in Basis(low_cochains[low_cochain_summand]) do
					image := high_cochains[i] ! 0; //summand of high_cochains corresponding to sigma
					
					for coset in stabiliser_cosets do
						image +:= module_action(level, module, gamma * coset);
					end for;
					
					//add into the matrix using the chosen bases of the cochains
					Append(~coord_list, Coordinates(high_cochains[i], image));
				end for;
				
				sigma_matrix[facet_rep_indices[i][j]] := VerticalJoin(coord_list);
			end if;
		end for;
		
		Append(~total_matrix, VerticalJoin(sigma_matrix));
	end for;
	
	return HorizontalJoin(total_matrix);
end function;

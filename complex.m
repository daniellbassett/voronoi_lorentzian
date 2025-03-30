//voronoi_lorentzian - Daniel Bassett, 2025
//--------------------------------------------------complex.m--------------------------------------------------
/*
	Calculates from an ideal tessellation of hyperbolic space a corresponding dual tessellation with no ideal vertices.
	Furthermore, hyperbolic space retracts by a strong deformation retraction onto this dual tessellation.
*/

import "symmetric_space.m" : facets, equivalent, barycentre;

//----------Full tessellation from top-cells----------
function tessellation(top_cells) //generates a full cell structure (*not* up to equivalence) given the top cells
	cells := [[cell : cell in top_cells]];
	
	//generate cells
	cell_vertex_count := #cells[#cells][1];
	while cell_vertex_count gt 1 do
		Append(~cells, []);
		for cell in cells[#cells-1] do //find facets of all cells of current lowest calculated dimension
			cell_facets := facets(cell);
			
			for facet in cell_facets do
				if not facet in cells[#cells] then
					Append(~cells[#cells], facet);
				end if;
			end for;
		end for;
	end while;
	
	cell_facet_indices := [];
	codim := 1; //actually codim+1
	
	//calculate cell inclusions
	while codim lt #cells do
		Append(~cell_facet_indices, []);
		
		for high_cell in cells[codim] do
			Append(~cell_facet_indices[#cell_facet_indices], []);
			
			for i in [1..#cells[codim+1]] do
				if IsSubsequence(cells[codim+1][i], high_cell : Kind := Setwise) then
					Append(~cell_facet_indices[#cell_facet_indices], i);
				end if;
			end for;
		end for;
	end while;
	
	return cells, cell_facet_indices;
end function;

function tessellation_representatives(top_cells)
	cell_reps := [[]];
	
	//find representatives for top dimension
	for cell in top_cells do
		new_class := true; //current cell is not equivalent to any previously checked top cells
		
		for rep in cell_reps do
			equiv, _ := equivalent(rep, cell);
			
			if equiv then
				new_class := false;
				break;
			end if;
		end for;
		
		if new_class then
			Append(~cell_reps[1], cell);
		end if;
	end for;
	
	
	cell_facets := [];
	facet_equiv_indices := [];
	facet_equiv_witnesses := [];
	
	cell_vertex_count := #cell_facets[1][1];
	while cell_vertex_count gt 1 do
		//generate facets of previous dimension's representatives
		Append(~cell_facets, []);
		
		for rep in cell_reps[#cell_reps] do
			new_facets := facets(rep);
			Append(~cell_facets[#cell_facets], new_facets);
		end for;
		
		//find representatives for the cells in the new dimension
		Append(~cell_reps, []);
		Append(~facet_equiv_indices, []);
		Append(~facet_equiv_witnesses, []);
		
		for higher_cell_facets in cell_facets[#cell_facets] do
			Append(~facet_equiv_indices[#facet_equiv_indices], []);
			Append(~facet_equiv_witnesses[#facet_equiv_indices], []);
			
			for facet in higher_cell_facets do
				if facet in cell_reps then
					Append(~facet_equiv_indices[#facet_equiv_indices][#facet_equiv_indices[#facet_equiv_indices]], i);
					Append(~facet_equiv_witnesses[#facet_equiv_witnesses][#facet_equiv_witnesses[#facet_equiv_witnesses]], MatrixRing(Rationals(), n+1) ! 1);
				else
					new_class := true;
					
					for i in [1..#cell_reps[#cell_reps]] do
						equiv, equivBy := equivalent(cell_reps[#cell_reps][i], facet);
						
						if equiv then
							new_class := false;
							
							Append(~facet_equiv_indices[#facet_equiv_indices][#facet_equiv_indices[#facet_equiv_indices]], i);
							Append(~facet_equiv_witnesses[#facet_equiv_witnesses][#facet_equiv_witnesses[#facet_equiv_witnesses]], equivBy);
							
							break;
						end if;
					end for;
					
					if new_class then
						Append(~cell_reps[#cell_reps], facet);
						
						Append(~facet_equiv_indices[#facet_equiv_indices][#facet_equiv_indices[#facet_equiv_indices]], #cell_reps[#cell_reps]);
						Append(~facet_equiv_witnesses[#facet_equiv_witnesses][#facet_equiv_witnesses[#facet_equiv_witnesses]], MatrixRing(Rationals(), n+1) ! 1);
					end if;
				end if;
			end for;
		end for;
		
		cell_vertex_count := #cell_facets[#cell_facets][1];
	end while;
	
	return cell_reps, cell_facets, facet_equiv_indices, facet_equiv_witnesses;
end function;

//----------Deformation retract----------
function retract_facets(cells, cell_facet_indices)
	vertices := [[barycentre(cell) : cell in cells[i]] : i in [1..#cells]]; //vertices of barycentric subdivision
	
	retract_top_cells := [];
	
	index_vector := [1 : i in [1..#cell_facet_indices]]; //flag of topdim, codim 1, ..., dim 2, dim 1. no dim 0 since we want cells without ideal vertices
	while index_vector[1] le #cells[1] do
		flag_indices := [index_vector[1]];
		while #flag_indices lt #cell_facet_indices do
			Append(~flag_indices, cell_facet_indices[#flag_indices][flag_indices[#flag_indices]][index_vector[#flag_indices+1]]);
		end while;
		
		Append(~retract_top_cells, [vertices[j][flag_indices[j]] : j in [1..#flag_indices]]); //top cell from flag
		
		//next top cell
		index_vector[#index_vector] +:= 1;
		for i in [#index_vector..2 by -1] do
			if index_vector[i] gt #cell_facet_indices[i-1][flag_indices[i-1]] then
				index_vector[i] := 1;
				index_vector[i-1] +:= 1;
			end if;
		end for;
	end while;
	
	return retract_top_cells;
end function;

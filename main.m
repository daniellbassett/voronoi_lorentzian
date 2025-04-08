//voronoi_lorentzian - Daniel Bassett, 2025
//--------------------------------------------------main.m--------------------------------------------------
/*	
TODO: write description
*/

import "voronoi.m" : voronoi;
import "complex.m" : tessellation, retractFacets, tessellationRepresentatives;
import "symmetric_space.m" : orientable, cellBasis, stabiliser, barycentre, compatibleOrientation, compatibleInducedOrientation;
import "congruence_modules.m" : coinduced_module;
import "cohomology.m" : cochain_modules, coboundary_map;

//----------Generate complex----------
//voronoi tessellation
perfect_points, perfect_cones := voronoi([]);
voronoi_tessellation_cells, voronoi_cell_inclusions := tessellation(perfect_cones); //generates the full voronoi tessellation for use in the well-rounded retract

print "voronoi complete";

//well-rounded retract cells up to equiv
retract_top_cells := retractFacets(voronoi_tessellation_cells, voronoi_cell_inclusions);
print "retract top cells generated";
cell_reps, cell_facets, facet_equiv_indices, facet_equiv_witnesses := tessellationRepresentatives(retract_top_cells);


//stabilisers and orientations of representatives
orientable_reps := []; //the reps that are orientable
orientability_indices := [];
rep_stabilisers := [];
rep_orientations := [];

for codim in [1..#cell_reps] do
	Append(~orientable_reps, []);
	Append(~orientability_indices, []);
	Append(~rep_stabilisers, []);
	Append(~rep_orientations, []);
	
	for rep in cell_reps[codim] do
		orientability, orientation, stab := orientable(rep, []);
		
		if orientability then
			Append(~orientable_reps[#orientable_reps], rep);
			Append(~orientability_indices[#orientability_indices], Index(cell_reps[codim], rep));
			Append(~rep_stabilisers[#rep_stabilisers], stab);
			Append(~rep_orientations[#rep_orientations], orientation);
		end if;
	end for;
end for;


//stabilisers and orientations of facets - truly the worst code you have ever seen
facet_stabilisers := [];
facet_orientations := [];
facet_rep_compatibility := []; //compatibility with orientation of equivalent representative
facet_high_compatibility := []; //compatibility with orientation of higher-dim cell

for codim in [1..#cell_facets] do
	Append(~facet_stabilisers, []);
	Append(~facet_orientations, []);
	Append(~facet_rep_compatibility, []);
	Append(~facet_high_compatibility, []);
	
	for i in [1..#cell_facets[codim]] do //data sorted per cell of higher dimension
		Append(~facet_stabilisers[#facet_stabilisers], []);
		Append(~facet_orientations[#facet_orientations], []);
		Append(~facet_rep_compatibility[#facet_rep_compatibility], []);
		Append(~facet_high_compatibility[#facet_high_compatibility], []);
		
		for j in [1..#cell_facets[codim][i]] do
			Append(~facet_stabilisers[#facet_stabilisers][#facet_stabilisers[#facet_stabilisers]], stabiliser(barycentre(cell_facets[codim][i][j]), cell_facets[codim][i][j]));
			
			facet_orientation := cellBasis(cell_facets[codim][i][j], Rank(VerticalJoin(cell_facets[codim][i][j])));
			Append(~facet_orientations[#facet_orientations][#facet_orientations[#facet_orientations]], facet_orientation);
			
			if facet_equiv_indices[codim][i][j] in orientability_indices[codim+1] then
				//compatibility with representative
				if compatibleOrientation(rep_orientations[codim+1][facet_equiv_indices[codim][i][j]], facet_orientation, facet_equiv_witnesses[codim][i][j]) then
					Append(~facet_rep_compatibility[#facet_rep_compatibility][#facet_rep_compatibility[#facet_rep_compatibility]], 1);
				else
					Append(~facet_rep_compatibility[#facet_rep_compatibility][#facet_rep_compatibility[#facet_rep_compatibility]], -1);
				end if;
				
				//compatibility with higherdim cell
				if i in orientability_indices[codim] then
					if compatibleInducedOrientation(facet_orientation, rep_orientations[codim][i]) then
						Append(~facet_high_compatibility[#facet_high_compatibility][#facet_high_compatibility[#facet_high_compatibility]], 1);
					else
						Append(~facet_high_compatibility[#facet_high_compatibility][#facet_high_compatibility[#facet_high_compatibility]], -1);
					end if;
				else //highdim cell not orientable
					Append(~facet_high_compatibility[#facet_high_compatibility][#facet_high_compatibility[#facet_high_compatibility]], 0);
				end if;
			else //facet not orientable
				Append(~facet_rep_compatibility[#facet_rep_compatibility][#facet_rep_compatibility[#facet_rep_compatibility]], 0);
				Append(~facet_high_compatibility[#facet_high_compatibility][#facet_high_compatibility[#facet_high_compatibility]], 0);
			end if;
		end for;
	end for;
end for;

//----------Cohomology calculations----------

//profiling
SetProfile(true);

L := 2; //level
cohomology_degrees := [2]; //degrees to calculate cohomology in; TODO: move to init

while true do
	if IsPrime(L) then //for now only dealing with prime levels
		print "";
		print "level",L;
		//construction of appropriate modules
		level := FiniteField(L, 1);
		M := coinduced_module(level);
		print "module has size", #M;
		
		cochains_per_rep := cochain_modules(level, M, orientable_reps, rep_stabilisers);
		
		//required coboundary maps
		coboundary_degrees := []; //degree is the dimension of the higher-dim cells used
		for d in cohomology_degrees do
			if d gt 0 then
				if d notin coboundary_degrees then
					Append(~coboundary_degrees, d);
				end if;
			end if;
			
			if d+1 lt #cell_reps then
				Append(~coboundary_degrees, d+1);
			end if;
		end for;
		
		//coboundary map calculation
		coboundary_codims := [#cell_reps - d : d in coboundary_degrees];
		print "calculating coboundaries in codims", coboundary_codims;
		
		coboundary_maps := [**];
		
		for codim in coboundary_codims do
			map := coboundary_map(level, M, orientable_reps[codim], rep_stabilisers[codim], rep_orientations[codim], cell_facets[codim], cell_reps[codim+1], orientability_indices[codim+1], facet_equiv_indices[codim], facet_equiv_witnesses[codim], facet_orientations[codim], facet_stabilisers[codim], facet_rep_compatibility[codim], facet_high_compatibility[codim], cochains_per_rep[codim+1], cochains_per_rep[codim]);
			
			Append(~coboundary_maps, map);
		end for;
		
		//cohomology calculation
		for d in cohomology_degrees do
			if d eq 0 then
				print "h ^ 0 =", Dimension(Kernel(coboundary_maps[1]));
			elif d+1 eq #cell_reps then
				print "h ^", d, "=", NumberOfColumns(coboundary_maps[#coboundary_maps]) - Rank(coboundary_maps[#coboundary_maps]);
				
				integral_matrix := ChangeRing(coboundary_maps[#coboundary_maps], Integers());
				print "torsion:", &*ElementaryDivisors(integral_matrix);
			else
				i := Index(coboundary_degrees, d);
				if not coboundary_maps[i] * coboundary_maps[i+1] eq 0 then
					print "error";
				end if;
				
				print "h ^", d, "=", Dimension(Kernel(coboundary_maps[i+1])) - Rank(coboundary_maps[i]);
				integral_matrix := ChangeRing(coboundary_maps[i], Integers());
				print "torsion:", &*ElementaryDivisors(integral_matrix);
			end if;
			
			
		end for;
		
		//profiling output
		G := ProfileGraph();
		ProfilePrintByTotalTime(G);
	end if;
	
	L +:= 1; //next level
end while;

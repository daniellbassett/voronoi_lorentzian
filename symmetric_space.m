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

import "init.m" : n, standard_form, q;//, multithread, num_threads;


//----------Lorentzian cone construction----------
//Ambient inner product space
if standard_form then //use [1, ..., 1] bilinear form
	B := MatrixRing(Rationals(), n+1) ! 1; //Identity matrix
else
	B := DiagonalMatrix(Rationals(), n+1, q cat [1]);
end if;
B_int := MatrixRing(Integers(), n+1) ! B; //for when working with integral types e.g. lattices

V := VectorSpace(Rationals(), n+1, B); //Equipped with bilinear form B

//Minkowski form
signature_matrix := DiagonalMatrix(Rationals(), [-1 : i in [1..n]] cat [1]); //signature n,1

J := B * signature_matrix;
J_int := MatrixRing(Integers(), n+1) ! J;

function minkowskiForm(v,w)
	return InnerProduct(v, w*signature_matrix);
end function;

function minkowskiNorm(v) //minkowskiForm(v,v)
	return InnerProduct(v, v*signature_matrix);
end function;

//Lorentzian cone
function isInteriorPoint(v)
	if v[n+1] lt 0 then
		return false;
	else
		return minkowskiNorm(v) gt 0;
	end if;
end function;

//----------Polytopes in Lorentz cone----------
function cellBasis(vertices, dim) //Finds a subset of vertices that is a basis for the space they span
	vertex_set := {v : v in vertices};
	for S in Subsets(vertex_set, dim) do
		if IsIndependent([s : s in S]) then
			return [s : s in S];
		end if;
	end for;
end function;

function cellNormal(vertices) //Finds a vector normal to all of vertices with respect to the minkowski form
	M := J * Transpose(VerticalJoin(vertices));
	orthBasis := KernelMatrix(M); //basis for orthogonal subspace
	
	return V ! RowSubmatrix(orthBasis, 1); //first element of basis as a vector
end function;

function facets(vertices)
	poly := Polytope([Eltseq(v) : v in vertices]);
	
	indices := FacetIndices(poly);
	facet_list := [];
	
	for facet in indices do
		Append(~facet_list, [vertices[i] : i in facet]);
	end for;
	
	return facet_list;
end function;

function barycentre(vertices)
	return (&+vertices)/#vertices;
end function;

//----------Integral boundary points----------
/*We iterate over the first n-1 coordinates, and check if this gives a possible solution for the nth coordinate,
given a fixed last coordinate ('height')*/

if q[n] eq 1 then //can avoid division logic, ~5% faster
	function boundaryPoints(height)
		V_low := VectorSpace(Rationals(), n, DiagonalMatrix(Rationals(), n, [-q[i] : i in [1..n-1]] cat [1]));
		
		points := [];
		
		v := V_low ! 0;
		v[n] := height;
		
		height_squared := height^2;
		last_coord_bound := q[n-1] * v[n-1]^2; //to check when the last coordinate being varied is too big
		
		while last_coord_bound le height_squared do
			norm := InnerProduct(v,v);
			square, coord := IsSquare(norm);
			
			if square then
				Append(~points, V ! ([v[i] : i in [1..n-1]] cat [coord, height]));
			end if;
			
			v[1] +:= 1;
			for i in [1..n-3] do
				if q[i] * v[i]^2 gt height_squared then
					v[i+1] +:= 1;
					v[i] := 0; //if this is too big then it just got increased, so the previous one has been set to 0; induction.
				end if;
			end for;
			
			if q[n-2] * v[n-2]^2 gt height_squared then //increasing last varied coord; remember to update last_coord_bound
				v[n-1] +:= 1;
				v[n-2] := 0;
				
				last_coord_bound := q[n]*v[n-1]^2;
			end if;
		end while;
		
		return points;
	end function;
else
	function boundaryPoints(height)
		V_low := VectorSpace(Rationals(), n, DiagonalMatrix(Rationals(), n, [-q[i] : i in [1..n-1]] cat [1]));
		
		points := [];
		
		v := V_low ! 0;
		v[n] := height;
		
		height_squared := height^2;
		last_coord_bound := q[n-1] * v[n-1]^2; //to check when the last coordinate being varied is too big
		
		while last_coord_bound le height_squared do
			norm := Integers() ! InnerProduct(v,v);
			denom := Integers() ! q[n];
			
			quotient, remainder := Quotrem(norm, denom);
			if remainder eq 0 then
				square, coord := IsSquare(quotient);
				
				if square then
					Append(~points, V ! ([v[i] : i in [1..n-1]] cat [coord, height]));
				end if;
			end if;

			v[1] +:= 1;
			for i in [1..n-3] do
				if q[i] * v[i]^2 gt height_squared then
					v[i+1] +:= 1;
					v[i] := 0; //if this is too big then it just got increased, so the previous one has been set to 0; induction.
				end if;
			end for;
			
			if q[n-2] * v[n-2]^2 gt height_squared then //increasing last varied coord; remember to update last_coord_bound
				v[n-1] +:= 1;
				v[n-2] := 0;
				
				last_coord_bound := q[n]*v[n-1]^2;
			end if;
		end while;
		
		return points;
	end function;
end if;

function signOrbit(v, sign_flip_indices) //orbit of v under a C_2 acting on each index in sign_flip_indices
	vectors := [];
	
	for i in [2^#sign_flip_indices .. 2^(#sign_flip_indices+1)-1] do
		bits := Prune(Intseq(i,2));
		signs := [2*bit - 1 : bit in bits];
		
		w := v;
		for j in [1..#sign_flip_indices] do
			w[sign_flip_indices[j]] *:= signs[j];
		end for;
		
		Append(~vectors, w);
	end for;
	
	return vectors;
end function;

//----------Equivalence testing----------
/*Testing for equivalence of Voronoi facets under the integral points of the isometry group
We may equivalently test for equivalence of the corresponding perfect vectors, and there is
a much faster method for this when standard_form is true, making use of lattice isometry testing.
*/

function naive_generators(G) //naive algorithm for finding a generating set of a rational matrix group of rank n+1
	generators := [];
	
	H := MatrixGroup<n+1, Rationals() | generators>;
	for g in G do
		if g notin H then
			Append(~generators, g);
			
			H := MatrixGroup<n+1, Rationals() | generators>;
			if #H eq #G then
				break;
			end if;
		end if;
	end for;
	
	return generators;
end function;

if standard_form then
	function posRep(v) //Associates a symmetric matrix to v; if minkowskiNorm(v) > 1/2 (e.g. for integral points, which we can get by rescaling) then this will be positive-definite.
		return (MatrixRing(Rationals(), n+1) ! Eltseq(TensorProduct(v,v))) - J/2;
	end function;
		
	function clearDenoms(v,w) //rescale consistently to integral vectors
		denoms := [];
		
		for i in [1..n+1] do
			Append(~denoms, Denominator(v[i]));
			Append(~denoms, Denominator(w[i]));
		end for;
		
		denom := LCM(denoms);
		v *:= denom;
		w *:= denom;
		
		return v,w;
	end function;
	 
	function equivalent(v, min_vecs_v, w, min_vecs_w)
		if #min_vecs_v eq #min_vecs_w then
			if minkowskiNorm(v) eq minkowskiNorm(w) then
				v,w := clearDenoms(v,w);
				
				L1 := LatticeWithGram(2*posRep(v));
				L2 := LatticeWithGram(2*posRep(w));
				
				equiv, equivBy := IsIsometric(L1, [J_int], L2, [J_int]); //Isometries of L1 and L2 that preserve J_int i.e. O(n,1; Z)
				if equiv then
					equivBy := MatrixRing(Rationals(), n+1) ! Transpose(equivBy); //tells us that v equiv^T = +/- w, with equiv in O(n,1; Z). if standard_form, this implies equiv^T in O^(n,1; Z)
					if v * equivBy eq w then
						return true, equivBy;
					else
						return true, -equivBy;
					end if;
				else
					return false, false;
				end if;
			else
				return false, false;
			end if;
		else
			return false, false;
		end if;
	end function;
	
	function stabiliser(v, min_vecs_v)
		v, _ := clearDenoms(v,v);
		
		L := LatticeWithGram(2*posRep(v));
		G := AutomorphismGroup(L, [J_int]);
		
		//find elements of G that genuinely fix v rather than just posRep(v)
		stab := [];
		for g in G do
			gamma := MatrixRing(Rationals(), n+1) ! Transpose(g);
			if v * gamma eq v then
				if Determinant(gamma) eq 1 then
					Append(~stab, gamma);
				end if;
			end if;
		end for;
		
		return naive_generators(stab); //TODO: find a way to treat this as a group within magma
	end function;
else
	function equivalent(v, min_vecs_v, w, min_vecs_w) //naive method: check all ordered subsets of n+1 elements of min_vecs_w and create the matrix
		if #min_vecs_v eq #min_vecs_w then
			if minkowskiNorm(v) eq minkowskiNorm(w) then //naive geometric invariants do not distinguish
				v_cell_basis := cellBasis(min_vecs_v, n+1);
				
				M := VerticalJoin(v_cell_basis);
				M_inv := M^-1;
				
				vertex_set_w := {p : p in min_vecs_w};
				sym_group := Sym({1..n+1});
				for S in Subsets(vertex_set_w, n+1) do
					subset_list := [s : s in S];
					for sigma in sym_group do
						sigma_seq := ElementToSequence(sigma);
						
						N := VerticalJoin([subset_list[sigma_seq[i]] : i in [1..n+1]]);
						gamma := M_inv * N;
						
						integral := true;
						for i in [1..n+1] do
							for j in [1..n+1] do
								if not IsIntegral(gamma[i,j]) then
									integral := false;
									break;
								end if;
							end for;
						end for;
						
						if integral then
							if gamma * J * Transpose(gamma) eq J then
								return true, gamma;
							end if;
						end if;
					end for;
				end for;
				
				return false, false; //no isometry found
			else
				return false, false;
			end if;
		else
			return false, false;
		end if;
	end function;
	
	function stabiliser(v, min_vecs_v);
		v_cell_basis := cellBasis(min_vecs_v, n+1);
		
		M := VerticalJoin(v_cell_basis);
		M_inv := M^-1;
		
		stab := [];
		
		vertex_set_v := {p : p in min_vecs_v};
		sym_group := Sym({1..n+1});
		for S in Subsets(vertex_set_v, n+1) do
			subset_list := [s : s in S];
			
			for sigma in sym_group do
				sigma_seq := ElementToSequence(sigma);
				
				N := VerticalJoin([subset_list[sigma_seq[i]] : i in [1..n+1]]);
				gamma := M_inv * N;
				
				integral := true;
				for i in [1..n+1] do
					for j in [1..n+1] do
						if not IsIntegral(gamma[i,j]) then
							integral := false;
							break;
						end if;
					end for;
				end for;
				
				if integral then
					if gamma * J * Transpose(gamma) eq J then
						Append(~stab, gamma);
					end if;
				end if;
			end for;
		end for;
		
		return stab; //TODO: find a way to treat this as a group within magma - look into MatrixGroup
	end function;
	
	function equivalentTwo(v, min_vecs_v, w, min_vecs_w) //makes use of (minkowski) inner products between elements of v_cell_basis. sometimes slower if there are very few different inner products
		if #min_vecs_v eq #min_vecs_w then
			if minkowskiNorm(v) eq minkowskiNorm(w) then //naive geometric invariants do not distinguish
				v_cell_basis := cellBasis(min_vecs_v, n+1);
				
				M := VerticalJoin(v_cell_basis);
				M_inv := M^-1;
				
				v_form_invariants := [minkowskiForm(v_cell_basis[1], v_cell_basis[i]) : i in [2..n+1]];
				
				
				for vec1 in min_vecs_w do
					w_form_invariants := [minkowskiForm(vec1, vec) : vec in min_vecs_w];
					
					possible := true;
					invariant_match_indices := [[] : i in [1..n]]; //find which minimal vectors of w give the correct inner product to match with those of v_cell_basis
					for i in [1..n] do
						invariant_match_indices[i] := [j : j in [1..#w_form_invariants] | w_form_invariants[j] eq v_form_invariants[i]];
						
						if #invariant_match_indices[i] eq 0 then
							possible := false;
							break;
						end if;
					end for;
					
					if not possible then
						continue;
					end if;
					
					//iterate over candidates
					multi_index := [1 : i in [1..n]];
					while multi_index[n] le #invariant_match_indices[n] do
						//construct
						N := VerticalJoin([vec1] cat [min_vecs_w[invariant_match_indices[i][multi_index[i]]] : i in [1..n]]);
						gamma := M_inv * N;
						
						//check
						integral := true;
						for i in [1..n+1] do
							for j in [1..n+1] do
								if not IsIntegral(gamma[i,j]) then
									integral := false;
									break;
								end if;
							end for;
						end for;
						
						if integral then
							if gamma * J * Transpose(gamma) eq J then
								return true, gamma;
							end if;
						end if;
						
						//next candidate
						multi_index[1] +:= 1;
						for i in [1..n-1] do
							if multi_index[i] gt #invariant_match_indices[i] then
								multi_index[i] := 1;
								multi_index[i+1] +:= 1;
							end if;
						end for;
					end while;
				end for;
				
				return false, false; //none found
			else
				return false, false;
			end if;
		else
			return false, false;
		end if;
	end function;
	
	//TODO: may be able to make use of all of the inner products, instead of just those with the first element
end if;

//----------Cell orientations----------
function orientable(vertices, stab) //checks if a cell is orientable under its stabiliser, and if so returns a choice of orientation
	v := barycentre(vertices);
	
	if #stab eq 0 then
		stab := stabiliser(v, vertices); //not sure whether this works outside of standard form - the alternative stabiliser routine probably doesn't work for lower dim cells
	end if;
	
	basis := cellBasis(vertices, Rank(VerticalJoin(vertices))); //determines an orientation
	
	V := VectorSpaceWithBasis(basis);
	W := VectorSpace(Rationals(), Dimension(V));
	for gamma in stab do
		images := [];
		
		for b in basis do
			Append(~images, W ! Coordinates(V, b*gamma));
		end for;
		
		M := VerticalJoin(images);
		if Sign(Determinant(M)) lt 0 then
			return false, false, false;
		end if;
	end for;
	
	return true, basis, stab;
end function;


function compatibleOrientation(basis1, basis2, gamma) //gamma takes basis1 to vertices in Span(basis2)
	acted_basis1 := [v * gamma : v in basis1];
	
	V := VectorSpace(Rationals(), #basis2);
	W := VectorSpaceWithBasis(basis2);
	images := [];
	for v in acted_basis1 do
		Append(~images, V ! Coordinates(W, v));
	end for;
	
	M := VerticalJoin(images);

	return Determinant(M) gt 0;
end function;


function compatibleInducedOrientation(vertices1, vertices2) //basis1 subset of Span(basis2); spans a space of codimension 1
	for v in vertices2 do
		if v notin vertices1 then
			extra_vertex := v;
			break;
		end if;
	end for;
	
	V := VectorSpace(Rationals(), #vertices2);
	W := VectorSpaceWithBasis(vertices2);
	
	coords := [];
	for v in (vertices1 cat [extra_vertex]) do
		Append(~coords, V ! Coordinates(W, v));
	end for;
	
	M := VerticalJoin(coords);
	
	return Determinant(M) gt 0;
end function;

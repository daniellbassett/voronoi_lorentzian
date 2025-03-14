//voronoi_lorentzian - Daniel Bassett, 2025
//--------------------------------------------------congruence_modules.m--------------------------------------------------
/*
TODO: write description
*/

import "init.m" : n, q;


//----------Coinduced module----------
/*For prime p, the coset space SO(n,1) / Gamma^0(p) has a model as the isotropic points in P^n(F_p)
We treat elements of the coinduced module as row vectors, indexed by these isotropic points.
Thus the coinduced_action_matrix A has a_ij = 1 when it takes basis vector i to basis vector j
*/

function projective_standard_form(v) //normalises so that the first non-zero coordinate is 1
	for i in [1..n+1] do
		if not v[i] eq 0 then
			return v / v[i];
		end if;
	end for;
end function;

function coinduced_module(p) //calculates the projective standard forms of the isotropic points in P^n(F_p)
	F_p := FiniteField(p, 1);
	V_p := VectorSpace(F_p, n+1, DiagonalMatrix(F_p, n+1, -q cat [1]));
	
	isotropic_points := [];
	for v in Vp do
		if v eq 0 then
			continue;
		end if;
		
		if projective_standard_form(v) eq v then //take only one point on each line
			if Norm(v) eq 0 then
				Append(~isotropic_points, v);
			end if;
		end if;
	end for;
	
	return isotropic_points;
end function;

function module_action(isotropic_point, gamma) //calculates the action of a group element gamma on a point isotropic_point
	return projective_standard_form(isotropic_point * gamma)
end function;

function coinduced_action_matrix(module, gamma) //calculates a permutation matrix representing the action of gamma on the coinduced module
	M := MatrixRing(Rationals(), #module) ! 0;
	
	for i in [1..#module] do //find the image of the ith basis vector
		image_vector := module_action(module[i], gamma);
		j := Index(module, image_vector);
		
		M[i,j] := 0;
	end for;
	
	return M;
end function;

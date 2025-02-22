//voronoi_lorentzian - Daniel Bassett, 2025
//--------------------------------------------------init.m--------------------------------------------------
/*
	Sets user-defined parameters: hyperbolic space dimension, quadratic form coefficients
n defines the dimension of the hyperbolic space. Recommended n <= 6 
q defines the coefficients of the quadratic form. Must be a list of n positive integers.
If standard_form is set to true, the value of q will be ignored, and replaced with a list of all 1s.
This allows the use of a faster algorithm for calculating orbits in this case.

	Sets practical parameters of the computation: multithreading settings and output filenames
One thread by default unless multithread is set to true, in which case num_threads processes are created
It is recommended that num_threads is at most the number of cores of your machine.
Multithreading currently only for "disgustingly parallel" portions of the computation e.g. finding boundary points.
In higher dimensions this tends to be a very significant proportion of computation time.
If save_outputs is set to true, the program will output reusable data to file.*/

//----------Hyperbolic lattice setup----------
n := 4;
standard_form := false;
q := [1, 1, 1, 1]; //Positive coefficients of quadratic form


//----------Computational settings----------
multithread := false;
num_threads := 1;


//----------Output filenames----------
save_outputs := true;
boundary_pt_file := "_boundary_points.txt" //For a database of integral boundary points of the Lorentzian cone
voronoi_data_file := "_voronoi_data.txt" //For representatives of the classes of perfect points, and their minimal vectors
complex_file := "_complex.txt" //For the cell structure of the Voronoi tessellation.


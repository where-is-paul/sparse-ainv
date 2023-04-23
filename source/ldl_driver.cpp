#include <solver.h>

#include <iostream>
#include <cassert>
#include <cstring>
#include <gflags.h>

DEFINE_string(filename, "", "The filename of the matrix to be factored"
		"(in matrix-market format).");

//DEFINE_double(fill, 3.0, "A parameter to control memory usage. Each column is guaranteed"
//		"to have fewer than fill*nnz(A) elements.");

DEFINE_double(tol, 0.001, "A parameter to control agressiveness of dropping. If 'relative', then in each column k,"
		"elements less than tol*||L(k+1:n,k)|| (1-norm) are dropped. If 'absolute', each element less than tol is dropped.");

DEFINE_double(beta, 1.0, "A parameter to aggressiveness of Bunch-Kaufman pivoting (BKP). "
		"When pp_tol >= 1, full BKP is used. When pp_tol is 0, BKP is faster"
		"but chooses poorer pivots. Values between 0 and 1 varies the aggressiveness of"
		"BKP in a continuous manner.");

DEFINE_string(pivot, "bunch", "Determines what kind of pivoting algorithm will be used"
		" during the factorization. Choices are 'wmn' and 'bunch'. The default is 'bunch'.");

DEFINE_string(reordering, "amd", "Determines what sort of preordering will be used"
		" on the matrix. Choices are 'amd', 'rcm', 'metis', 'mc64', and 'none'.");

DEFINE_string(equil, "bunch", "Decides if the matrix should be equilibriated before factoring is done. "
		"Options are 'bunch' and 'none'. If the option is 'bunch', the matrix is equilibrated "
		"with Bunch's algorithm in the max norm. The default is 'bunch'.");

DEFINE_bool(save, true, "If yes, saves the factors (in matrix-market format) into a folder "
		"called output_matrices/ in the same directory as ldl_driver.");

#ifdef SYM_ILDL_DEBUG
DEFINE_bool(display, false, "If yes, outputs a human readable version of the factors onto"
		" standard out. Generates a large amount of output if the "
		"matrix is big.");
#endif

DEFINE_int32(max_iters, -1, "If >= 0 and supplied with a right hand side, SYM-ILDL will attempt "
		"to use the preconditioner generated with the chosen solver to "
		"solve the system. This parameter controls the max iterations of "
		"the solver (has no effect if choosing the full solve).");

DEFINE_string(solver, "sqmr", "The solver used if supplied a right-hand side. The "
		"solution will be written to output_matrices/ in matrix-market "
		"format. Choices are 'sqmr', 'minres', and 'full'");

DEFINE_double(solver_tol, 1e-6, "A tolerance for the iterative solver used. When the iterate x satisfies ||Ax-b||/||b|| < solver_tol, the solver is terminated. Has no effect when doing a full solve.");

DEFINE_string(rhs_file, "", "The filename of the right hand side (in matrix-market format).");

DEFINE_string(drop_type, "absolute", "The type of drop tolerance used for removing small elements. Options are 'relative' or 'absolute'.");

int main(int argc, char* argv[])
{
	std::string usage("Performs an incomplete LDL factorization of a given matrix.\n"
			"Sample usage:\n"
			"\t./ldl_driver -filename=test_matrices/testmat1.mtx "
			"-display=true -save=false\n"
			"Additionally, these flags can be loaded from a single file "
			"with the option -flagfile=[filename].");

	google::SetUsageMessage(usage);
	google::ParseCommandLineFlags(&argc, &argv, true);

	if (FLAGS_filename.empty()) {
		std::cerr << "No file specified! Type ./ldl_driver --help for a description of the program parameters." << std::endl;
		return 0;
	}

	typedef long double precision_type;
	sparse_ainv::solver<precision_type> solv;

	//default is statistics output
	solv.set_message_level("statistics");

	//load matrix
	solv.load(FLAGS_filename);

	//default reordering scheme is AMD
	solv.set_reorder_scheme(FLAGS_reordering.c_str());

	//default is equil on
	solv.set_equil(FLAGS_equil.c_str()); 

	//default solver is SQMR
	solv.set_solver(FLAGS_solver.c_str());
	if (FLAGS_max_iters > 0 || !FLAGS_rhs_file.empty()) {
		if (FLAGS_max_iters <= 0) {
			printf("Using SQMR (200 max iterations) as default solver since RHS was loaded.\n");
			FLAGS_max_iters = 200;
		}

		vector<precision_type> rhs;
		if (!FLAGS_rhs_file.empty()) {
			sparse_ainv::read_vector(rhs, FLAGS_rhs_file);
		} else {
			// for testing purposes only
			rhs.resize(solv.A.n_cols(), 1);
		}

		if (rhs.size() != solv.A.n_cols()) {
			std::cout << "The right hand side dimensions do not match the dimensions of A." << std::endl;
			return 1;
		}
		solv.set_rhs(rhs);
	}

	solv.set_drop_type(FLAGS_drop_type.c_str());
	solv.set_pivot(FLAGS_pivot.c_str());
	sparse_ainv::solver_params par;
	par.ainv_tol = FLAGS_tol;
	par.ainv_beta = FLAGS_beta;
	par.max_iters = FLAGS_max_iters;
	par.solver_tol = FLAGS_solver_tol;

	solv.solve(par);

	if (FLAGS_save) {
		solv.save();
	}

#ifdef SYM_ILDL_DEBUG
	if (FLAGS_display) {
		solv.display();
		std::cout << endl;
	}
#endif

	std::cout << "Factorization Complete. ";
	if (FLAGS_save) {
		std::cout << "All output written to /output_matrices directory.";
	}
	std::cout << std::endl;

	return 0;
}

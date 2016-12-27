#ifndef _LILC_MATRIX_AINV_H_
#define _LILC_MATRIX_AINV_H_

#include <default_pivoter.h>
#include <bkp_pivoter.h>
#include <set_unioner.h>

#include <map>

#if 0
#include <float.h>
unsigned int fp_control_state = _controlfp(_EM_INEXACT, _MCW_EM);
#endif

using std::endl;
using std::cout;
using std::abs;

template <class el_type>
void lilc_matrix<el_type> :: ainv(lilc_matrix<el_type>& L, block_diag_matrix<el_type>& D, idx_vector_type& perm, params par) {
	//---------------------------- initialize temporary variables ------------------------------//
	const int ncols = n_cols(); //number of cols in A.

	elt_vector_type work(ncols, 0);
	idx_vector_type curr_nnzs;
	curr_nnzs.reserve(ncols); //reserves space for worse case (entire col is non-zero)
	
	vector<int> p(ncols, 0), pinv(ncols, 0);
	for (int k = 0; k < ncols; k++) {
		p[k] = k;
		pinv[k] = k;
	}

	int count = 0; //the total number of nonzeros stored in L.

	double normA = 0; // frobenius norm of A
	// Ensure nnz pattern of A is sorted
	vector<std::pair<int, el_type>> colval;
	for (int k = 0; k < ncols; k++) {
		colval.clear();
		for (int j = 0; j < m_x[k].size(); j++) {
			colval.push_back(std::make_pair(m_idx[k][j], m_x[k][j]));
		}
		sort(colval.begin(), colval.end());

		for (int j = 0; j < m_x[k].size(); j++) {
			m_idx[k][j] = colval[j].first;
			m_x[k][j] = colval[j].second;
		}

		double tnorm = norm(m_x[k], 2.0);
		normA += tnorm * tnorm;
	}
	normA = sqrt(normA);


	//------------------------------- allocate memory for L and D -------------------------------//
	L.resize(ncols, ncols); //allocate a vector of size n for Llist as well
	D.resize(ncols);

	// Initialize L to be identity
	for (int k = 0; k < ncols; k++) {
		L.m_x[k].push_back(1.0);
		L.m_idx[k].push_back(k);
		L.m_list[k].insert(k);
	}

	// choose pivoting strategy
	pivot_strategy<el_type>* pivoter;
	if (par.piv_type == pivot_type::BKP) {
		pivoter = new bkp_pivoter<el_type>(this, &L, &p, &pinv, par.beta);
	} else {
		pivoter = new default_pivoter<el_type>(this, &L, &p, &pinv);
	}
	pivoter->set_regularization(normA);

	// Some structs for set union operations
	vector<bool> seen1(ncols, 0);
	set_unioner seen0;
	seen0.reset(ncols);

	// where the index value stored is > k
	vector<int> A1_idx, Ar_idx, pvt_idx;
	vector<el_type> A1(ncols), Ar(ncols);

	// ------------ functions for ordering by number of non zeros -------------//
	// Tracks number of non-zeros per col so we can sort by sparsity
	vector<int> num_nz(ncols, 1);
	auto comp = [&](const int i, const int j) {
		if (num_nz[i] == num_nz[j]) return p[i] > p[j];
		return num_nz[i] < num_nz[j];
	};

	set<int, decltype(comp)> q(comp);
	for (int i = 0; i < ncols; i++) {
		q.insert(i);
	}

	//------------------- convenience functions for pivoting ------------------//
	auto pivot_vec = [&](int k, int r, vector<el_type>& v, vector<int>& idx) {
		for (int i = 0; i < idx.size(); i++) {
			if (idx[i] == k) {
				idx[i] = r;
			} else if (idx[i] == r) {
				idx[i] = k;
			}
		}
		std::swap(v[k], v[r]);
	};

	auto pivot = [&](int k, int r) {
		//std::cerr << "swapping " << k << " and " << r << std::endl;
		if (k == r) return;
		q.erase(k);
		q.erase(r);
#if 0
		std::cerr << "-----------pre pivot---------- " << std::endl;
		// output m_list and make sure im not crazy
		for (int i = 0; i < ncols; i++) {
			for (int j : L.m_list[i]) {
				std::cerr << j << " ";
			}
			std::cerr << std::endl;
		}
		std::cerr << std::endl;
#endif

		// Swap values in L
		// Swap k and r values in L.m_list
		for (int j : L.m_idx[k]) {
			if (L.m_list[j].count(r)) {
				continue;
			}
			L.m_list[j].erase(k);
			L.m_list[j].insert(r);
		}

		for (int j : L.m_idx[r]) {
			if (L.m_list[j].count(k)) {
				continue;
			}
			L.m_list[j].erase(r);
			L.m_list[j].insert(k);
		}

		std::swap(L.m_list[k], L.m_list[r]);
		L.m_x[k].swap(L.m_x[r]);
		L.m_idx[k].swap(L.m_idx[r]);
		// Because nnz indices are sorted, we can just change the last
		// elts indices since theyre guaranteed to be the diagonal
		L.m_idx[k][L.m_idx[k].size()-1] = k;
		L.m_idx[r][L.m_idx[r].size()-1] = r;

		// Swap values in A
		std::swap(p[k], p[r]);

		// Dont ask.
		//q.erase(p[k]);
		//q.erase(p[r]);
		std::swap(pinv[p[k]], pinv[p[r]]);

		// Swap values in perm vector
		std::swap(perm[k], perm[r]);

#if 0
		std::cerr << "-----------post pivot---------- " << std::endl;
		// output m_list and make sure im not crazy
		for (int i = 0; i < ncols; i++) {
			for (int j : L.m_list[i]) {
				std::cerr << j << " ";
			}
			std::cerr << std::endl;
		}
		std::cerr << "------------------------------- " << std::endl;
#endif

		num_nz[k] = static_cast<int>(L.m_idx[r].size());
		num_nz[r] = static_cast<int>(L.m_idx[k].size());
		
		q.insert(k);
		q.insert(r);
		//if (p[k] >= k) q.insert(p[k]);
		//if (p[r] >= k) q.insert(p[r]);
	};

	auto advance_list = [&](int k) {
		// increment row_first array if needed
		for (int j : L.m_idx[k]) {
			while (!L.m_list[j].empty() && *L.m_list[j].begin() <= k) {
				L.m_list[j].erase(L.m_list[j].begin());
			}
		}

		// drop elements of col k
		q.erase(k);
		drop_tol(L.m_x[k], L.m_idx[k], par.tol, k);
		
		// increase non-zero count of L
		count += static_cast<int>(L.m_x[k].size());
	};

	auto update_schur = [&](double coef, int j, int k) {
		//std::cerr << "altering " << j << std::endl;
		col_wrapper<el_type> lk(L.m_x[k].data(), L.m_idx[k].data(), L.m_x[k].size()),
							 lj(L.m_x[j].data(), L.m_idx[j].data(), L.m_x[j].size());

		if (coef != coef /*|| std::abs(coef) > normA*/) {
				std::cerr << j << " " << k << " wait wtffff bad coef " << coef << " " << normA << std::endl;
				//_sleep(100000);
			}

		sparse_vec_add(1.0, lj, coef, lk, work, curr_nnzs);

		for (auto x : work) {
			if (x != x) {
				//std::cerr << "wait wtffff bad computed valuees" << std::endl;
				break;
			}
		}

		for (int i = 0; i < (int)curr_nnzs.size()-1; i++) {
			if (curr_nnzs[i] > curr_nnzs[i+1]) {
				std::cerr << "waitwtfff" << std::endl;
			}
		}

		// Add existing nnz into set
		for (int i : L.m_idx[j]) {
			seen1[i] = true;
		}

		L.m_x[j].swap(work);
		L.m_idx[j].swap(curr_nnzs);

		// Add to m_list if needed
		for (int i : L.m_idx[j]) {
			if (!seen1[i]) {
				L.m_list[i].insert(j);
			}
			seen1[i] = false;
		}

		// Undo set additions
		for (int i : curr_nnzs) {
			if (seen1[i]) {
				L.m_list[i].erase(j);
			}
			seen1[i] = false;
		}

		q.erase(j);
		num_nz[j] = static_cast<int>(L.m_idx[j].size());
		q.insert(j);
	};

	auto apply_dropping_rules = [&](const vector<int>& pvts, int k) {
		// Figure out droptol
		double min_M_colsum = std::numeric_limits<double>::max();
		for (int i : pvts) {
			if (i <= k) continue;
			min_M_colsum = std::min(min_M_colsum, pow(norm<el_type>(L.m_x[i], 2.0), 2.0));
		}

		for (int i : pvts) {
			// Apply dropping rules to schur complement
			if (i <= k) continue;

			for (int j : L.m_idx[i]) {
				seen1[j] = true;
			}
			curr_nnzs = L.m_idx[i];
			drop_tol(L.m_x[i], L.m_idx[i], 0.1 * par.tol * min_M_colsum, i, true);
			
			for (int j : L.m_idx[i]) {
				seen1[j] = false;
			}
			for (int j : curr_nnzs) {
				if (seen1[j]) {
					seen1[j] = false;
					L.m_list[j].erase(i);
				}
			}
			curr_nnzs.clear();

			// Fix ordering
			q.erase(i);
			num_nz[i] = static_cast<int>(L.m_x[i].size());
			q.insert(i);
		}
	};

	//------------------- main loop: factoring begins -------------------------//
	for (int k = 0; k < ncols; k++) {
		//std::vector<int> rest(q.begin(), q.end());
		//int nk = *q.begin();
		//if (num_nz[nk] < 0.5 * num_nz[k]) {
		//	pivot(k, nk);
		//}

		static std::map<int, bool> done;
		int pcnt = (100*(k+1))/ncols;
		if (pcnt%10 == 0 && !done[pcnt]) {
			std::cerr << k << std::endl;
			done[pcnt] = 1;
			std::cerr << pcnt << " percent complete. ";
			int cnt = 0;
			for (int i = 0; i < ncols; i++) cnt += static_cast<int>(L.m_x[i].size());
			std::cerr << "nnz(M) = " << cnt << ", ";

			// Compute frobenius norm
			double fro = 0;
			for (int i = 0; i < ncols; i++) fro += pow(norm<el_type>(L.m_x[i], 2), 2);
			std::cerr << "Frobenius norm: " << sqrt(fro) << std::endl;
		}

		//std::cerr << "on iteration " << k << std::endl;
#if 0
		std::cerr << "on iteration " << k << std::endl;
		// output L and make sure im not crazy
		for (int i = 0; i < ncols; i++) {
			for (int j = 0; j < ncols; j++) {
				std::cerr << L.coeff(i, j) << " ";
			}
			std::cerr << std::endl;
		}
		std::cerr << std::endl << std::endl;
#endif

#if 0
		bool sane = true;
		for (int i = 0; i < ncols; i++) {
			for (int j = 0; j < int(L.m_idx[i].size())-1; j++) {
				if (L.m_idx[i][j] > L.m_idx[i][j+1]) {
					sane = false;
				}
			}
			for (int x : L.m_list[i]) {
				if (x < k) std::cerr << "something went terribly wrong" << std::endl;
			}
		}
		if (!sane) std::cerr << "sanity check: " << sane << endl;
#endif

#if 0
		bool sane = true;
		for (int i = 0; i < ncols; i++) {
			for (int j = 0; j < int(L.m_idx[i].size())-1; j++) {
				if (L.m_idx[i][j] > L.m_idx[i][j+1]) {
					sane = false;
				}
			}
			if (L.m_list[i].size() && *L.m_list[i].begin() < k) {
				sane = false;
				std::cerr << "list sanity failed!" << std::endl;
			}
		}
		assert(sane);
#endif

#if 0
			std::cerr << "Current matrix: " << std::endl;
			for (int i = 0; i < ncols; i++) {
				for (int j = 0; j < ncols; j++) {
					std::cerr << coeff(p[i], p[j]) << " ";
				}
				std::cerr << std::endl;
			}
			std::cerr << endl;
#endif

#if 0
		std::cerr << "on iteration " << k << std::endl;
		// output m_list and make sure im not crazy
		for (int i = 0; i < ncols; i++) {
			for (int j : L.m_list[i]) {
				std::cerr << j << " ";
			}
			std::cerr << std::endl;
		}
		std::cerr << std::endl;
#endif

		// need to pivot on ((A*M(:,p(dd)))'*M(:,p(dd:end)))';
		// or, M'(:, p(dd)) * (A M(:, p(dd:end))) = l_{pk}' * [A l_{pk} A l_{p{k+1}) ... A l_pn]
		// this is equivalent to pivoting on the D[j..n] vector
		pivot_struct piv_info = pivoter->find_pivot(k);

		//std::cerr << k << std::endl;
		//if (k == 7917) {
		//	std::cerr << "at problem point. pivot info: " << piv_info.r << " " << piv_info.size_two << std::endl;
		//}

		// swap necessary rows and columns
		if (piv_info.size_two) {
			//static int np2 = 0;
			//std::cerr << ++np2 << " size two pivots." << std::endl;
			
			// figure out list of D[k]'s to compute and update
			pivoter->flush_col(A1, A1_idx, 0);
			pivoter->flush_col(Ar, Ar_idx, 1);

			if (piv_info.r != k+1) {
				//static int np2 = 0;
				//std::cerr << ++np2 << " 2x2 pivots so far. " << std::endl;
				pivot(k+1, piv_info.r);
				pivot_vec(k+1, piv_info.r, A1, A1_idx);
				pivot_vec(k+1, piv_info.r, Ar, Ar_idx);

				// Swap indices in Ar and A1... this is not necessary in
				// the other case since A1 was guaranteed to have indices k
				// and r
				/*
				for (int& j : A1_idx) {
					if (j == piv_info.r) {
						j = k+1;
					} else if (j == k+1) {
						j = piv_info.r;
					}
				}

				for (int& j : Ar_idx) {
					if (j == piv_info.r) {
						j = k+1;
					} else if (j == k+1) {
						j = piv_info.r;
					}
				}*/
			}

#if 0
		std::cerr << "on iteration " << k << std::endl;
		// output L and make sure im not crazy
		for (int i = 0; i < ncols; i++) {
			for (int j = 0; j < ncols; j++) {
				std::cerr << L.coeff(i, j) << " ";
			}
			std::cerr << std::endl;
		}
		std::cerr << std::endl << std::endl;

		bool sane = true;
		for (int i = 0; i < ncols; i++) {
			for (int j = 0; j < int(L.m_idx[i].size())-1; j++) {
				if (L.m_idx[i][j] > L.m_idx[i][j+1]) {
					sane = false;
				}
			}
		}
		std::cerr << "sanity check: " << sane << endl;
#endif

			// merge non-zero indices, compute accumlation depending on even or odd
			seen0.add_set(A1_idx);
			seen0.add_set(Ar_idx);
			seen0.flush(pvt_idx);

#if 0
			for (int j : pvt_idx) {
				std::cerr << j << " ";
			}
			std::cerr << endl;
#endif

			// now i have the two columns.. merge their non-zero patterns
			el_type Di[2][2];
			Di[0][0] = A1[k]; 
			Di[0][1] = Di[1][0] = (A1[k+1] + Ar[k]) / 2;
			Di[1][1] = Ar[k+1];

			D[k] = Di[0][0];
			D.off_diagonal(k) = Di[0][1];
			D[k+1] = Di[1][1];

			//if (k == 7917) {
			//	for (int i = 0; i < 2; i++) {
			//		for (int j =0 ; j < 2; j++) {
			//			std::cerr << Di[i][j] << " ";
			//		}
			//		std::cerr << std::endl;
			//	}
			//}

#if 0
			for (int i = k; i < 10; i++) {
				std::cerr << A1[i] << " " << Ar[i] << endl;
			}

			std::cerr << "heres the diagonal: " << std::endl;
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					std::cerr << Di[i][j] << " ";
				}
				std::cerr << std::endl;
			}

			std::cerr << endl;
			assert(std::abs(A1[k+1] - Ar[k]) < 1e-6);
#endif

			for (int j : pvt_idx) {
				if (j <= k+1) continue;

				// Z(:, [j0, j1]) -= Z(:, [k, k+1]) * [a b; c d]
				//[a * zk + c * z_k+1, b * zk + d * zk+1]
				el_type DinvZ[2];
				el_type det = Di[0][0] * Di[1][1] - Di[0][1] * Di[1][0];
				DinvZ[0] =  A1[j] * Di[1][1] - Ar[j] * Di[0][1], 
				DinvZ[1] = -A1[j] * Di[1][0] + Ar[j] * Di[0][0];

				if (det == 0 || det != det) {
					std::cerr << "bad det" << std::endl;
					return;
				}

				update_schur(-DinvZ[0] / det, j, k);
				update_schur(-DinvZ[1] / det, j, k+1);

#if 0
		std::cerr << "on iteration " << k << std::endl;
		// output L and make sure im not crazy
		for (int i = 0; i < ncols; i++) {
			for (int j = 0; j < ncols; j++) {
				std::cerr << L.coeff(i, j) << " ";
			}
			std::cerr << std::endl;
		}
		std::cerr << std::endl << std::endl;

		bool sane = true;
		for (int i = 0; i < ncols; i++) {
			for (int j = 0; j < int(L.m_idx[i].size())-1; j++) {
				if (L.m_idx[i][j] > L.m_idx[i][j+1]) {
					sane = false;
				}
			}
		}
		std::cerr << "sanity check: " << sane << endl;
#endif
			}
			apply_dropping_rules(pvt_idx, k+1);

			// Extra iter accounted for
			advance_list(k);
			advance_list(k+1);
			k++;
		} else {
			// figure out list of D[k]'s to compute and update
			pivoter->flush_col(A1, A1_idx, 0);

			// do swaps for pivots if necessary
			if (piv_info.r != k) {
				//static int np1 = 0;
				//std::cerr << ++np1 << " 1x1 pivots so far." << std::endl;

				pivot(k, piv_info.r);
				pivot_vec(k, piv_info.r, A1, A1_idx);

#if 0
				std::cerr << "On iter " << k << " we chose diagonal pivot " << piv_info.r << std::endl;
				std::cerr << "Diagonal pivot value: " << D[k] << std::endl;
#endif
			}
#if 0	
			std::cerr << "Diagonal pivot value: " << D[k] << std::endl;
			for (int i = k; i < ncols; i++) {
				std::cerr << D[i] << " ";
			}
			std::cerr << endl;
#endif

			for (int j : A1_idx) {
				if (j < k) continue;
				D[j] = A1[j];
			}

			for (int j : A1_idx) {
				if (j <= k) continue;
				//if (std::abs(D[j]) < eps) continue;

				if (A1[k] == 0 || A1[k] != A1[k]) {
					std::cerr << "bad dkk" << std::endl;
					return;
				}

				// Compute new schur complement
				update_schur(-A1[j] / A1[k], j, k);
			}
			apply_dropping_rules(A1_idx, k);
			
			advance_list(k);
		}

#if 0
		std::cerr << "current permutation" << std::endl;
		for (int i = 0; i < ncols; i++) {
			std::cerr << p[i] << " ";
		}
		std::cerr << std::endl;
		std::cerr << "rest of ordering" << std::endl;
		for (int x : q) {
			std::cerr << pinv[x] << " ";
		}
		std::cerr << endl;

		
		std::cerr << "Current matrix: " << std::endl;
		for (int i = 0; i < ncols; i++) {
			for (int j = 0; j < ncols; j++) {
				std::cerr << coeff(p[i], p[j]) << " ";
			}
			std::cerr << std::endl;
		}
		std::cerr << endl;

		std::cerr << "Current factor: " << std::endl;
		for (int i = 0; i < ncols; i++) {
			for (int j = 0; j < ncols; j++) {
				std::cerr << L.coeff(i, j) << " ";
			}
			std::cerr << std::endl;
		}
		std::cerr << endl;
#endif

		//std::cerr << "iteration complete" << std::endl;
	}

	// assign number of non-zeros in L to L.nnz_count
	L.nnz_count = count;

	// apply permutations accumulated from pivoting
	sym_perm(p);
}

#endif

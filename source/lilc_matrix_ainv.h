#ifndef _LILC_MATRIX_AINV_H_
#define _LILC_MATRIX_AINV_H_

#include <default_pivoter.h>
#include <bkp_pivoter.h>
#include <set_unioner.h>

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

	double normA = 0;
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

		normA += norm(m_x[k]);
	}


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

	// Some structs for set union operations
	vector<bool> seen1(ncols, 0);

	// where the index value stored is > k
	vector<int> A1_idx, Ar_idx;
	vector<el_type> A1(ncols), Ar(ncols);

	// ------------ functions for ordering by number of non zeros -------------//
	// Tracks number of non-zeros per col so we can sort by sparsity
	vector<int> num_nz(ncols, 1);
	auto comp = [&](const int i, const int j) {
		if (num_nz[i] == num_nz[j]) return i < j;
		return num_nz[i] < num_nz[j];
	};

	set<int, decltype(comp)> q(comp);
	for (int i = 0; i < ncols; i++) {
		q.insert(i);
	}

	//------------------- convenience functions for pivoting ------------------//
	auto pivot = [&](int k, int r) {
		if (k == r) return;

		if (q.count(k)) {
			q.erase(k);
			num_nz[k] = static_cast<int>(L.m_idx[r].size());
			q.insert(k);
		}

		if (q.count(r)) {
			q.erase(r);
			num_nz[r] = static_cast<int>(L.m_idx[k].size());
			q.insert(r);
		}
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

	};

	auto advance_list = [&](int k) {
		// increase non-zero count of L
		count += static_cast<int>(L.m_x[k].size());

		// increment row_first array if needed
		for (int j : L.m_idx[k]) {
			while (!L.m_list[j].empty() && *L.m_list[j].begin() <= k) {
				L.m_list[j].erase(L.m_list[j].begin());
			}
		}
	};

	//------------------- main loop: factoring begins -------------------------//
	for (int k = 0; k < ncols; k++) {
		pivot(k, *q.begin());
		q.erase(q.begin());

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
		
		piv_info.size_two = false;
		// swap necessary rows and columns
		if (piv_info.size_two) {
			pivot(k+1, piv_info.r);
		} else {
			// figure out list of D[k]'s to compute and update
			pivoter->flush_col(A1, A1_idx, 0);

			// TODO: Apply dropping rules to pivot col
			//drop_tol_dense(A1, A1_idx, par.tol, piv_info.r);

			for (int j : A1_idx) {
				D[j] = A1[j];
			}

			// do swaps for pivots if necessary
			if (piv_info.r != k) {
				pivot(k, piv_info.r);
				std::swap(D[k], D[piv_info.r]);

#if 0
				std::cerr << "On iter " << k << " we chose diagonal pivot " << piv_info.r << std::endl;
				std::cerr << "Diagonal pivot value: " << D[k] << std::endl;
#endif

#if 0
				std::cerr << "current permutation" << std::endl;
				for (int i = 0; i < ncols; i++) {
					std::cerr << p[i] << " ";
				}
				std::cerr << std::endl;
#endif
			}

#if 0	
			std::cerr << "Diagonal pivot value: " << D[k] << std::endl;
			for (int i = k; i < ncols; i++) {
				std::cerr << D[i] << " ";
			}
			std::cerr << endl;
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
			// half-epsilon regularization
			if (std::abs(D[k]) < eps) {
#if 0
				std::cerr << "Diagonal pivot value: " << D[k] << std::endl;
				std::cerr << "Pivot value too low. Using half-epsilon regularization..." << std::endl;
#endif
				D[k] = eps * normA;
			}

			for (int j : A1_idx) {
				if (j == k) continue;
				if (std::abs(D[j]) < eps) continue;

				// Compute new schur complement
				col_wrapper<el_type> lk(L.m_x[k].data(), L.m_idx[k].data(), L.m_x[k].size()),
									 lj(L.m_x[j].data(), L.m_idx[j].data(), L.m_x[j].size());
				sparse_vec_add(1.0, lj, -D[j] / D[k], lk, work, curr_nnzs);

				// Apply dropping rules to schur complement
				drop_tol(work, curr_nnzs, par.tol, j);

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
			}
			
			advance_list(k);
		}
	}

	// assign number of non-zeros in L to L.nnz_count
	L.nnz_count = count;

	// apply permutations accumulated from pivoting
	sym_perm(p);
}

#endif

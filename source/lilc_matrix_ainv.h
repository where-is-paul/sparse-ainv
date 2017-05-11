#ifndef _LILC_MATRIX_AINV_H_
#define _LILC_MATRIX_AINV_H_

#include <default_pivoter.h>
#include <bkp_pivoter.h>
#include <wmn_pivoter.h>
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
	} else if (par.piv_type == pivot_type::WMN) {
		pivoter = new wmn_pivoter<el_type>(this, &L, &p, &pinv, par.beta);
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
		if (num_nz[i] == num_nz[j]) return i < j;
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
		if (k == r) return;
		q.erase(k);
		q.erase(r);

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

		L.m_list[k].swap(L.m_list[r]);
		L.m_x[k].swap(L.m_x[r]);
		L.m_idx[k].swap(L.m_idx[r]);
		// Because nnz indices are sorted, we can just change the last
		// elts indices since theyre guaranteed to be the diagonal
		assert(L.m_idx[k].back() == r);
		assert(L.m_idx[r].back() == k);
		L.m_idx[k][L.m_idx[k].size()-1] = k;
		L.m_idx[r][L.m_idx[r].size()-1] = r;

		// Swap values in A
		std::swap(p[k], p[r]);
		std::swap(pinv[p[k]], pinv[p[r]]);

		// Swap values in perm vector
		std::swap(perm[k], perm[r]);

		num_nz[k] = static_cast<int>(L.m_idx[k].size());
		num_nz[r] = static_cast<int>(L.m_idx[r].size());
		
		q.insert(k);
		q.insert(r);
	};

	auto advance_list = [&](int k) {
		// increment row_first array if needed
		for (int j : L.m_idx[k]) {
			L.m_list[j].erase(k);
		}

		// drop elements of col k
		q.erase(k);
		drop_tol(L.m_x[k], L.m_idx[k], par.tol, k, nullptr, true);
		
		// increase non-zero count of L
		count += static_cast<int>(L.m_x[k].size());
	};

	concurrent_unordered_map<int, concurrent_vector<int>> to_change;
	auto update_schur = [&](el_type coef, int j, int k) {
		col_wrapper<el_type> lk(L.m_x[k].data(), L.m_idx[k].data(), L.m_x[k].size()),
							 lj(L.m_x[j].data(), L.m_idx[j].data(), L.m_x[j].size());
	
		vector<int> curr_nnzs, extra;
		vector<el_type> work;
		sparse_vec_add(el_type(1.0), lj, coef, lk, work, curr_nnzs, &extra);

		L.m_x[j].swap(work);
		L.m_idx[j].swap(curr_nnzs);

		// Add to m_list if needed
		/*
		parallel_for(blocked_range<size_t>(0, extra.size(), 100), 
			[&](const blocked_range<size_t>& r) {
				for (int id = r.begin(); id != r.end(); id++) {
					int i = extra[id];
					to_change[i].push_back(j);
				}
			}
		);
		//*/
		///*
		for (int i : extra) {
			to_change[i].push_back(j);
		}
		//*/
	};

	auto update_schur_cleanup = [&](int j) {
		if (to_change.size() > 0) {
			parallel_for_each(to_change.begin(), to_change.end(),
				[&](const std::pair<int, concurrent_vector<int>>& p) {
					for (auto q : p.second) {
						L.m_list[p.first].insert(q);
					}
				}
			);
			to_change.clear();
		}

		q.erase(j);
		num_nz[j] = static_cast<int>(L.m_idx[j].size());
		q.insert(j);
	};

	auto apply_dropping_rules = [&](const vector<int>& pvts, int k) {
		double tol = 0;
		if (par.drop_type == drop_type::ABSOLUTE) {
			tol = par.tol / 4;
		} else if (par.drop_type == drop_type::RELATIVE) {
			// Figure out droptol
			double min_M_colsum = std::numeric_limits<double>::max();
			for (int i : pvts) {
				if (i <= k) continue;
				min_M_colsum = std::min(min_M_colsum, norm(L.m_x[i], 2.0));
			}
			tol = 0.1 * par.tol * min_M_colsum;
		}

		parallel_for_each(pvts.begin(), pvts.end(),
			[&](int i) {
				// Apply dropping rules to schur complement
				if (i <= k) return;

				vector<int> dropped;
				drop_tol(L.m_x[i], L.m_idx[i], tol, i, &dropped, true);
			
				for (int j : dropped) {
					to_change[j].push_back(i);
				}
			}
		);

		parallel_for_each(to_change.begin(), to_change.end(),
			[&](const std::pair<int, concurrent_vector<int>>& p) {
				for (auto q : p.second) {
					L.m_list[p.first].erase(q);			
				}
			}
		);
		to_change.clear();

		for (int i : pvts) {
			q.erase(i);
			num_nz[i] = static_cast<int>(L.m_x[i].size());
			q.insert(i);
		}
	
	};

	auto apply_dropping_rules_cleanup = [&](const vector<int>& pvts, int k) {
		for (int i : pvts) {
			if (i <= k) continue;

			// Fix ordering
			q.erase(i);
			num_nz[i] = static_cast<int>(L.m_x[i].size());
			q.insert(i);
		}
	};

	// Element dropping period
	int period = 1;//(int) (sqrt(ncols) + 1);
	int lastT = -1;
	int np0 = 0, np1 = 0, np2 = 0;
	std::map<int, bool> done;
	//------------------- main loop: factoring begins -------------------------//
	for (int k = 0; k < ncols; k++) {
		int nscr = int(0.5 + ceil(log10(num_nz[*(--q.end())])));

		int best = *q.begin();
		double bnorm = 0;
		auto it = q.begin();
		while (it != q.end() && nscr--) {
			double cnorm = norm(L.m_x[*it], 100.0); // approx of max norm, replace later
			if (cnorm > bnorm) {
				bnorm = cnorm;
				best = *it;
			}
		}
		pivot(k, best);

		int pcnt = (100*(k+1))/ncols;
		if (pcnt%10 == 0 && !done[pcnt]) {
			done[pcnt] = 1;
			std::cerr << k << " " << pcnt << " percent complete. ";
			int cnt = 0;
			for (int i = 0; i < ncols; i++) cnt += static_cast<int>(L.m_x[i].size());
			std::cerr << "nnz(M) = " << cnt << ", ";

			// Compute frobenius norm
			double fro = 0;
			for (int i = 0; i < ncols; i++) fro += pow(norm(L.m_x[i], 2.0), 2.0);
			std::cerr << "Frobenius norm: " << sqrt(fro) << ". Pivots: [" << np0 << " " << np1 << " " << np2 << "]. " << std::endl;
		}

		// need to pivot on ((A*M(:,p(dd)))'*M(:,p(dd:end)))';
		// or, M'(:, p(dd)) * (A M(:, p(dd:end))) = l_{pk}' * [A l_{pk} A l_{p{k+1}) ... A l_pn]
		// this is equivalent to pivoting on the D[j..n] vector
		pivot_struct piv_info = pivoter->find_pivot(k);
		
		// swap necessary rows and columns
		if (piv_info.size_two) {
			// figure out list of D[k]'s to compute and update
			pivoter->flush_col(A1, A1_idx, 0);
			pivoter->flush_col(Ar, Ar_idx, 1);

			if (piv_info.r != k+1) {
				pivot(k+1, piv_info.r);
				pivot_vec(k+1, piv_info.r, A1, A1_idx);
				pivot_vec(k+1, piv_info.r, Ar, Ar_idx);
			}
			np2++;

			// merge non-zero indices, compute accumlation depending on even or odd
			seen0.add_set(A1_idx);
			seen0.add_set(Ar_idx);
			seen0.flush(pvt_idx);

			// now i have the two columns.. merge their non-zero patterns
			el_type Di[2][2];
			Di[0][0] = A1[k]; 
			Di[0][1] = Di[1][0] = (A1[k+1] + Ar[k]) / 2;
			Di[1][1] = Ar[k+1];

			D[k] = Di[0][0];
			D.off_diagonal(k) = Di[0][1];
			D[k+1] = Di[1][1];

			// Z(:, [j0, j1]) -= Z(:, [k, k+1]) * [a b; c d]
			//[a * zk + c * z_k+1, b * zk + d * zk+1]
			
			el_type det = Di[0][0] * Di[1][1] - Di[0][1] * Di[1][0];
			parallel_for_each(pvt_idx.begin(), pvt_idx.end(),
				[&](int j) {
					if (j <= k+1) return;
					el_type DinvZ[2];
					DinvZ[0] =  A1[j] * Di[1][1] - Ar[j] * Di[0][1], 
					DinvZ[1] = -A1[j] * Di[1][0] + Ar[j] * Di[0][0];

					update_schur(-DinvZ[0] / det, j, k);
					update_schur(-DinvZ[1] / det, j, k+1);
				}
			);

			for (int j : pvt_idx) {
				update_schur_cleanup(j);
			}

			if (k/period > lastT) {
				apply_dropping_rules(pvt_idx, k+1);
				lastT = k/period;
			}

			advance_list(k);
			advance_list(k+1);
			// Extra iter accounted for
			k++;
		} else {
			// figure out list of D[k]'s to compute and update
			pivoter->flush_col(A1, A1_idx, 0);

			// do swaps for pivots if necessary
			if (piv_info.r != k) {
				np1++;
				pivot(k, piv_info.r);
				pivot_vec(k, piv_info.r, A1, A1_idx);
			} else {
				np0++;
			}


			for (int j : A1_idx) {
				D[j] = A1[j];
			}

			parallel_for_each(A1_idx.begin(), A1_idx.end(),
				[&](int j) {
					if (j <= k) return;

					// Compute new schur complement
					update_schur(-A1[j] / A1[k], j, k);
				}
			);

			for (int j : A1_idx) {
				update_schur_cleanup(j);
			}

			if (k/period > lastT) {
				apply_dropping_rules(A1_idx, k);
				lastT = k/period;
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

#ifndef _LILC_MATRIX_AINV_H_
#define _LILC_MATRIX_AINV_H_


using std::endl;
using std::cout;
using std::abs;

template <class el_type>
void lilc_matrix<el_type> :: ainv(lilc_matrix<el_type>& L, block_diag_matrix<el_type>& D, idx_vector_type& perm, params par) {
	//----------------- initialize temporary variables --------------------//
	const int ncols = n_cols(); //number of cols in A.

	elt_vector_type work(ncols, 0);
	idx_vector_type curr_nnzs;
	curr_nnzs.reserve(ncols); //reserves space for worse case (entire col is non-zero)

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
	//--------------- allocate memory for L and D ------------------//
	L.resize(ncols, ncols); //allocate a vector of size n for Llist as well
	D.resize(ncols);

	// Copy the nnz pattern of A into L.
	for (int k = 0; k < ncols; k++) {
		L.m_x[k].push_back(1.0);
		L.m_idx[k].push_back(k);
		L.m_list[k].insert(k);
	}

	vector<bool> seen(ncols, 0), seen2(ncols, 0);

	// where the index value stored is > k
	vector<int> to_update;
	//------------------- main loop: factoring begins -------------------------//
	for (int k = 0; k < ncols; k++) {
		// clear update queue
		to_update.clear();

		// figure out list of D[k]'s to compute and update
		for (int j : m_idx[k]) {
			// should union everything in list past row_first
			//for (int i = row_first[j]; i < L.m_list[j].size(); i++) {
			for (int idx : L.m_list[j]) {				
				if (!seen[idx]) {
					seen[idx] = true;
					to_update.push_back(idx);
				}
			}
		}
		sort(to_update.begin(), to_update.end());

		for (int j : to_update) {
			col_wrapper<el_type> ak(m_x[k].data(), m_idx[k].data(), m_x[k].size()),
								 lj(L.m_x[j].data(), L.m_idx[j].data(), L.m_x[j].size());
			D[j] = sparse_dot_prod(ak, lj);

			// half-epsilon regularization
			if (std::abs(D[j]) < eps) {
				D[j] = eps * normA;
			}
		}
		// half-epsilon regularization for k, since it might not be in to_update
		if (std::abs(D[k]) < eps) {
			D[k] = eps * normA;
		}

		for (int j : to_update) {
			if (j == k) continue;
			// Compute new schur complement
			col_wrapper<el_type> lk(L.m_x[k].data(), L.m_idx[k].data(), L.m_x[k].size()),
								 lj(L.m_x[j].data(), L.m_idx[j].data(), L.m_x[j].size());
			sparse_vec_add(1.0, lj, -D[j] / D[k], lk, work, curr_nnzs);

			// Apply dropping rules
			drop_tol(work, curr_nnzs, par.tol, j);

			// Add existing nnz into set
			for (int i : L.m_idx[j]) {
				seen2[i] = true;
			}

			L.m_x[j].swap(work);
			L.m_idx[j].swap(curr_nnzs);


			// Add to m_list if needed
			for (int i : L.m_idx[j]) {
				if (!seen2[i]) {
					//L.m_list[i].push_back(j);
					L.m_list[i].insert(j);
				}
				seen2[i] = false;
			}

			// Undo set additions
			for (int i : curr_nnzs) {
				if (seen2[i]) {
					L.m_list[i].erase(j);
				}
				seen2[i] = false;
			}
		}
		// increase non-zero count of L
		count += static_cast<int>(L.m_x[k].size());

		// reset seen array
		for (int j : to_update) {
			seen[j] = false;
		}

		// increment row_first array if needed
		for (int j : L.m_idx[k]) {
			while (!L.m_list[j].empty() && *L.m_list[j].begin() <= k) {
				L.m_list[j].erase(L.m_list[j].begin());
			}
		}
	}

	//assign number of non-zeros in L to L.nnz_count
	L.nnz_count = count;
}

#endif

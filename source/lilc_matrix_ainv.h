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
		L.m_list[k].push_back(k);
	}

	//------------------- main loop: factoring begins -------------------------//
	for (int k = 0; k < ncols; k++) {
		for (int j = k; j < ncols; j++) {
			col_wrapper<el_type> ak(m_x[k].data(), m_idx[k].data(), m_x[k].size()),
								 lj(L.m_x[j].data(), L.m_idx[j].data(), L.m_x[j].size());
			D[j] = sparse_dot_prod(ak, lj);
		}
		
		// half-epsilon regularization
		if (std::abs(D[k]) < eps * normA) {
			D[k] = eps * normA;
		}

		for (int j = k+1; j < ncols; j++) {
			if (D[j] == 0) continue;

			// Compute new schur complement
			col_wrapper<el_type> lk(L.m_x[k].data(), L.m_idx[k].data(), L.m_x[k].size()),
								 lj(L.m_x[j].data(), L.m_idx[j].data(), L.m_x[j].size());
			sparse_vec_add(1.0, lj, -D[j] / D[k], lk, work, curr_nnzs);

			// Apply dropping rules
			drop_tol(work, curr_nnzs, par.tol, j);

			L.m_x[j].swap(work);
			L.m_idx[j].swap(curr_nnzs);
		}

		count += static_cast<int>(L.m_x[k].size());
	}

	//assign number of non-zeros in L to L.nnz_count
	L.nnz_count = count;
}

#endif

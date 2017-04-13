#ifndef _LILC_MATRIX_SYM_METIS_H_
#define _LILC_MATRIX_SYM_METIS_H_

extern "C" {
#include <metis.h>
}

#include <vector>

// Change include file above to change support to floats, complex, complex doubles, etc.
template<class el_type> 
inline void lilc_matrix<el_type> :: sym_metis(std::vector<int>& perm) {
	int m = m_n_cols;
	std::vector<el_type> row_val;
	std::vector<idx_t> row_ind, col_ptr;
	std::vector<std::vector<idx_t>> adj(m);
	
	to_csc(row_val, row_ind, col_ptr);
	for (int i = 0; i < m; i++) {
		for (int j = col_ptr[i]; j < col_ptr[i+1]; j++) {
			if (i == row_ind[j]) continue;
			adj[i].push_back(row_ind[j]);
			adj[row_ind[j]].push_back(i);
		}
	}

	int last = 0;
	std::vector<idx_t> fullrow, fullcol;
	fullrow.push_back(last);
	for (int i = 0; i < m; i++) {
		last += adj[i].size();
		for (int j = 0; j < adj[i].size(); j++) {
			fullcol.push_back(adj[i][j]);
		}
		fullrow.push_back(last);
	}

	vector<idx_t> t_iperm, t_perm;
	t_perm.resize(m);
	t_iperm.resize(m);

	int output_error = METIS_NodeND(&m, fullrow.data(), fullcol.data(), NULL, NULL, t_perm.data(), t_iperm.data());
	perm.resize(m);
	for (int i = 0; i < t_perm.size(); i++) {
		perm[i] = t_perm[i];
	}
     
	if(output_error != METIS_OK) {
		std::cerr << "ERROR WHILE CALLING THE METIS PACKAGE" << std::endl; 
	}
	return;
}

#endif
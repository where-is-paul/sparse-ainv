//-*- mode: c++ -*-
#ifndef _LIL_MATRIX_SYM_EQUIL_H_
#define _LIL_MATRIX_SYM_EQUIL_H_

using std::abs;

template<class el_type>
void lilc_matrix<el_type> :: sym_equil() {

	//find termination points for loops with binary search later.
	int i, ncols = n_cols();
	// this is required since we do S[i] = max(S[i], ...)
	S.resize(ncols, 1.0);
	
	for (i = 0; i < ncols; i++) {
		int j = 0;
		if (!m_idx[i].empty()) {
			while (j < m_idx[i].size() && m_idx[i][j] < i) j++;
			
			if (j < m_idx[i].size() && m_idx[i][j] == i) {
				S[i] = sqrt(abs(m_x[i][j]));
			}
		}
		
		//assumes indices are ordered. since this procedure is run
		//before factorization pivots matrix, this is a fair assumption
		//for most matrix market matrices.
		j--;
		for (; j >= 0; j--) {
			S[i] = std::max(S[i], abs(m_x[i][j]));
		}
		
		//S[i] > 0 since its the square root of a +ve number
		if (S[i] > eps) { 
			for (int j = 0; j < m_idx[i].size(); j++) {
				m_x[i][j] /= S[i];
			}

			std::pair<idx_it, elt_it> elem_its;
			for (int j = 0; j < m_idx[i].size(); j++) {
				coeffRef(i, m_idx[i][j], elem_its);
				*elem_its.second /= S[i];
			}
		} 
	}
	
	for (i = 0; i < ncols; i++) {
		S[i]  = 1.0/S[i];
	}
}

#endif // _LIL_MATRIX_SYM_EQUIL_H_

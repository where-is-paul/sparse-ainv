#ifndef _DEFAULT_PIVOTER_H_
#define _DEFAULT_PIVOTER_H_

#include <pivot_strategy.h>

template<class el_type>
class default_pivoter : public pivot_strategy<el_type> {
public:
	default_pivoter() {}

	default_pivoter(lilc_matrix<el_type>* A_, lilc_matrix<el_type>* L_, const vector<int>* p_, const vector<int>* pinv_) 
		: pivot_strategy<el_type>(A_, L_, p_, pinv_) {
	}

	pivot_struct find_pivot(int col) { 
		this->update_col(this->A1, this->A1_idx, col);

		// regularization
		if (std::abs(this->A1[col]) <= this->m_bound) {
			this->A1[col] = (this->A1[col] >= 0 ? 1 : -1) * this->m_bound;
			return pivot_struct(false, col);
		}

		return pivot_struct(false, col);
	}
};

#endif
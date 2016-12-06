#ifndef _DEFAULT_PIVOTER_H_
#define _DEFAULT_PIVOTER_H_

#include <pivot_strategy.h>

template<class el_type>
class default_pivoter : public pivot_strategy<el_type> {
public:
	default_pivoter() {}

	default_pivoter(lilc_matrix<el_type>* A_, lilc_matrix<el_type>* L_, const vector<int>* p_) 
		: pivot_strategy<el_type>(A_, L_, p_) {
	}

	pivot_struct find_pivot(int col) { 
		update_col(A1, A1_idx, col);
		return pivot_struct(false, col);
	}
};

#endif
#ifndef _DEFAULT_PIVOTER_H_
#define _DEFAULT_PIVOTER_H_

#include <pivot_strategy.h>

class default_pivoter : public pivot_strategy {
public:
	default_pivoter() {}

	template<class el_type>
	pivot_struct find_pivot(const lilc_matrix<el_type>* mat, const vector<int>& p, int col) { 
		return pivot_struct(false, col);
	}
};

#endif
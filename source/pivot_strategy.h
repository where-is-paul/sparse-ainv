#ifndef _PIVOT_STRATEGY_H_
#define _PIVOT_STRATEGY_H_

using std::vector;

struct pivot_struct {
	bool size_two;
	int r;
	pivot_struct() {}
	pivot_struct(bool size_two, int r) : size_two(size_two), r(r) {}
};

class pivot_strategy {
public:
	pivot_strategy() {}

	template<class el_type>
	pivot_struct find_pivot(const lilc_matrix<el_type>* mat, const vector<int>& p, int col) { 
		return pivot_struct(false, col);
	}
};

#endif
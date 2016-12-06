#ifndef _PIVOT_STRATEGY_H_
#define _PIVOT_STRATEGY_H_

#include <set_unioner.h>

using std::vector;

struct pivot_struct {
	bool size_two;
	int r;
	pivot_struct() {}
	pivot_struct(bool size_two, int r) : size_two(size_two), r(r) {}
};

template<class el_type>
class pivot_strategy {
public:
	pivot_strategy() {}
	pivot_strategy(lilc_matrix<el_type>* A_, lilc_matrix<el_type>* L_, const vector<int>* p_) {
		A1.resize(A_->n_cols());
		Ar.resize(A_->n_cols());
		seen0.reset(A_->n_cols());
		A = A_;
		L = L_;
		p = p_;
	}

	pivot_struct find_pivot(int col) { 
		update_col(A1, A1_idx, col);
		return pivot_struct(false, col);
	}

	void flush_col(vector<el_type>& v, vector<int>& idx, int col = 0) {
		if (col == 0) {
			v.swap(A1);
			idx.swap(A1_idx);
		} else if (col == 1) {
			v.swap(Ar);
			idx.swap(Ar_idx);
		}
	}

protected:
	void update_col(vector<el_type>& v, vector<int>& idx, int k) {
		for (int j : A->m_idx[k]) {
			seen0.add_set(L->m_list[j]);
		}
		seen0.add_single(k);
		seen0.flush(idx);

		for (int j : idx) {
			col_wrapper<el_type> ak(A->m_x[k].data(), A->m_idx[k].data(), A->m_x[k].size()),
								 lj(L->m_x[j].data(), L->m_idx[j].data(), L->m_x[j].size());
			v[j] = sparse_dot_prod(ak, lj);
		}
	}

	vector<el_type> A1, Ar;
	vector<int> A1_idx, Ar_idx;
	lilc_matrix<el_type> *A, *L;
	const vector<int>* p;

	set_unioner seen0;
};

#endif
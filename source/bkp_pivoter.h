#ifndef _BKP_PIVOTER_H_
#define _BKP_PIVOTER_H_

#include <pivot_strategy.h>
#include <cmath>

class bkp_pivoter : public pivot_strategy {
public:
	bkp_pivoter() : m_beta(1.0) {
		m_alpha = (1 + std::sqrt(17.0)) / 8;
	}
	
	bkp_pivoter(double beta) : m_beta(beta) {
		// Compute this with bisection
		m_alpha = (1 + std::sqrt(17.0)) / 8;
	}

	// in: unpermuted index.
	// out: unpermuted index
	// PRECONDITION: ptr is sorted
	template<class el_type>
	pivot_struct find_pivot(const lilc_matrix<el_type>* mat, const vector<int>& p, int col) { 
		// r: unpermuted index
		double a11 = 0, w1 = 0;
		int r = 0;

		vector<el_type>& val = mat->m_x[p[col]];
		vector<el_type>& ptr = mat->m_idx[p[col]];
		for (int i = 0; i < ptr.size(); i++) {
			if (ptr[i] == p[col]) {
				a11 = std::abs(val[i]);
			}

			if (std::abs(val[i]) > w1) {
				w1 = std::abs(val[i]);
				r = ptr[i];
			}
		}

		if (a11 > m_alpha * m_beta * w1) {
			return pivot_struct(false, col);
		} else {
			// r_: unpermuted index
			double wr = 0, arr = 0;
			int r_ = 0;
			val = mat->m_x[p[r]];
			ptr = mat->m_idx[p[r]];
			for (int i = 0; i < ptr.size(); i++) {
				if (std::abs(val[i]) > m_beta * w1) {
					wr = std::abs(val[i]);
					r_ = ptr[i];
				}

				if (ptr[i] == p[r]) {
					arr = std::abs(val[i]);
				}
			}

			if (a11 * wr >= m_alpha * pow(m_beta * w1, 2.0)) {
				return pivot_struct(false, col);
			} else if (arr >= m_alpha * m_beta * wr) {
				return pivot_struct(false, r);
			} else {
				return pivot_struct(true, r);
			}
		}
	}

private:
	double m_beta, m_alpha;
};

#endif
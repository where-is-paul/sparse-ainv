#ifndef _BKP_PIVOTER_H_
#define _BKP_PIVOTER_H_

#include <pivot_strategy.h>
#include <cmath>

template<class el_type>
class bkp_pivoter : public pivot_strategy<el_type> {
public:
	bkp_pivoter(lilc_matrix<el_type>* A_, lilc_matrix<el_type>* L_, const vector<int>* p_, const vector<int>* pinv_) 
		: m_beta(1.0), pivot_strategy<el_type>(A_, L_, p_, pinv_) {
		m_alpha = (1 + std::sqrt(17.0)) / 8;
	}
	
	bkp_pivoter(lilc_matrix<el_type>* A_, lilc_matrix<el_type>* L_, const vector<int>* p_, const vector<int>* pinv_, double beta) 
		: pivot_strategy<el_type>(A_, L_, p_, pinv_) {
		m_beta = beta;
		m_alpha = (1 + std::sqrt(17.0)) / 8;
	}

	pivot_struct find_pivot(int col) {
		update_col(A1, A1_idx, col);

		double a11 = 0, w1 = 0;
		int r = -1;

		for (int j : A1_idx) {
			if (j < col) {
				continue;
			}
			el_type el = std::abs(A1[j]);
			if (j == col) {
				a11 = el;
			}

			if (el > w1) {
				w1 = el;
			}
		}

#if 0
		w1 = 1e9;
		for (int j : A1_idx) {
			el_type el = std::abs(A1[j]);
			if (j == col) {
				a11 = el;
				continue;
			}

			if (el < 1e-8) continue;
			if (el < w1) {
				w1 = el;
				r = j;
			}
		}

		std::cerr << "choosing " << r << " as pivot on col " << col << " with val " << w1 << std::endl;
		for (int j : A1_idx) {
			std::cerr << j << ":" << A1[j] << " ";
		}
		std::cerr << std::endl;

		update_col(Ar, Ar_idx, r);
		return pivot_struct(true, r);
#endif

		// regularization
		double bound = 64 * m_eps * m_reg;
		if (std::max(a11, w1) <= bound) {
			A1[col] = bound;
			return pivot_struct(false, col);
		} else if (a11 > m_alpha * m_beta * w1 + m_eps) {
			return pivot_struct(false, col);
		} else {
			double wr = 0, arr = 0;
			for (int j : A1_idx) {
				if (j < col) {
					continue;
				}
				el_type el = std::abs(A1[j]);
				if (el >= m_beta * w1 - m_eps) {
					r = j;
					break;
				}
			}

			update_col(Ar, Ar_idx, r);
			for (int j : Ar_idx) {
				if (j < col) continue;
				el_type el = std::abs(Ar[j]);
				if (j == r) {
					arr = el;
				} else if (el > wr) {
					wr = el;
				}
			}
		
			if (a11 * wr > m_alpha * pow(m_beta * w1, 2.0) + m_eps) {
				return pivot_struct(false, col);
			} else if (arr > m_alpha * m_beta * wr + m_eps) {
				A1.swap(Ar);
				A1_idx.swap(Ar_idx);
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
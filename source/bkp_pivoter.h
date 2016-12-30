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
		this->update_col(this->A1, this->A1_idx, col);

		double a11 = 0, w1 = 0;
		int r = -1;

		for (int j : this->A1_idx) {
#if 1
			if (j < col) {
				std::cerr << "this should be removed" << std::endl;
				continue;
			}
#endif
			el_type el = this->A1[j];
			if (j == col) {
				a11 = el;
			} else if (std::abs(el) > w1) {
				w1 = std::abs(el);
			}
		}

		// regularization
		if (std::max(std::abs(a11), w1) <= this->m_bound) {
			this->A1[col] = (a11 >= 0 ? 1 : -1) * this->m_bound;
			return pivot_struct(false, col);
		} else if (std::abs(a11) >= m_alpha * m_beta * w1 - this->m_eps) {
			return pivot_struct(false, col);
		} else {
			double wr = 0, arr = 0;
			// TODO: remove this second scan after ensuring bit-accuracy
			for (int j : this->A1_idx) {
#if 1
				if (j < col) {
					//TODO: Remove this
					std::cerr << "this should be removed" << std::endl;
					continue;
				}
#endif 

				if (j == col) {
					continue;
				}
				el_type el = std::abs(this->A1[j]);
				if (el >= m_beta * w1 - this->m_eps) {
					r = j;
					break;
				}
			}

			this->update_col(this->Ar, this->Ar_idx, r);
			for (int j : this->Ar_idx) {
				if (j < col) {
					std::cerr << "this shouldnt happen" << std::endl;
					continue;
				}
				el_type el = std::abs(this->Ar[j]);
				if (j == r) {
					arr = el;
				} else if (el > wr) {
					wr = el;
				}
			}
		
			if (std::abs(a11) * wr >= m_alpha * pow(m_beta * w1, 2.0) - this->m_eps) {
				return pivot_struct(false, col);
			} else if (arr >= m_alpha * m_beta * wr - this->m_eps) {
				this->A1.swap(this->Ar);
				this->A1_idx.swap(this->Ar_idx);
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
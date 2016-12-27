#ifndef _WMN_PIVOTER_H_
#define _WMN_PIVOTER_H_

#include <pivot_strategy.h>
#include <cmath>

template<class el_type>
class wmn_pivoter : public pivot_strategy<el_type> {
public:
	wmn_pivoter(lilc_matrix<el_type>* A_, lilc_matrix<el_type>* L_, const vector<int>* p_, const vector<int>* pinv_) 
		: m_beta(1.0), pivot_strategy<el_type>(A_, L_, p_, pinv_) {
		m_alpha = (1 + std::sqrt(17.0)) / 8;
	}
	
	wmn_pivoter(lilc_matrix<el_type>* A_, lilc_matrix<el_type>* L_, const vector<int>* p_, const vector<int>* pinv_, double beta) 
		: pivot_strategy<el_type>(A_, L_, p_, pinv_) {
		m_beta = beta;
		m_alpha = (1 + std::sqrt(17.0)) / 8;
	}

	pivot_struct find_pivot(int col) {
		this->update_col(this->A1, this->A1_idx, col);

		double a11 = 0, w1 = 0;
		int r = -1;

		for (int j : this->A1_idx) {
			el_type el = this->A1[j];
			if (j == col) {
				a11 = el;
			} else if (std::abs(el) > w1) {
				w1 = std::abs(el);
			}
		}

		// regularization
		double bound = this->m_eps * this->m_reg;
		if (std::max(std::abs(a11), w1) <= bound) {
			this->A1[col] = (a11 >= 0 ? 1 : -1) * bound;
			return pivot_struct(false, col);
		} else if (std::abs(a11) > m_alpha * m_beta * w1 + this->m_eps) {
			return pivot_struct(false, col);
		} else {
			double arr = 0;
			double ga = 0, gb = 0, gc = 0;
			for (int j : this->A1_idx) {
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
				if (j < col) continue;
				el_type el = std::abs(this->Ar[j]);
				if (j == r) {
					arr = el;
					break;
				}
			}
		
			ga = norm(this->A1, this->A1_idx, 2.0) / std::max(std::abs(a11), bound);
			gb = norm(this->Ar, this->Ar_idx, 2.0) / std::max(std::abs(arr), bound);
			seen0.add_set(this->A1_idx);
			seen0.add_set(this->Ar_idx);
			
			el_type d00, d11, d10;
			d00 = this->A1[col];
			d11 = this->Ar[r];
			d10 = this->A1[r];
			el_type det = d00 * d11 - d10 * d10;
			// [d11  -d10
			//  -d10  d00] / det
			
			for (int i : pvt_idx) {
				double x0 = ( this->A1[i] * d11 - this->Ar[i] * d10) / det,
					   x1 = (-this->A1[i] * d10 + this->Ar[i] * d00) / det;
				gc += x0 * x0 + x1 * x1;
			}
			gc = sqrt(gc);

			double eta = 0.1, gamma = 0.01;
			if (eta * ga < gc && gamma * ga < gb) {
				return pivot_struct(false, col);
			} else if (gc < eta * ga) {
				return pivot_struct(true, r);
			} else {
				this->A1.swap(this->Ar);
				this->A1_idx.swap(this->Ar_idx);
				return pivot_struct(false, r);
			}
		}
	}

private:
	double m_beta, m_alpha;
};

#endif
#ifndef _PIVOT_STRATEGY_H_
#define _PIVOT_STRATEGY_H_

#include <cfloat>
#include <set_unioner.h>

using std::vector;

struct pivot_struct {
  bool size_two;
  int r;
  pivot_struct() {}
  pivot_struct(bool size_two, int r) : size_two(size_two), r(r) {}
};

template <class el_type> class pivot_strategy {
public:
  pivot_strategy() {}
  pivot_strategy(lilc_matrix<el_type> *A_, lilc_matrix<el_type> *L_,
                 const vector<int> *p_, const vector<int> *pinv_) {
    A1.resize(A_->n_cols());
    Ar.resize(A_->n_cols());
    seen0.reset(A_->n_cols());

    A = A_;
    L = L_;
    p = p_;
    pinv = pinv_;
    m_eps = sqrt(std::numeric_limits<el_type>::epsilon());
  }

  void set_regularization(double reg) {
    m_reg = reg;
    m_bound = 64 * m_eps * m_reg;
  }

  virtual pivot_struct find_pivot(int col) {
    update_col(A1, A1_idx, col);
    return pivot_struct(false, col);
  }

  void flush_col(vector<el_type> &v, vector<int> &idx, int col = 0) {
    if (col == 0) {
      v.swap(A1);
      idx.swap(A1_idx);
      clean(A1, A1_idx);
    } else if (col == 1) {
      v.swap(Ar);
      idx.swap(Ar_idx);
      clean(Ar, Ar_idx);
    }
  }

protected:
  void update_col(vector<el_type> &v, vector<int> &idx, int k) {
    clean(v, idx);

//#define AINV_MODE 1
#ifdef AINV_MODE
    int pk = (*p)[k];
    for (int j : A->m_idx[pk]) {
      j = (*pinv)[j];
      seen0.add_set(L->m_list[j]);
    }
    v[k] = 0;
    seen0.add_single(k);
    seen0.flush(idx);

    for (int j : idx) {
      col_wrapper<el_type> ak(A->m_x[pk].data(), A->m_idx[pk].data(),
                              A->m_x[pk].size()),
          lj(L->m_x[j].data(), L->m_idx[j].data(), L->m_x[j].size());
      v[j] = sparse_dot_prod(ak, lj, pinv);
    }
#else
    // Compute a sparse-sparse matrix vector product
    static vector<el_type> Az_x, work;
    static vector<int> Az_idx;
    if (work.empty()) {
      work.resize(A->n_cols());
    }

    Az_x.clear();
    Az_idx.clear();

    // Want A * z_i --> sum z_ij Aj
    // need accumulation vector, and non-zero merger
    for (int i = 0; i < L->m_idx[k].size(); i++) {
      int pi = (*p)[L->m_idx[k][i]];
      el_type ci = L->m_x[k][i];
      seen0.add_set(A->m_idx[pi]);
      for (int j = 0; j < A->m_idx[pi].size(); j++) {
        int r = A->m_idx[pi][j];
        work[r] += ci * A->m_x[pi][j];
      }
    }
    seen0.flush(Az_idx);

    for (int i : Az_idx) {
      Az_x.push_back(work[i]);
      work[i] = 0;
    }

    for (int j : Az_idx) {
      j = (*pinv)[j];
      seen0.add_set(L->m_list[j]);
    }
    v[k] = 0;
    seen0.add_single(k);
    seen0.flush(idx);

    for (int j : idx) {
      col_wrapper<el_type> ak(Az_x.data(), Az_idx.data(), Az_x.size()),
          lj(L->m_x[j].data(), L->m_idx[j].data(), L->m_x[j].size());
      v[j] = sparse_dot_prod(ak, lj, pinv);
    }
#endif
  }

  void clean(vector<el_type> &v, vector<int> &idx) {
    for (int j : idx)
      v[j] = 0;
    idx.clear();
  }

  vector<el_type> A1, Ar;
  vector<int> A1_idx, Ar_idx, pvt_idx;
  lilc_matrix<el_type> *A, *L;
  const vector<int> *p, *pinv;

  set_unioner seen0;

  double m_eps, m_reg, m_bound;
};

#endif
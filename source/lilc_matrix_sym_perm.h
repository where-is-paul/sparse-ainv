//-*- mode: c++ -*-
#ifndef _LIL_MATRIX_SYM_PERM_H_
#define _LIL_MATRIX_SYM_PERM_H_

template <class el_type>
void lilc_matrix<el_type>::sym_perm(std::vector<int> &perm) {
  vector<idx_vector_type> m_idx_new(m_n_cols);
  vector<elt_vector_type> m_x_new(m_n_cols);

  int i, j, pi, pj;
  el_type px;
  vector<int> pinv(m_n_cols);
  for (i = 0; i < m_n_cols; i++) {
    pinv[perm[i]] = i;
    m_list[i].clear();
  }

  for (j = 0; j < m_n_cols;
       j++) { // no need to use function call n_cols() every iter
    pj = pinv[j];

    for (i = 0; i < (int)m_idx[j].size(); i++) {
      pi = pinv[m_idx[j][i]];
      px = m_x[j][i];

      m_idx_new[pi].push_back(pj);
      m_x_new[pi].push_back(px);
      m_list[pj].insert(pi);
      if (pi != pj) {
        m_list[pi].insert(pj);
      }
    }
  }

  m_idx.swap(m_idx_new);
  m_x.swap(m_x_new);
}

#endif // _LIL_MATRIX_SYM_PERM_H_
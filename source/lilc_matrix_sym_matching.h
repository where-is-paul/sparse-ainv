#ifndef _LILC_MATRIX_SYM_MC64_H_
#define _LILC_MATRIX_SYM_MC64_H_

#include <vector>

namespace {

typedef int ll;
const int N = 100, M = 100, INF = 0x3f3f3f3f;
// data structures and helper functions common to all flow routines
int par[N], first[N], nxt[2 * M], ep[2 * M], m;
ll flo[2 * M], cap[2 * M], cost[2 * M];
void init() {
  m = 0;
  memset(first, -1, sizeof first);
  memset(flo, 0, sizeof flo);
}
void add_edge(int a, int b, ll c = 1,
              ll p = 0) { // a,b - nodes, c, p - cap, price
  nxt[m] = first[ep[m] = a], first[ep[m]] = m, cap[m] = c, cost[m] = p, ++m;
  nxt[m] = first[ep[m] = b], first[ep[m]] = m, cap[m] = 0, cost[m] = -p, ++m;
}

// Minimum cost maximum flow, assuming there are no negative cost cycles
// USAGE: 1) init(); 2) add edges, 3) mcf_pot_init(numV) and 4)
// ll price=0,flow=0; while (ll df = mcf_update(s, t, price, numV)) flow += df;
//!
//! for sparse graphs, may help to change Dijkstra from O(N^2) to O(M*lgN)
//! code is provided in comments

bool vis[N];
ll pot[N], dist[N];
struct dsort {
  bool operator()(int a, int b) {
    return (dist[a] == dist[b]) ? (a < b) : (dist[a] < dist[b]);
  }
};
void mcf_pot_init(int s, int n) {
  memset(pot, 0, sizeof pot);
  // if all edge costs >= 0, we can safely return before the Bellman-Ford here
  for (int k = 1; k < n; ++k)
    for (int e = 0; e < m; ++e)
      if (cap[e])
        pot[ep[e ^ 1]] = std::min(pot[ep[e ^ 1]], pot[ep[e]] + cost[e]);
}
ll mcf_update(int s, int t, ll &price, int n) {
  memset(vis, 0, sizeof vis);
  memset(dist, INF, sizeof dist);
  dist[s] = 0;
  std::set<int, dsort> q;
  q.insert(s);
  while (q.size() > 0) {
    int u = *q.begin();
    q.erase(q.begin());
    if (vis[u] == 1)
      continue;
    vis[u] = 1;
    for (int e = first[u]; e != -1; e = nxt[e]) {
      int v = ep[e ^ 1];
      if (flo[e] < cap[e] && dist[v] > dist[u] + pot[u] - pot[v] + cost[e]) {
        q.erase(v);
        dist[v] = dist[u] + pot[u] - pot[v] + cost[e], par[v] = e;
        q.insert(v);
      }
    }
  }
  if (dist[t] >= INF)
    return 0;
  ll df = INF;
  for (int e, i = t; i != s; i = ep[e])
    e = par[i], df = std::min(df, cap[e] - flo[e]);
  for (int e, i = t; i != s; i = ep[e])
    e = par[i], flo[e] += df, flo[e ^ 1] -= df, price += df * cost[e];
  for (int i = 0; i < n; ++i)
    pot[i] = std::min(INF, dist[i] + pot[i]);
  return df;
}

std::vector<int> weighted_matching(const std::vector<double> &edge_weights,
                                   const std::vector<int> &edge_u,
                                   const std::vector<int> &edge_v) {
  return std::vector<int>();
}
}

// Change include file above to change support to floats, complex, complex
// doubles, etc.
template <>
inline std::vector<double> lilc_matrix<>::sym_matching(std::vector<int> &perm) {
  std::vector<double> row_val;
  std::vector<int> row_ind, col_ptr;

  to_csc(row_val, row_ind, col_ptr);

  // Assign permutation matrix
  perm.resize(m_n_cols);
  for (int i = 0; i < m_n_cols; i++) {
    perm[i] = i;
  }

  // Return scaling factors
  std::vector<double> res(m_n_cols, 1);
  return res;
}

#endif
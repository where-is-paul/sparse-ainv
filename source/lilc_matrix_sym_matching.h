#ifndef _LILC_MATRIX_SYM_MC64_H_
#define _LILC_MATRIX_SYM_MC64_H_

#include <algorithm>
#include <map>
#include <queue>
#include <set>
#include <vector>

class MinCostMaxFlow {
public:
  constexpr static double INF = 1e18;
  constexpr static double EPS = 1e-2;

  // There's extra allocation here because we use 2 * E as a proxy for
  // the number of vertices. TODO: Optimize this
  MinCostMaxFlow(const std::vector<double> &weights,
                 const std::vector<int> &leftEnd,
                 const std::vector<int> &rightEnd)
      : E(leftEnd.size()), par(2 * E), first(2 * E, -1), nxt(2 * E), ep(2 * E),
        m(0), flow(2 * E), capacity(2 * E), cost(2 * E), potential(2 * E) {
    for (size_t i = 0; i < weights.size(); i++) {
      add_edge(leftEnd[i], rightEnd[i], weights[i]);
    }
    assert(m == 2 * E);
  }

  int mcf_update(int s, int t, double &price) {
    dist.clear();
    inqueue.clear();

    dist[s] = 0.0;
    std::queue<int> q({s});
    while (!q.empty()) {
      int u = q.front();
      q.pop();
      inqueue.erase(u);

      // std::cerr << "Visiting " << u << std::endl;

      assert(dist.count(u));
      for (int e = first[u]; e != -1; e = nxt[e]) {
        // std::cerr << ep[e] << " " << ep[e^1] << " " << e << " " << dist[u] +
        // cost[e] << std::endl;
        int v = ep[e ^ 1];

        double dv = dist.count(v) ? dist[v] : INF;
        double du = dist[u]; // We know for sure u has dist
        if (flow[e] < capacity[e] &&
            // dv > du + cost[e]) {
            dv >= (1 + EPS) * (du + cost[e]) + EPS) {
          dist[v] = du + cost[e];
          par[v] = e;
          if (!inqueue.count(v)) {
            inqueue.insert(v);
            q.push(v);
          }
        }
      }
    }

    if (!dist.count(t)) {
      return 0;
    }

    std::cerr << "Found " << t << " in distance " << dist[t] << std::endl;

    double df = INF;
    int counter = 0;
    for (int e, i = t; i != s; i = ep[e]) {
      std::cerr << "Backtracking flow " << i << " " << ep[par[i]] << " " << s
                << " " << t << std::endl;
      if (counter++ > 10)
        throw std::runtime_error("test");
      e = par[i];
      df = std::min(df, capacity[e] - flow[e]);
    }
    for (int e, i = t; i != s; i = ep[e]) {
      // std::cerr << "Adding flow " << i << " " << par[i] << std::endl;
      e = par[i];
      flow[e] += df;
      flow[e ^ 1] -= df;
      price += df * cost[e];
    }
    // throw std::runtime_error("test");
    return df;
  }

  std::vector<std::pair<int, int>> run_flow(int source, int sink) {
    double price = 0;
    int nflow = 0;
    while (int df = mcf_update(source, sink, price)) {
      nflow += df;
      std::cerr << "Added delta flow " << df << " " << nflow << std::endl;
    }
    std::cerr << "Final flow value: " << nflow << std::endl;
    // Extract out edges
    std::vector<std::pair<int, int>> edges;
    for (size_t i = 0; i < flow.size(); i++) {
      if (flow[i] > 0) {
        edges.push_back({ep[i], ep[i ^ 1]});
      }
    }

    return edges;
  }

private:
  int m, N, E;
  std::vector<int> par, first, nxt, ep;
  std::vector<double> flow, capacity, cost, potential;
  std::map<int, double> dist;
  std::set<int> inqueue;

  void add_edge(int a, int b, double c) {
    // std::cerr << "Adding edge " << a << " " << b << " with cost " << c <<
    // std::endl;
    ep[m] = a;
    nxt[m] = first[a];
    first[a] = m;
    capacity[m] = 1;
    cost[m] = c;
    ++m;

    ep[m] = b;
    nxt[m] = first[b];
    first[b] = m;
    capacity[m] = 0;
    cost[m] = -c;
    ++m;
  }
};

// Change include file above to change support to floats, complex, complex
// doubles, etc.
template <>
inline std::vector<double> lilc_matrix<>::sym_matching(std::vector<int> &perm) {
  std::vector<double> row_val;
  std::vector<int> row_ind, col_ptr;

  to_csc(row_val, row_ind, col_ptr);

  std::vector<int> node_u, node_v;
  std::vector<double> weights;

  auto tform = [](double x) { return exp(x); };

  assert(col_ptr.size() == m_n_cols + 1);
  // Add a small diagonal
  std::vector<bool> has_diag(m_n_cols);
  for (int i = 0; i < col_ptr.size() - 1; i++) {
    for (int j = col_ptr[i]; j < col_ptr[i + 1]; j++) {
      // Apply transformations here to weight as needed
      weights.push_back(tform(row_val[node_u.size()]));
      node_u.push_back(i);
      node_v.push_back(m_n_cols + row_ind[j]);
      if (i == row_ind[j]) {
        has_diag[i] = true;
      }
    }
  }

  // Add some small diagonal elements
  for (int i = 0; i < m_n_cols; i++) {
    if (!has_diag[i]) {
      weights.push_back(tform(0));
      node_u.push_back(i);
      node_v.push_back(m_n_cols + i);
    }
  }

  // Now add a sink and a source that connects to the left and right
  int source = 2 * m_n_cols, sink = 2 * m_n_cols + 1;
  for (int i = 0; i < 2 * m_n_cols; i++) {
    weights.push_back(tform(0));
    if (i < m_n_cols) {
      node_u.push_back(source);
      node_v.push_back(i);
    } else {
      node_u.push_back(i);
      node_v.push_back(sink);
    }
  }

  MinCostMaxFlow flow(weights, node_u, node_v);
  auto flow_edges = flow.run_flow(source, sink);

  // Assign permutation matrix
  perm.resize(m_n_cols);
  for (const auto &e : flow_edges) {
    if (std::max(e.first, e.second) < 2 * m_n_cols) {
      int u = e.first % m_n_cols;
      int v = e.second % m_n_cols;
      perm[u] = v;
      std::cerr << "Matched edge " << u << " " << v << std::endl;
    }
  }

  // Return scaling factors (TODO)
  std::vector<double> res(m_n_cols, 1);
  return res;
}

#endif
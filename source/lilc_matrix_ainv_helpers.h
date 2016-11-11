#ifndef _LILC_MATRIX_AINV_HELPERS_H
#define _LILC_MATRIX_AINV_HELPERS_H

using std::abs;
using std::vector;

typedef vector<int>::iterator idx_it;

template<class el_type>
struct col_wrapper {
	el_type* val;
	int* ptr;
	size_t len;

	col_wrapper(el_type* val, int* ptr, size_t len) : val(val), ptr(ptr), len(len) {}
};

/*! \brief Computes the dot product of v and w. Only works when el_type is real right now.
	\param v the first vector whose dot product we wish to compute
	\param w the second vector whose dot product we wish to compute
*/
template <class el_type>
inline double dot_product(const vector<el_type>& v, const vector<el_type>& w) {
	double res = 0;
	for (int i = 0; i < v.size(); i++) {
		res += v[i]*w[i];
	}
	return res;
}

/*! \brief Computes the sum of two vectors: u = a*v + b*w. a and b are scalars, u, v, and w are vectors.
	\param u the storage vector for the result
*/
template <class el_type>
inline void vector_sum(double a, vector<el_type>& v, double b, vector<el_type>& w, vector<el_type>& u) {
	for (int i = 0; i < v.size(); i++) {
		u[i] = a*v[i] + b*w[i];
	}
}

/*! \brief Computes the maximum (in absolute value) element of v(curr_nnzs) and it's index.
	\param v the vector whose max element is to be computed.
	\param curr_nnzs a list of indices representing non-zero elements in v.
	\param r the index of the maximum element of v
	
	\return the max element of v.
*/
template <class el_type>
inline double max(vector<el_type>& v, vector<int>& curr_nnzs, int& r) { 
	double res = 0;
	for (idx_it it = curr_nnzs.begin(), end = curr_nnzs.end(); it != end; ++it) {
		if (abs(v[*it]) > res) {
			res = abs(v[*it]);
			r = *it;
		}
	}
	
	return res;
}

/*! \brief Computes the norm of v(curr_nnzs).
	\param v the vector whose norm is to be computed.
	\param curr_nnzs a list of indices representing non-zero elements in v.
	\param p The norm number.
	\return the norm of v.
*/
template <class el_type>
inline double norm(vector<el_type>& v, vector<int>& curr_nnzs, el_type p = 1) { 
	el_type res = 0;
	for (idx_it it = curr_nnzs.begin(), end = curr_nnzs.end(); it != end; ++it) {
		res += pow(abs(v[*it]), p);  
	}
	
	return pow(res, 1/p);
}

/*! \brief Computes the norm of v.
	\param v the vector whose norm is to be computed.
	\param p The norm number.
	\return the norm of v.
*/
template <class el_type>
inline double norm(vector<el_type>& v, el_type p = 1) { 
	el_type res = 0;
	for (int i = 0; i < v.size(); i++) {
		res += pow(abs(v[i]), p);
	}
	return pow(res, 1/p);
}


//-------------Dropping rules-------------//
namespace {
/*! \brief Functor for comparing elements by value (in decreasing order) instead of by index.
	\param v the vector that contains the values being compared.
*/
template <class el_type>
struct by_value {
	const vector<el_type>& v; 
	by_value(const vector<el_type>& vec) : v(vec) {}
	inline bool operator()(int const &a, int const &b) const { 
		// Not needed if we're using sort. If using this comparator
		// in a set, then uncomment the line below.
		//if (abs(v[a]) == abs(v[b])) return a < b;
		return abs(v[a]) > abs(v[b]);
	}
};

/*! \brief Functor for determining if a variable is below the tolerance given.
    \param v the vector that contains the values being checked.
    \param eps the tolerance given.
*/
template <class el_type>
struct by_tolerance {
  const vector<el_type>& v; 
  double eps;
	by_tolerance(const vector<el_type>& vec, const double& eps) : v(vec), eps(eps) {}
	inline bool operator()(int const &i) const { 
		return abs(v[i]) < eps;
	}
};
}

/*! \brief Performs the dual-dropping criteria outlined in Li & Saad (2005).
	\param v the vector that whose elements will be selectively dropped.
	\param curr_nnzs the non-zeros in the vector v.
	\param tol a parameter to control agressiveness of dropping. Elements less than tol*norm(v) are dropped.
*/
template <class el_type>
inline void drop_tol(vector<el_type>& vals, vector<int>& curr_nnzs, const double& tol, int keep) { 
	static vector<el_type> work;
	static vector<int> nnzs;

	//determine dropping tolerance. all elements with value less than tolerance = tol * norm(v) is dropped.
	el_type tolerance = tol * norm<el_type>(vals);

	work.clear();
	nnzs.clear();
	for (int i = 0; i < curr_nnzs.size(); i++) {
		if (curr_nnzs[i] != keep && std::abs(vals[i]) < tolerance) continue;
		work.push_back(vals[i]);
		nnzs.push_back(curr_nnzs[i]);
	}

	vals.swap(work);
	curr_nnzs.swap(nnzs);
}

// PRECONDITION: a and b have sorted indices
template<class el_type>
el_type sparse_dot_prod(const col_wrapper<el_type>& a, const col_wrapper<el_type>& b) {
	el_type res = 0;

#ifdef USE_BLAS
#else
	int i = 0, j = 0;
	while (i < a.len && j < b.len) {
		if (a.ptr[i] == b.ptr[j]) {
			res += a.val[i] * b.val[j];
			i++;
			j++;
		} else if (a.ptr[i] < b.ptr[j]) {
			i++;
		} else {
			j++;
		}
	}
#endif
	return res;
}

// PRECONDITION: a and b have sorted indices
template<class el_type>
int sparse_vec_add(el_type c0, const col_wrapper<el_type>& a, el_type c1, const col_wrapper<el_type>& b, vector<el_type>& val, vector<int>& ptr) {
	val.clear();
	ptr.clear();

#ifdef USE_BLAS
#else
	int i = 0, j = 0;
	el_type res;
	int res_idx;
	while (i < a.len || j < b.len) {
		if (i < a.len && j < b.len) {
			if (a.ptr[i] == b.ptr[j]) {
				res = c0 * a.val[i] + c1 * b.val[j];
				res_idx = a.ptr[i];
				i++;
				j++;
			} else if (a.ptr[i] < b.ptr[j]) {
				res = c0 * a.val[i];
				res_idx = a.ptr[i];
				i++;
			} else {
				res = c1 * b.val[j];
				res_idx = b.ptr[j];
				j++;
			}
		} else if (i < a.len) {
			res = c0 * a.val[i];
			res_idx = a.ptr[i];
			i++;
		} else {
			res = c1 * b.val[j];
			res_idx = b.ptr[j];
			j++;
		}

		val.push_back(res);
		ptr.push_back(res_idx);
	}
#endif
	return static_cast<int>(val.size());
}

#endif
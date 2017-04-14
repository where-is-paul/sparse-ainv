#ifndef _SET_UNIONER_H_
#define _SET_UNIONER_H_

#include <tbb/parallel_for_each.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>

#include <vector>

using std::vector;

using namespace tbb;

struct set_unioner {
	concurrent_vector<bool> in_set;
	concurrent_vector<int> indices;

	void reset(int sz) {
		in_set.clear();
		in_set.resize(sz, 0);
		indices.clear();
	}

	// PRECONDITION: s contains no duplicates
	template<class Container>
	void add_set(const Container& s) {
		parallel_for_each(s.begin(), s.end(),
			[&](int x) {
				if (!in_set[x]) {
					indices.push_back(x);
					in_set[x] = true;
				}
			}
		);
	}

	void add_single(const int& x) {
		if (!in_set[x]) {
			indices.push_back(x);
			in_set[x] = true;
		}
	}

	void flush(vector<int>& res) {
		parallel_for_each(indices.begin(), indices.end(),
			[&](int x) {
				in_set[x] = false;
			}
		);

		res.clear();
		for (int x : indices) {
			res.push_back(x);
		}
		indices.clear();
	}
};

#endif
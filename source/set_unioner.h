#ifndef _SET_UNIONER_H_
#define _SET_UNIONER_H_

#include <vector>

using std::vector;

struct set_unioner {
	vector<bool> in_set;
	vector<int> indices;

	void reset(int sz) {
		in_set.clear();
		in_set.resize(sz, 0);
		indices.clear();
	}

	template<class Container>
	void add_set(const Container& s) {
		for (int x : s) {
			if (!in_set[x]) {
				indices.push_back(x);
				in_set[x] = true;
			}
		}
	}

	void add_single(const int& x) {
		if (!in_set[x]) {
			indices.push_back(x);
			in_set[x] = true;
		}
	}

	void flush(vector<int>& res) {
		for (int x : indices) {
			in_set[x] = false;
		}
		res.swap(indices);
		indices.clear();
	}
};

#endif
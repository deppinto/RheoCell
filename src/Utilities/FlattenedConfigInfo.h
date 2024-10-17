#ifndef SRC_UTILITIES_FLATTENEDCONFIGINFO_H_
#define SRC_UTILITIES_FLATTENEDCONFIGINFO_H_

#include "../defs.h"

#include <vector>
#include <memory>

class BaseField;

struct FlattenedVectorArray {
	int rows();
	int cols();
	size_t size();

	void set(int p_idx, std::vector<number> &v);

	std::vector<number> data;
};

struct FlattenedConfigInfo {
	void update(long long int step, const std::vector<BaseField *> &fields);

	FlattenedVectorArray positions;
	std::vector<int> types;

	long long int last_updated = -1;
};

#endif /* SRC_UTILITIES_FLATTENEDCONFIGINFO_H_ */

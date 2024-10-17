#include "FlattenedConfigInfo.h"

#include "../Fields/BaseField.h"

int FlattenedVectorArray::rows() {
	return size() / cols();
}

int FlattenedVectorArray::cols() {
	return 3;
}

size_t FlattenedVectorArray::size() {
	return data.size();
}

void FlattenedVectorArray::set(int p_idx, std::vector<number> &v) {
	int base_idx = p_idx * cols();
	data[base_idx] = v[0];
	data[base_idx + 1] = v[1];
	data[base_idx + 2] = v[2];
}

void FlattenedConfigInfo::update(long long int step, const std::vector<BaseField *> &fields) {
	if(step == last_updated) {
		return;
	}

	last_updated = step;
	int N = fields.size();

	positions.data.resize(N * positions.cols());
	types.resize(N);

	for(int i = 0; i < N; i++) {
		BaseField *p = fields[i];

		positions.set(i, p->CoM);
		types[i] = p->type;
	}
}

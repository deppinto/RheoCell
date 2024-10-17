#include "BaseForce.h"

#include "../Boxes/BaseBox.h"

#include <vector>
#include <string>
#include <tuple>

BaseForce::BaseForce() {
	p_ptr = P_VIRTUAL;
}

BaseForce::~BaseForce() {

}

std::tuple<std::vector<int>, std::string> BaseForce::init(input_file &inp) {
	getInputString(&inp, "group_name", group_name, 0);
	getInputString(&inp, "id", id, 0);
	getInputString(&inp, "type", type, 1);

	return std::make_tuple(std::vector<int>(), "BaseForce");
}

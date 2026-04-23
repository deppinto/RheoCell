#include "BaseField.h"
#include "../Boxes/BaseBox.h"

BaseField::BaseField() :
				index(-1), LsubX(30), LsubY(30), border(4), area(0.), sumF(0.) {
	CoM = std::vector<number> {0., 0.};
	CoM_old = std::vector<number> {0., 0.};
}
				

void BaseField::copy_from(const BaseField &p) {
	index = p.index;
	CoM = p.CoM;
	CoM_old = p.CoM;
	for(int i=0; i<subSize; i++){
		velocityX[i] = p.velocityX[i];
		velocityY[i] = p.velocityY[i];
	}
}

void BaseField::get_interaction_values(int R) {

}

BaseField::~BaseField() {

}

bool BaseField::add_ext_force(BaseForce *f) {
	ext_forces.push_back(f);
	return true;
}

void BaseField::init() {
}

void BaseField::check() {
	assert(index >= 0);
}

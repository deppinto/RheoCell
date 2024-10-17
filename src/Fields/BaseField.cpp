#include "BaseField.h"
#include "../Boxes/BaseBox.h"

BaseField::BaseField() :
				index(-1), LsubX(25), LsubY(25), area(0.), sumF(0.) {
	ext_potential = (number) 0.;
	force = std::vector<number> {0., 0.};
	CoM = std::vector<number> {0., 0.};
}
				

void BaseField::copy_from(const BaseField &p) {
	index = p.index;
	CoM = p.CoM;
	velocityX = p.velocityX;
	velocityY = p.velocityY;
	force = p.force;
	ext_potential = p.ext_potential;
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
	force = std::vector<number> {0., 0.};
	check();
}

void BaseField::check() {
	assert(index >= 0);
}

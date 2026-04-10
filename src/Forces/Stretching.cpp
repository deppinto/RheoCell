#include "Stretching.h"
#include "../Utilities/RCexception.h"
#include "../Fields/BaseField.h"

Stretching::Stretching() :
				BaseForce() {
	//lambda_wall = 1.;
	//lambda = 2.;
	//kappa_wall = 0.;
	shear_rate = 0.;
	shear_rate_active = 0.;
	generate_inside = false;
}

std::tuple<std::vector<int>, std::string> Stretching::init(input_file &inp) {
	BaseForce::init(inp);

	std::string field_string;
	getInputString(&inp, "field", field_string, 1);

	//getInputNumber(&inp, "lambda_wall", &lambda_wall, 0);
	//getInputNumber(&inp, "kappa_wall", &kappa_wall, 1);
	//getInputNumber(&inp, "lambda", &lambda, 1);
	getInputNumber(&inp, "shear_rate", &shear_rate_active, 1);

	getInputBool(&inp, "generate_inside", &generate_inside, 0);

	auto field_ids = Utils::get_fields_from_string(CONFIG_INFO->fields(), field_string, "Stretching");
	std::string description = Utils::sformat("Stretching (Rate=%g)", shear_rate_active);

	return std::make_tuple(field_ids, description);
}


void Stretching::apply_changes_after_equilibration(){
	shear_rate = shear_rate_active;
}

number Stretching::free_energy(number phi, number walls, number laplacianSquare) {

	return 0.;
}

number Stretching::free_energy_density(number phi, number walls, number laplacianSquare) {

	return 0.;
}

number Stretching::potential(int k, number walls) {
	return 0.;
}

std::vector<number> Stretching::velocity_profile(int k, int lx, int ly) {
	//int x = (k-(int(k/lx)*lx));
	int y = (k/lx);
	int x = k - y * lx;
	if(x==0) return std::vector<number> { - shear_rate , 0.};
	else if(x==lx-1) return std::vector<number> { shear_rate , 0.};
	else return std::vector<number> {0., 0.};
}

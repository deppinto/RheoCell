#include "KolmogorovFlowChannel.h"
#include "../Utilities/RCexception.h"
#include "../Fields/BaseField.h"

KolmogorovFlowChannel::KolmogorovFlowChannel() :
				BaseForce() {
	lambda_wall = 0.;
	lambda = 2.;
	kappa_wall = 0.;
	shear_rate = 0.;
	shear_rate_active = 0.;
	generate_inside = false;
}

std::tuple<std::vector<int>, std::string> KolmogorovFlowChannel::init(input_file &inp) {
	BaseForce::init(inp);

	std::string field_string;
	getInputString(&inp, "field", field_string, 1);

	getInputNumber(&inp, "lambda_wall", &lambda_wall, 0);
	getInputNumber(&inp, "kappa_wall", &kappa_wall, 1);
	getInputNumber(&inp, "lambda", &lambda, 1);
	getInputNumber(&inp, "shear_rate", &shear_rate_active, 1);

	//lambda_wall += 2 * lambda;

	getInputBool(&inp, "generate_inside", &generate_inside, 0);

	auto field_ids = Utils::get_fields_from_string(CONFIG_INFO->fields(), field_string, "KolmogorovFlowChannel");
	std::string description = Utils::sformat("KolmogorovFlowChannel (kappa_wall=%g)", kappa_wall);

	return std::make_tuple(field_ids, description);
}


void KolmogorovFlowChannel::apply_changes_after_equilibration(){
	shear_rate = shear_rate_active;
}

number KolmogorovFlowChannel::free_energy(number phi, number walls, number laplacianSquare) {

	return 4.0*kappa_wall/lambda*phi*walls*walls;
}

number KolmogorovFlowChannel::free_energy_density(number phi, number walls, number laplacianSquare) {

	return 2.0*kappa_wall/lambda*phi*phi*walls*walls;
}

number KolmogorovFlowChannel::potential(int k, number walls) {
	return walls;
}

std::vector<number> KolmogorovFlowChannel::velocity_profile(int k, int lx, int ly) {
	//int x = (k-(int(k/lx)*lx));
	int y = (k/lx);
	if(y > lambda_wall && y < ly - lambda_wall) return std::vector<number> { shear_rate * sin(2*PI*y/(ly/2)), 0.};
	else return std::vector<number> {0., 0.};
}

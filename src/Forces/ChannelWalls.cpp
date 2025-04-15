#include "ChannelWalls.h"
#include "../Utilities/RCexception.h"
#include "../Fields/BaseField.h"

ChannelWalls::ChannelWalls() :
				BaseForce() {
	//lambda_wall = 1.;
	lambda = 2.;
	kappa_wall = 0.;
	omega_wall = 0.;
	generate_inside = false;
}

std::tuple<std::vector<int>, std::string> ChannelWalls::init(input_file &inp) {
	BaseForce::init(inp);

	std::string field_string;
	getInputString(&inp, "field", field_string, 1);

	//getInputNumber(&inp, "lambda_wall", &lambda_wall, 0);
	getInputNumber(&inp, "kappa_wall", &kappa_wall, 0);
	getInputNumber(&inp, "omega_wall", &omega_wall, 0);
	getInputNumber(&inp, "lambda", &lambda, 0);

	getInputBool(&inp, "generate_inside", &generate_inside, 0);

	auto field_ids = Utils::get_fields_from_string(CONFIG_INFO->fields(), field_string, "ChannelWalls");
	std::string description = Utils::sformat("ChannelWalls (kappa_wall=%g, omega_wall=%g)", kappa_wall, omega_wall);

	return std::make_tuple(field_ids, description);
}

number ChannelWalls::free_energy(number phi, number walls, number laplacianSquare) {

	return 2.0*kappa_wall/lambda*phi*walls*walls - 4.0*omega_wall*lambda*phi; //((laplacianSquare)/sqrt(1+(laplacianSquare)*(laplacianSquare)));
}

number ChannelWalls::potential(int k, number walls) {
	return walls;
}

#include "HardWall.h"
#include "../Fields/BaseField.h"

HardWall::HardWall() : BaseForce(), position(-1.), sigma(1.) {
}


HardWall::~HardWall() {

}

/*std::tuple<std::vector<int>, std::string> HardWall::init(input_file &inp) {
	BaseForce::init(inp);

	std::string particles_string;
	getInputString(&inp, "particle", particles_string, 1);

	getInputNumber(&inp, "position", &_position, 1);
	getInputNumber(&inp, "sigma", &_sigma, 0);

	int tmpi;
	number tmpf[3];
	std::string strdir;
	getInputString(&inp, "dir", strdir, 1);
	tmpi = sscanf(strdir.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
	if(tmpi != 3) {
		throw ("Could not parse dir %s in external forces file. Aborting", strdir.c_str());
	}
	direction = vector<number> ((number) tmpf[0], (number) tmpf[1]);
	direction.normalize();

	auto particle_ids = Utils::get_particles_from_string(CONFIG_INFO->particles(), particles_string, "HardWall");
	std::string description = Utils::sformat("HardWall (position=%g, dir=%g,%g,%g, sigma=%g)", _position, _direction.x, _direction.y, _direction.z, _sigma);

	return std::make_tuple(particle_ids, description);
}*/

std::vector<number> HardWall::value(llint step, std::vector<number> &pos) {
	throw ("HardWall cannot be used");
}

number HardWall::potential(llint step, std::vector<number> &pos) {
	number distance = (direction[0] * pos[0] + direction[1] * pos[1] + position) / sigma; // distance from the plane in units of sigma
	if(distance < 1.) return 10e8;
	return 0.;
}

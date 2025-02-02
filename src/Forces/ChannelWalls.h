#ifndef CHANNELWALLS_H_
#define CHANNELWALLS_H_

#include "BaseForce.h"

class ChannelWalls: public BaseForce {
private:
	bool generate_inside;
	//number lambda_wall;
	number kappa_wall;
	number omega_wall;
	//std::vector<double> walls;
	number lambda;

public:
	ChannelWalls();
	virtual ~ChannelWalls() {
	}

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	virtual number free_energy(number phi, number walls, number laplacianSquare);
	virtual number potential(int k, number walls);
};

#endif // CHANNELWALLS_H

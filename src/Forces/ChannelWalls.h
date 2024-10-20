#ifndef CHANNELWALLS_H_
#define CHANNELWALLS_H_

#include "BaseForce.h"

/**
 * @brief External force field that confines cells into a channel. It can feature a short-range attraction.
 *
 * Example section in the external forces file:

 \n
 {\n
 particle = -1          # acts on all particles\n
 type = LJ_plane\n
 position = 7.          # in simulation unit lengths\n
 stiff = 50.            # quite stiff. Good for MC, not for MD\n
 sigma = 1.5            # diameter of the wall
 generate_inside = true # useful when generating the initial configuration
 }\n\n

 * @verbatim
 position = <float> (defines the position of the plane along the direction identified by the plane normal.)
 particle = <int> (index of the particle on which the force shall be applied. If -1, the force will be exerted on all the particles.)
 [stiff = <float> (stiffness of the repulsion. Defaults to 1.)]
 [lambda = <float> ("Diameter" of the wall. It effectively rescales the distance between particle and wall. Defaults to 1.)]
 [only_repulsive = <bool> (If true, the interactio between particle and wall gets cut-off at the minimum, resulting in a purely-repulsive potential. Defaults to false.)]
 [generate_inside = <bool> (If true the wall-particle interaction may not diverge, even for negative distances. Useful when generating the starting configuration. Defaults to false)]
 @endverbatim
 */

class ChannelWalls: public BaseForce {
private:
	bool generate_inside;
	//number lambda_wall;
	number kappa_wall;
	//std::vector<double> walls;
	number lambda;

public:
	ChannelWalls();
	virtual ~ChannelWalls() {
	}

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	virtual number free_energy(int k, number phi, number walls);
	virtual number potential(int k, number walls);
};

#endif // CHANNELWALLS_H

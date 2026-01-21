#ifndef SHEARFLOWCHANNEL_H_
#define SHEARFLOWCHANNEL_H_

#include "BaseForce.h"

class ShearFlowChannel: public BaseForce {
private:
	bool generate_inside;
	//number lambda_wall;
	number kappa_wall;
	number shear_rate;
	number shear_rate_active;
	//std::vector<double> walls;
	number lambda;

public:
	ShearFlowChannel();
	virtual ~ShearFlowChannel() {}

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	virtual number free_energy(number phi, number walls, number laplacianSquare);
	virtual number potential(int k, number walls);
	virtual std::vector<number> velocity_profile(int k, int lx, int ly);
	virtual void apply_changes_after_equilibration();
};

#endif // SHEARFLOWCHANNEL_H

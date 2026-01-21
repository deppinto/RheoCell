#ifndef POISEUILLEFLOWCHANNEL_H_
#define POISEUILLEFLOWCHANNEL_H_

#include "BaseForce.h"

class PoiseuilleFlowChannel: public BaseForce {
private:
	bool generate_inside;
	//number lambda_wall;
	number kappa_wall;
	number shear_rate;
	number shear_rate_active;
	//std::vector<double> walls;
	number lambda;

public:
	PoiseuilleFlowChannel();
	virtual ~PoiseuilleFlowChannel() {}

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	virtual number free_energy(number phi, number walls, number laplacianSquare);
	virtual number potential(int k, number walls);
	virtual std::vector<number> velocity_profile(int k, int lx, int ly);
	virtual void apply_changes_after_equilibration();
};

#endif // POISEUILLEFLOWCHANNEL_H

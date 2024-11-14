#include "BoxFactory.h"
#include "SquareBox.h"
#include "OrthogonalBox.h"
#include "LeesEdwardsSquareBox.h"
#include "Channel.h"
#include "SquareWalls.h"

#include "../Utilities/RCexception.h"

BoxPtr BoxFactory::make_box(input_file &inp) {
	// the default box is the square one
	std::string box_type("square");
        getInputString(&inp, "box", box_type, 0);
	bool lees_edwards = false;
	getInputBool(&inp, "lees_edwards", &lees_edwards, 0); 	

	std::string type_str("none");
	getInputString(&inp, "type", type_str, 0);


	if(box_type.compare("square") == 0) {
		if(lees_edwards) return std::make_shared<LeesEdwardsSquareBox>();
			
		if(type_str.compare("channel_walls") == 0) return std::make_shared<SquareWalls>();
		else return std::make_shared<SquareBox>();
	}
	if(box_type.compare("orthogonal") == 0) {
		if(lees_edwards) throw RCexception("Lees-Edwards boundary conditions are not compatible with Orthogonal Box!");

		if(type_str.compare("channel_walls") == 0) return std::make_shared<Channel>();
		else return std::make_shared<OrthogonalBox>();
	}
	else throw RCexception("Unsupported box");
}

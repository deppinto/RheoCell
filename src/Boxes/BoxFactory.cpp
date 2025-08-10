#include "BoxFactory.h"
#include "SquareBox.h"
#include "OrthogonalBox.h"
#include "LeesEdwardsBox.h"
#include "Channel.h"
#include "SquareWalls.h"
#include "CircleWalls.h"
#include "TriangleWalls.h"
#include "MovingWalls.h"

#include "../Utilities/RCexception.h"

BoxPtr BoxFactory::make_box(input_file &inp) {
	// the default box is the square one
	std::string box_type("square");
        getInputString(&inp, "box", box_type, 0);

	std::string type_str("none");
	getInputString(&inp, "type", type_str, 0);

	if(box_type.compare("square") == 0) {
		if(type_str.compare("channel_walls") == 0) return std::make_shared<SquareWalls>();
		else if(type_str.compare("circle_walls") == 0) return std::make_shared<CircleWalls>();
		else if(type_str.compare("triangle_walls") == 0) return std::make_shared<TriangleWalls>();
		else return std::make_shared<SquareBox>();
	}
	else if(box_type.compare("orthogonal") == 0) {
		if(type_str.compare("channel_walls") == 0) return std::make_shared<Channel>();
		else if(type_str.compare("moving_walls") == 0) return std::make_shared<MovingWalls>();
		else if(type_str.compare("shear_flow_channel") == 0) return std::make_shared<Channel>();
		else if(type_str.compare("pois_flow_channel") == 0) return std::make_shared<Channel>();
		else return std::make_shared<OrthogonalBox>();
	}
	else if(box_type.compare("lees_edwards") == 0) {
		if(type_str.compare("channel_walls") == 0)throw RCexception("Channel walls are not compatible with Lees-Edwards boundary conditions!"); 
		else return std::make_shared<LeesEdwardsBox>();
	}
	else throw RCexception("Unsupported box");
}

#include "BoxFactory.h"
#include "SquareBox.h"
#include "LeesEdwardsSquareBox.h"

BoxPtr BoxFactory::make_box(input_file &inp) {
	// the default box is the square one
	std::string box_type("square");
        getInputString(&inp, "box", box_type, 0);
	bool lees_edwards = false;
	getInputBool(&inp, "lees_edwards", &lees_edwards, 0); 	

	if(box_type.compare("square") == 0) {
		if(!lees_edwards) return std::make_shared<SquareBox>();
		else return std::make_shared<LeesEdwardsSquareBox>();
	}
	else throw("Unsupported box");
}
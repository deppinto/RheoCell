#include "LeesEdwardsSquareBox.h"
#include "../Utilities/ConfigInfo.h"

LeesEdwardsSquareBox::LeesEdwardsSquareBox() :
				SquareBox() {

}


LeesEdwardsSquareBox::~LeesEdwardsSquareBox() {

}


void LeesEdwardsSquareBox::get_settings(input_file &inp) {
	number dt, shear_rate;
	getInputNumber(&inp, "dt", &dt, 1);
	getInputNumber(&inp, "lees_edwards_shear_rate", &shear_rate, 1);
	factor = dt * shear_rate;
}


void LeesEdwardsSquareBox::init(number Lx, number Ly) {
	SquareBox::init(Lx, Ly);

	factor *= Ly;
}


int LeesEdwardsSquareBox::getElementX(int site, int distX){
	int x = (site%this->side)+distX;
	while(x<0)x+=this->side;
	while(x>=this->side)x-=this->side;
	return x;
}


int LeesEdwardsSquareBox::getElementY(int site, int distY){
	int y = (site/this->side)+distY;
        while(y<0)y+=this->side;
        while(y>=this->side)y-=this->side;
	return y;
}


std::vector<number> LeesEdwardsSquareBox::min_image(const std::vector<number> &v1, const std::vector<number> &v2) const {
	static llint last_step = -1;
	static number delta_x = 0.;

	llint curr_step = CONFIG_INFO->curr_step;
	if(curr_step != last_step) {
		delta_x = factor * curr_step;
		last_step = curr_step;
	}

	number ny = v2[1] - v1[1];
	number cy = rint(ny / this->side);
	number nx = v2[0] - v1[0] - cy * delta_x;

	return std::vector<number> (nx - rint(nx / this->side) * this->side, ny - cy * this->side);
}


number LeesEdwardsSquareBox::sqr_min_image_distance(const std::vector<number> &v1, const std::vector<number> &v2) const {
	std::vector<number> v = min_image(v1, v2);
	return v[0] * v[0] + v[1] * v[1]; 
}

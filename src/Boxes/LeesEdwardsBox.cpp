#include "LeesEdwardsBox.h"
#include "../Fields/BaseField.h"

LeesEdwardsBox::LeesEdwardsBox() {
lees_edwards=true;
}


LeesEdwardsBox::~LeesEdwardsBox() {
lees_edwards=true;
}


void LeesEdwardsBox::get_settings(input_file &inp) {
	number dt, shear_rate;
	getInputNumber(&inp, "dt", &dt, 1);
	getInputNumber(&inp, "lees_edwards_shear_rate", &shear_rate, 1);
	factor = dt * shear_rate;
}


void LeesEdwardsBox::init(int Lx, int Ly) {
	if(Lx != Ly) throw RCexception("The box in the configuration file is not a square");

	side = Lx;
	sides[0] = sides[1]= Lx;

	neighbors.resize(Lx*Ly*9);
	BaseBox::setNeighborsPeriodic(Lx, Ly);

	factor *= Ly;
}


int LeesEdwardsBox::getElementX(int site, int distX){
	int x = (site%this->side)+distX;
	while(x<0)x+=this->side;
	while(x>=this->side)x-=this->side;
	return x;
}


int LeesEdwardsBox::getElementY(int site, int distY){
	int y = (site/this->side)+distY;
        while(y<0)y+=this->side;
        while(y>=this->side)y-=this->side;
	return y;
}

int LeesEdwardsBox::getXsize(){return side;}
int LeesEdwardsBox::getYsize(){return side;}
int LeesEdwardsBox::getBoxsize(){return side*side;}

std::vector<number> LeesEdwardsBox::get_abs_pos(BaseField *p) {
	return std::vector<number> {p->CoM[0] + p->pos_shift[0], p->CoM[1] + p->pos_shift[1]};
}

void LeesEdwardsBox::shift_particle(BaseField * p, std::vector<number> &amount) {
        p->pos_shift[0] += amount[0];
        p->pos_shift[1] += amount[1];
        p->CoM[0] += p->pos_shift[0];
        p->CoM[1] += p->pos_shift[1];
}


std::vector<number> LeesEdwardsBox::normalised_in_box(const std::vector<number> &v) {
        return std::vector<number> {
                (v[0] / side - floor(v[0] / side)),
                (v[1] / side - floor(v[1] / side)) 
        };
}


std::vector<number> LeesEdwardsBox::min_image(const std::vector<number> &v1, const std::vector<number> &v2) const {
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


number LeesEdwardsBox::sqr_min_image_distance(const std::vector<number> &v1, const std::vector<number> &v2) const {
	std::vector<number> v = min_image(v1, v2);
	return v[0] * v[0] + v[1] * v[1]; 
}

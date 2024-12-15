#include "OrthogonalBox.h"

#include "../Fields/BaseField.h"

using namespace std;

OrthogonalBox::OrthogonalBox() {
	sidex = -1.0;
	sidey = -1.0;
	sides.resize(2);
}

OrthogonalBox::~OrthogonalBox() {

}

void OrthogonalBox::get_settings(input_file &inp) {

}


void OrthogonalBox::init(int Lx, int Ly) {
	if(Lx == Ly) throw RCexception("Using squared orthogonal box...just use square box...");

	sidex = Lx;
	sidey = Ly;
	sides[0] = Lx;
	sides[1] = Ly;

	neighbors.resize(Lx*Ly*9);
	BaseBox::setNeighborsPeriodic(Lx, Ly);
}

int OrthogonalBox::getElementX(int site, int distX){
	int x = (site-(int(site/sidex)*sidex))+distX;
	while(x<0)x+=sidex;
	while(x>=sidex)x-=sidex;
	return x;
//(((site/side)+distX)/side - floor( ((site/side)+distX) / side) ) * side;
}

int OrthogonalBox::getElementY(int site, int distY){
	int y = (site/sidex)+distY;
        while(y<0)y+=sidey;
        while(y>=sidey)y-=sidey;
	return y;
//(((site%side)+distY)/side - floor( ((site%side)+distY) / side) ) * side;
}

int OrthogonalBox::getXsize(){return sidex;}
int OrthogonalBox::getYsize(){return sidey;}
int OrthogonalBox::getBoxsize(){return sidex*sidey;}

std::vector<number> OrthogonalBox::get_abs_pos(BaseField *p) {
	return vector<number> {p->CoM[0] + p->pos_shift[0], p->CoM[1] + p->pos_shift[1]};
}

void OrthogonalBox::shift_particle(BaseField * p, std::vector<number> &amount) {
        p->pos_shift[0] += amount[0];
        p->pos_shift[1] += amount[1];
        p->CoM[0] += p->pos_shift[0];
        p->CoM[1] += p->pos_shift[1];
}


std::vector<number> OrthogonalBox::normalised_in_box(const std::vector<number> &v) {
        return std::vector<number> {
                (v[0] / sidex - floor(v[0] / sidex)) * sidex,
                (v[1] / sidey - floor(v[1] / sidey)) * sidey
        };
}


inline std::vector<number> OrthogonalBox::min_image(const std::vector<number> &v1, const std::vector<number> &v2) const {
        return std::vector<number> {
                v2[0] - v1[0] - rint((v2[0] - v1[0]) / sidex) * sidex,
                v2[1] - v1[1] - rint((v2[1] - v1[1]) / sidey) * sidey
        };
}

number OrthogonalBox::sqr_min_image_distance(const std::vector<number> &v1, const std::vector<number> &v2) const {
	std::vector<number> v = min_image(v1, v2);
        return v[0]*v[0]+v[1]*v[1];
}



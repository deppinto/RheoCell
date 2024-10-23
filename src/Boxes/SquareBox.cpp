#include "SquareBox.h"

#include "../Fields/BaseField.h"

using namespace std;

SquareBox::SquareBox() {
	side = -1.0;
	sides.resize(2);
}

SquareBox::~SquareBox() {

}

void SquareBox::get_settings(input_file &inp) {

}


void SquareBox::init(int Lx, int Ly) {
	if(Lx != Ly) throw RCexception("The box in the configuration file is not a square");

	side = Lx;
	sides[0] = sides[1]= Lx;

	neighbors.resize(Lx*Ly*9);
	BaseBox::setNeighborsPeriodic(Lx, Ly);
}

int SquareBox::getElementX(int site, int distX){
	int x = (site-(int(site/side)*side))+distX;
	while(x<0)x+=side;
	while(x>=side)x-=side;
	return x;
//(((site/side)+distX)/side - floor( ((site/side)+distX) / side) ) * side;
}

int SquareBox::getElementY(int site, int distY){
	int y = (site/side)+distY;
        while(y<0)y+=side;
        while(y>=side)y-=side;
	return y;
//(((site%side)+distY)/side - floor( ((site%side)+distY) / side) ) * side;
}

int SquareBox::getXsize(){return side;}
int SquareBox::getYsize(){return side;}
int SquareBox::getBoxsize(){return side*side;}

std::vector<number> SquareBox::get_abs_pos(BaseField *p) {
	return vector<number> {p->CoM[0] + p->pos_shift[0], p->CoM[1] + p->pos_shift[1]};
}

void SquareBox::shift_particle(BaseField * p, std::vector<number> &amount) {
        p->pos_shift[0] += amount[0];
        p->pos_shift[1] += amount[1];
        p->CoM[0] += p->pos_shift[0];
        p->CoM[1] += p->pos_shift[1];
}


std::vector<number> SquareBox::normalised_in_box(const std::vector<number> &v) {
        return std::vector<number> {
                (v[0] / side - floor(v[0] / side)),
                (v[1] / side - floor(v[1] / side)) 
        };
}


inline std::vector<number> SquareBox::min_image(const std::vector<number> &v1, const std::vector<number> &v2) const {
        return std::vector<number> {
                v2[0] - v1[0] - rint((v2[0] - v1[0]) / side) * side,
                v2[1] - v1[1] - rint((v2[1] - v1[1]) / side) * side
        };
}

number SquareBox::sqr_min_image_distance(const std::vector<number> &v1, const std::vector<number> &v2) const {
	std::vector<number> v = min_image(v1, v2);
        return v[0]*v[0]+v[1]*v[1];
}



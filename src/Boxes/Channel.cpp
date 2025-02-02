#include "Channel.h"

#include "../Fields/BaseField.h"

using namespace std;

Channel::Channel() {
	sidex = -1.0;
	sidey = -1.0;
	sides.resize(2);
}

Channel::~Channel() {

}

void Channel::get_settings(input_file &inp) {

	getInputNumber(&inp, "lambda_wall", &lambda_wall, 0);
}


void Channel::init(int Lx, int Ly) {
	if(Lx == Ly) throw RCexception("Using squared orthogonal box...just use square box...");

	sidex = Lx;
	sidey = Ly;
	sides[0] = Lx;
	sides[1] = Ly;

	walls.resize(Lx*Ly);
	for(int y=0; y<Ly; y++){
		for(int x=0; x<Lx; x++){
			int k=x+y*Lx;
			walls[k] = exp(-double(y)/lambda_wall) + exp(-double(Ly-y-1)/lambda_wall);  
		}
	}

	neighbors.resize(Lx*Ly*9);
	BaseBox::setNeighborsPeriodic(Lx, Ly);

	laplacian_walls.resize(Lx*Ly);
	for(int y=0; y<Ly; y++){
		for(int x=0; x<Lx; x++){
			int k=x+y*Lx;
			laplacian_walls[k] = walls[neighbors[5+k*9]] + walls[neighbors[7+k*9]] + walls[neighbors[3+k*9]] + walls[neighbors[1+k*9]] - 4.*walls[k];
		}
	}
}


number Channel::getWalls(int k) {return walls[k];}
number Channel::getLaplacianWalls(int k) {return laplacian_walls[k];}


int Channel::getElementX(int site, int distX){
	int x = (site-(int(site/sidex)*sidex))+distX;
	while(x<0)x+=sidex;
	while(x>=sidex)x-=sidex;
	return x;
//(((site/side)+distX)/side - floor( ((site/side)+distX) / side) ) * side;
}

int Channel::getElementY(int site, int distY){
	int y = (site/sidex)+distY;
        while(y<0)y+=sidey;
        while(y>=sidey)y-=sidey;
	return y;
//(((site%side)+distY)/side - floor( ((site%side)+distY) / side) ) * side;
}

int Channel::getXsize(){return sidex;}
int Channel::getYsize(){return sidey;}
int Channel::getBoxsize(){return sidex*sidey;}

std::vector<number> Channel::get_abs_pos(BaseField *p) {
	return vector<number> {p->CoM[0] + p->pos_shift[0], p->CoM[1] + p->pos_shift[1]};
}

void Channel::shift_particle(BaseField * p, std::vector<number> &amount) {
        p->pos_shift[0] += amount[0];
        p->pos_shift[1] += amount[1];
        p->CoM[0] += p->pos_shift[0];
        p->CoM[1] += p->pos_shift[1];
}


std::vector<number> Channel::normalised_in_box(const std::vector<number> &v) {
        return std::vector<number> {
                (v[0] / sidex - floor(v[0] / sidex))*sidex,
                (v[1] / sidey - floor(v[1] / sidey))*sidey 
        };
}


inline std::vector<number> Channel::min_image(const std::vector<number> &v1, const std::vector<number> &v2) const {
        return std::vector<number> {
                v2[0] - v1[0] - rint((v2[0] - v1[0]) / sidex) * sidex,
                v2[1] - v1[1] - rint((v2[1] - v1[1]) / sidey) * sidey
        };
}

number Channel::sqr_min_image_distance(const std::vector<number> &v1, const std::vector<number> &v2) const {
	std::vector<number> v = min_image(v1, v2);
        return v[0]*v[0]+v[1]*v[1];
}



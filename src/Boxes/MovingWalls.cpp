#include "MovingWalls.h"

#include "../Fields/BaseField.h"

using namespace std;

MovingWalls::MovingWalls() {
	sidex = -1.0;
	sidey = -1.0;
	sides.resize(2);
}

MovingWalls::~MovingWalls() {

}

void MovingWalls::get_settings(input_file &inp) {
	number dt, shear_rate;
	getInputNumber(&inp, "dt", &dt, 1);
	getInputNumber(&inp, "lees_edwards_shear_rate", &shear_rate, 1);
	factor = dt * shear_rate;

	getInputNumber(&inp, "lambda_wall", &lambda_wall, 0);
}


void MovingWalls::init(int Lx, int Ly) {
	if(Lx == Ly) throw RCexception("Using squared orthogonal box...just use square box...");

	sidex = Lx;
	sidey = Ly;
	sides[0] = Lx;
	sides[1] = Ly;
	factor *= Ly;

	walls.resize(Lx*Ly);
	for(int y=0; y<Ly; y++){
		for(int x=0; x<Lx; x++){
			int k=x+y*Lx;
			walls[k] = exp(-double(x)/lambda_wall) + exp(-double(Lx-x-1)/lambda_wall);  
			//walls[k] += exp(-double(y)/16) + exp(-double(Ly-y-1)/16);  
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


number MovingWalls::getWalls(int k) {return walls[k];}
number MovingWalls::getLaplacianWalls(int k) {return laplacian_walls[k];}


int MovingWalls::getElementX(int site, int distX){
	int x = (site-(int(site/sidex)*sidex))+distX;
	while(x<0)x+=sidex;
	while(x>=sidex)x-=sidex;
	return x;
//(((site/side)+distX)/side - floor( ((site/side)+distX) / side) ) * side;
}

int MovingWalls::getElementY(int site, int distY){
	int y = (site/sidex)+distY;
        while(y<0)y+=sidey;
        while(y>=sidey)y-=sidey;
	return y;
//(((site%side)+distY)/side - floor( ((site%side)+distY) / side) ) * side;
}

int MovingWalls::getXsize(){return sidex;}
int MovingWalls::getYsize(){return sidey;}
int MovingWalls::getBoxsize(){return sidex*sidey;}

std::vector<number> MovingWalls::get_abs_pos(BaseField *p) {
	return vector<number> {p->CoM[0] + p->pos_shift[0], p->CoM[1] + p->pos_shift[1]};
}

void MovingWalls::shift_particle(BaseField * p, std::vector<number> &amount) {
        p->pos_shift[0] += amount[0];
        p->pos_shift[1] += amount[1];
        p->CoM[0] += p->pos_shift[0];
        p->CoM[1] += p->pos_shift[1];
}


std::vector<number> MovingWalls::normalised_in_box(const std::vector<number> &v) {
        return std::vector<number> {
                (v[0] / sidex - floor(v[0] / sidex))*sidex,
                (v[1] / sidey - floor(v[1] / sidey))*sidey 
        };
}


inline std::vector<number> MovingWalls::min_image(const std::vector<number> &v1, const std::vector<number> &v2) const {
        return std::vector<number> {
                v2[0] - v1[0] - rint((v2[0] - v1[0]) / sidex) * sidex,
                v2[1] - v1[1] - rint((v2[1] - v1[1]) / sidey) * sidey
        };
}

number MovingWalls::sqr_min_image_distance(const std::vector<number> &v1, const std::vector<number> &v2) const {
	std::vector<number> v = min_image(v1, v2);
        return v[0]*v[0]+v[1]*v[1];
}


void MovingWalls::UpdateWalls(bool wallflag) {

	llint curr_step = CONFIG_INFO->curr_step;
	if(curr_step != last_step) {
		delta_x = factor * curr_step;
		last_step = curr_step;
	}
	else return;

	number updated_lambda = lambda_wall - delta_x;
	if(updated_lambda<0)updated_lambda=0;

	int Lx = sides[0];
	int Ly = sides[1];

	for(int y=0; y<Ly; y++){
		for(int x=0; x<Lx; x++){
			int k=x+y*Lx;
			walls[k] = exp(-double(x)/updated_lambda) + exp(-double(Lx-x-1)/updated_lambda);  
		}
	}

	for(int y=0; y<Ly; y++){
		for(int x=0; x<Lx; x++){
			int k=x+y*Lx;
			laplacian_walls[k] = walls[neighbors[5+k*9]] + walls[neighbors[7+k*9]] + walls[neighbors[3+k*9]] + walls[neighbors[1+k*9]] - 4.*walls[k];
		}
	}


}

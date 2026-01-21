#include "SquareWalls.h"

#include "../Fields/BaseField.h"

using namespace std;

SquareWalls::SquareWalls() {
	sidex = -1.0;
	sidey = -1.0;
	sides.resize(2);
}

SquareWalls::~SquareWalls() {

}

void SquareWalls::get_settings(input_file &inp) {

	getInputNumber(&inp, "lambda_wall", &lambda_wall, 0);
}


void SquareWalls::init(int Lx, int Ly) {
	if(Lx != Ly) throw RCexception("Using orthogonal sides for square box!");

	sidex = Lx;
	sidey = Ly;
	sides[0] = Lx;
	sides[1] = Ly;

	walls.resize(Lx*Ly);
	int empty_area=Lx*Ly;
	for(int y=0; y<Ly; y++){
		for(int x=0; x<Lx; x++){
			int k=x+y*Lx;
			/*double wallsx=exp(-double(x)/lambda_wall) + exp(-double(Lx-x-1)/lambda_wall);
			double wallsy=exp(-double(y)/lambda_wall) + exp(-double(Ly-y-1)/lambda_wall);
			if(wallsx>wallsy)walls[k]=wallsx;
			else walls[k]=wallsy;*/
			if(x<lambda_wall-1 || x>Lx-lambda_wall){
				if(walls[k]==0)empty_area-=1;
				walls[k] = 1.; 
			}
			if(y<lambda_wall || y>Ly-lambda_wall){
				if(walls[k]==0)empty_area-=1;
				walls[k] = 1.; 
			}
		}
	}

	std::cout<<"TESTING: Walls empty area: "<<empty_area<<std::endl;
	neighbors.resize(Lx*Ly*9);
	BaseBox::setNeighborsPeriodic(Lx, Ly);
}


number SquareWalls::getWalls(int k) {return walls[k];}


int SquareWalls::getElementX(int site, int distX){
	int x = (site-(int(site/sidex)*sidex))+distX;
	while(x<0)x+=sidex;
	while(x>=sidex)x-=sidex;
	return x;
//(((site/side)+distX)/side - floor( ((site/side)+distX) / side) ) * side;
}

int SquareWalls::getElementY(int site, int distY){
	int y = (site/sidex)+distY;
        while(y<0)y+=sidey;
        while(y>=sidey)y-=sidey;
	return y;
//(((site%side)+distY)/side - floor( ((site%side)+distY) / side) ) * side;
}

int SquareWalls::getXsize(){return sidex;}
int SquareWalls::getYsize(){return sidey;}
int SquareWalls::getBoxsize(){return sidex*sidey;}

std::vector<number> SquareWalls::get_abs_pos(BaseField *p) {
	return vector<number> {p->CoM[0] + p->pos_shift[0], p->CoM[1] + p->pos_shift[1]};
}

void SquareWalls::shift_particle(BaseField * p, std::vector<number> &amount) {
        p->pos_shift[0] += amount[0];
        p->pos_shift[1] += amount[1];
        p->CoM[0] += p->pos_shift[0];
        p->CoM[1] += p->pos_shift[1];
}


std::vector<number> SquareWalls::normalised_in_box(const std::vector<number> &v) {
        return std::vector<number> {
                (v[0] / sidex - floor(v[0] / sidex))*sidex,
                (v[1] / sidey - floor(v[1] / sidey))*sidey 
        };
}


inline std::vector<number> SquareWalls::min_image(const std::vector<number> &v1, const std::vector<number> &v2) const {
        return std::vector<number> {
                v2[0] - v1[0] - rint((v2[0] - v1[0]) / sidex) * sidex,
                v2[1] - v1[1] - rint((v2[1] - v1[1]) / sidey) * sidey
        };
}

number SquareWalls::sqr_min_image_distance(const std::vector<number> &v1, const std::vector<number> &v2) const {
	std::vector<number> v = min_image(v1, v2);
        return v[0]*v[0]+v[1]*v[1];
}



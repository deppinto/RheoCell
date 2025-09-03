#include "TriangleWalls.h"

#include "../Fields/BaseField.h"

using namespace std;

/*
 * Sides for Vania experiment: Circle (8), Square (14), Triangle (15)
 */


TriangleWalls::TriangleWalls() {
	sidex = -1.0;
	sidey = -1.0;
	sides.resize(2);
}

TriangleWalls::~TriangleWalls() {

}

void TriangleWalls::get_settings(input_file &inp) {

	getInputNumber(&inp, "lambda_wall", &lambda_wall, 0);
}


void TriangleWalls::init(int Lx, int Ly) {
	if(Lx != Ly) throw RCexception("Using orthogonal sides for square box!");

	sidex = Lx;
	sidey = Ly;
	sides[0] = Lx;
	sides[1] = Ly;

	int side = Lx - lambda_wall * 2;
    	int height = int(round(side * sqrt(3.0) / 2.0));
	int center = Lx / 2;
	int top = Lx / 2 - height / 2;

	walls.resize(Lx*Ly);
	for(int y=0; y<Ly; y++){
		for(int x=0; x<Lx; x++){
			walls[int(x + y * Lx)] = 1.;
		}
	}

	int empty_area=0;
	for (int row = 0; row < height; ++row) {
		int halfWidth = int((side / 2.0) * (1.0 - (double)row / height));
		int left = center - halfWidth;
		int right = center + halfWidth;
		
		int y = top + row;
		if (y >= 0 && y < Ly) {
			for (int x = left; x <= right && x < Lx; ++x) {
				if (x >= 0) {
					if(walls[int(x + y * Lx)]==1)empty_area+=1;
					walls[int(x + y * Lx)] = 0.;
				}
			}
		}
	}


	std::cout<<"TESTING: Walls empty area: "<<empty_area<<std::endl;

	/*
	walls.resize(Lx*Ly);
	for(int y=0; y<Ly; y++){
		for(int x=0; x<Lx; x++){
			int k=x+y*Lx;
			double wallsx=exp(-double(x)/lambda_wall) + exp(-double(Lx-x-1)/lambda_wall);
			double wallsy=exp(-double(y)/lambda_wall) + exp(-double(Ly-y-1)/lambda_wall);
			if(wallsx>wallsy)walls[k]=wallsx;
			else walls[k]=wallsy;
		}
	}*/

	neighbors.resize(Lx*Ly*9);
	BaseBox::setNeighborsPeriodic(Lx, Ly);
}


number TriangleWalls::getWalls(int k) {return walls[k];}


int TriangleWalls::getElementX(int site, int distX){
	int x = (site-(int(site/sidex)*sidex))+distX;
	while(x<0)x+=sidex;
	while(x>=sidex)x-=sidex;
	return x;
//(((site/side)+distX)/side - floor( ((site/side)+distX) / side) ) * side;
}

int TriangleWalls::getElementY(int site, int distY){
	int y = (site/sidex)+distY;
        while(y<0)y+=sidey;
        while(y>=sidey)y-=sidey;
	return y;
//(((site%side)+distY)/side - floor( ((site%side)+distY) / side) ) * side;
}

int TriangleWalls::getXsize(){return sidex;}
int TriangleWalls::getYsize(){return sidey;}
int TriangleWalls::getBoxsize(){return sidex*sidey;}

std::vector<number> TriangleWalls::get_abs_pos(BaseField *p) {
	return vector<number> {p->CoM[0] + p->pos_shift[0], p->CoM[1] + p->pos_shift[1]};
}

void TriangleWalls::shift_particle(BaseField * p, std::vector<number> &amount) {
        p->pos_shift[0] += amount[0];
        p->pos_shift[1] += amount[1];
        p->CoM[0] += p->pos_shift[0];
        p->CoM[1] += p->pos_shift[1];
}


std::vector<number> TriangleWalls::normalised_in_box(const std::vector<number> &v) {
        return std::vector<number> {
                (v[0] / sidex - floor(v[0] / sidex))*sidex,
                (v[1] / sidey - floor(v[1] / sidey))*sidey 
        };
}


inline std::vector<number> TriangleWalls::min_image(const std::vector<number> &v1, const std::vector<number> &v2) const {
        return std::vector<number> {
                v2[0] - v1[0] - rint((v2[0] - v1[0]) / sidex) * sidex,
                v2[1] - v1[1] - rint((v2[1] - v1[1]) / sidey) * sidey
        };
}

number TriangleWalls::sqr_min_image_distance(const std::vector<number> &v1, const std::vector<number> &v2) const {
	std::vector<number> v = min_image(v1, v2);
        return v[0]*v[0]+v[1]*v[1];
}



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

	sidex = Lx;
	sidey = Ly;
	sides[0] = Lx;
	sides[1] = Ly;
	factor *= Ly;

	neighbors.resize(Lx*Ly*9);
	neighbors_next.resize(Lx*Ly*9);
	weight_site.resize(Lx*Ly*9);
	weight_site_next.resize(Lx*Ly*9);
	setNeighborsPeriodic(Lx, Ly);
}

int LeesEdwardsBox::getElementX(int site, int distX){
	int x = (site-(int(site/sidex)*sidex))+distX;
	while(x<0)x+=sidex;
	while(x>=sidex)x-=sidex;
	return x;
}

int LeesEdwardsBox::getElementY(int site, int distY){
	int y = (site/sidex)+distY;
        while(y<0)y+=sidey;
        while(y>=sidey)y-=sidey;
	return y;
}

int LeesEdwardsBox::getElement(int site, int distX, int distY){

	int cross=0;

	int x = (site-(int(site/sidex)*sidex))+distX;
	while(x<0)x+=sidex;
	while(x>=sidex)x-=sidex;

	int y = (site/sidex)+distY;
        while(y<0){y+=sidey;cross=1;}
        while(y>=sidey){y-=sidey;cross=2;}

	if(cross==1){
		x = int(x + delta_x);
		while(x>=sidex)x-=sidex;
	}
	else if(cross==2){
		x = int(x - delta_x);
		while(x<0)x+=sidex;
	}

	return x + y * sidex;
}

int LeesEdwardsBox::getXsize(){return sidex;}
int LeesEdwardsBox::getYsize(){return sidey;}
int LeesEdwardsBox::getBoxsize(){return sidex*sidey;}

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
                (v[0] / sidex - floor(v[0] / sidex)),
                (v[1] / sidey - floor(v[1] / sidey)) 
        };
}

std::vector<number> LeesEdwardsBox::min_image(const std::vector<number> &v1, const std::vector<number> &v2) const {

	number ny = v2[1] - v1[1];
	number cy = rint(ny / sidey);
	number nx = v2[0] - v1[0] - cy * delta_x;

	return std::vector<number> (nx - rint(nx / sidex) * sidex, ny - cy * sidey);
}

number LeesEdwardsBox::sqr_min_image_distance(const std::vector<number> &v1, const std::vector<number> &v2) const {
	std::vector<number> v = min_image(v1, v2);
	return v[0] * v[0] + v[1] * v[1]; 
}

void LeesEdwardsBox::setNeighborsPeriodic(int Lx, int Ly){

	llint curr_step = CONFIG_INFO->curr_step;
	if(curr_step != last_step) {
		delta_x = factor * curr_step;
		last_step = curr_step;
	}

	//initialize neighbors of lattice sites
	int x,y,xx,yy,ss;
	int x_next;
	number weight;
	number x_shifted;
	for(int i =0; i<Lx*Ly; i++) {
		y = i/Lx;
		x = i-(int(i/Lx)*Lx);
		ss=0;
		for(int j=-1; j<=1; j++){
			for(int k=-1; k<=1; k++){
				xx=x+k;
				yy=y+j;

				x_next=0;
				weight=0.;

				if(x==0 && k==-1){
					xx=Lx-1;
				}
				else if(x==Lx-1 && k==1){
					xx=0;
				}

				if(y==0 && j==-1){
					x_shifted = (double)xx + delta_x;
					while(x_shifted>=sidex)x_shifted-=sidex;
					xx = (int)x_shifted;
					x_next = xx + 1;
					if(x_next>=sidex)x_next-=sidex;
					weight = x_shifted - (double)xx;
					yy = Ly-1;
				}
				else if(y==Ly-1 && j==1){
					x_shifted = (double)xx - delta_x;
					while(x_shifted<0)x_shifted+=sidex;
					xx = (int)x_shifted;
					x_next = xx + 1;
					if(x_next>=sidex)x_next-=sidex;
					weight = x_shifted - (double)xx;
					yy = 0;
				}

				neighbors[ss+i*9]=xx+yy*Lx;
				neighbors_next[ss+i*9]=x_next+yy*Lx;
				weight_site[ss+i*9]=1-weight;
				weight_site_next[ss+i*9]=weight;

				ss++;
			}
		}
	}
}

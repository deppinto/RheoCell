#include "BaseBox.h"

#include "../Fields/BaseField.h"

BaseBox::BaseBox() {

}

BaseBox::~BaseBox() {

}

std::vector<number> BaseBox::min_image(const BaseField &p, const BaseField &q) {
        return min_image(p.CoM, q.CoM);
}

std::vector<number> BaseBox::min_image(const BaseField *p, const BaseField *q) {
        return min_image(p->CoM, q->CoM);
}

number BaseBox::sqr_min_image_distance(const BaseField &p, const BaseField &q) {
        return sqr_min_image_distance(p.CoM, q.CoM);
}

number BaseBox::sqr_min_image_distance(const BaseField *p, const BaseField *q) {
        return sqr_min_image_distance(p->CoM, q->CoM);
}

void BaseBox::setNeighborsPeriodic(int Lx, int Ly){
	//initialize neighbors of lattice sites
	int x,y,xx,yy,site,ss;
	for(int i =0; i<Lx*Ly; i++) {
		y = i/Lx;
		x = i-(int(i/Lx)*Lx);
		ss=0;
		for(int j=-1; j<=1; j++){
			for(int k=-1; k<=1; k++){
				xx=x+k;
				yy=y+j;
				if(x==0 && k==-1)site=Lx-1+yy*Lx;
				else if(x==Lx-1 && k==1)site=yy*Lx;
				else if(y==0 && j==-1)site=xx+(Ly-1)*Lx;
				else if(y==Ly-1 && j==1)site=xx;
				else site=xx+yy*Lx;
				neighbors[ss+i*9]=site;
				ss++;
			}
		}
	}
}

#include "MultiPhaseField.h"


MultiPhaseField::MultiPhaseField() :  BaseField() {

}


MultiPhaseField::MultiPhaseField(int Lx, int Ly) : BaseField() {

	LsubX=Lx;
	LsubY=Ly;
	subSize=Lx*Ly;

}

MultiPhaseField::~MultiPhaseField() {

}

void MultiPhaseField::get_interaction_values(int R) {
        BaseField::get_interaction_values(R);
	init_radius=R/2;
	init_radius2=init_radius*init_radius;
}


void MultiPhaseField::init(int Lx, int Ly) {
        LsubX = Lx;
	LsubY = Ly;
	subSize = Lx*Ly;
	fieldScalar.resize(subSize);
	fieldDX.resize(subSize);
	fieldDY.resize(subSize);	
	freeEnergy.resize(subSize);
	fieldScalar_old.resize(subSize);
	dfield_old.resize(subSize);
	neighbors_sub.resize(subSize*9);
	velocityX=0.;
        velocityY=0.;
	S00=0.;
        S01=0.;
	Fpressure= std::vector<number> {0.,0.};
	Fshape= std::vector<number> {0.,0.};
	area=0;
	sumF=0;
	offset.resize(2);
	offset[0]=0; offset[1]=0;
        //int center=LsubX/2+(LsubY/2)*LsubX;
        /*int x,y;
        for(int i=0; i<subSize; i++){
		freeEnergy[i]=0;
                x=LsubX/2-i%LsubX;
                y=LsubY/2-i/LsubY;
                if(x*x+y*y<init_radius2){fieldScalar[i]=1.;area+=1.;sumF+=1;}
                else fieldScalar[i]=0.;
        }*/
	setNeighborsSub();
}


void MultiPhaseField::init() {
        subSize = LsubX*LsubY;
        fieldScalar.resize(subSize);
	fieldDX.resize(subSize);
        fieldDY.resize(subSize);
	freeEnergy.resize(subSize);
	fieldScalar_old.resize(subSize);
        dfield_old.resize(subSize);
	neighbors_sub.resize(subSize*9);
	velocityX=0.;
        velocityY=0.;
	S00=0.;
        S01=0.;
	Fpressure= std::vector<number> {0.,0.};
        Fshape= std::vector<number> {0.,0.};
	area=0;
	offset.resize(2);
	offset[0]=0; offset[1]=0;
	//int center=LsubX/2+(LsubY/2)*LsubX;
	int x,y;
	for(int i=0; i<subSize; i++){
		x=LsubX/2-i%LsubX;
		y=LsubY/2-i/LsubY;
		if(x*x+y*y<init_radius2+0.5){fieldScalar[i]=1.;area+=1;sumF+=1;}
		else fieldScalar[i]=0.;
	}
	setNeighborsSub();
}

void MultiPhaseField::setNeighborsSub() {
	//setNeighborsSubSquareDirichlet();
	setNeighborsSubSquarePeriodic();
}


void MultiPhaseField::setNeighborsSubSquareDirichlet() {
	int x,y,xx,yy,site,ss;
	for(int i =0; i<subSize; i++) {
		x=i%LsubX;
		y=i/LsubY;
		ss=0;
		for(int j=-1; j<=1; j++){
			for(int k=-1; k<=1; k++){
				xx=x+k;
				yy=y+j;
				if(x==0 && k==-1)site=-1;
				else if(x==LsubX-1 && k==1)site=-1;
				else if(y==0 && j==-1)site=-1;
				else if(y==LsubY-1 && j==1)site=-1;
				else site=xx+yy*LsubX;
				neighbors_sub[ss+i*9]=site;
				ss++;
			}
		}
	}
}

void MultiPhaseField::setNeighborsSubSquarePeriodic() {
	int x,y,xx,yy,site,ss;
	for(int i =0; i<subSize; i++) {
		x=i%LsubX;
		y=i/LsubY;
		ss=0;
		for(int j=-1; j<=1; j++){
			for(int k=-1; k<=1; k++){
				xx=x+k;
				yy=y+j;
				if(x==0 && k==-1)site=LsubX-1+yy*LsubX;
				else if(x==LsubX-1 && k==1)site=yy*LsubX;
				else if(y==0 && j==-1)site=xx+(LsubY-1)*LsubX;
				else if(y==LsubY-1 && j==1)site=xx;
				else site=xx+yy*LsubX;
				neighbors_sub[ss+i*9]=site;
				ss++;
			}
		}
	}
}


void MultiPhaseField::set_positions_initial(BaseBox *box) {
	int x=CoM[0];
	int y=CoM[1];
	int site=x+y*box->getXsize();
	sub_corner_bottom_left=box->getElementX(site, (int)-LsubX/2) + box->getElementY(site, (int)-LsubY/2) * box->getXsize();
}

void MultiPhaseField::set_positions(BaseBox *box) {
	int x=CoM[0];
	int y=CoM[1];
	int site=x+y*box->getXsize();
	std::vector<number> corner_old = std::vector<number> { (number)box->getElementX(sub_corner_bottom_left, 0) , (number)box->getElementY(sub_corner_bottom_left, 0) };
	sub_corner_bottom_left=box->getElementX(site, (int)-LsubX/2) + box->getElementY(site, (int)-LsubY/2) * box->getXsize();
	
	std::vector<number> corner_new = std::vector<number> { (number)box->getElementX(sub_corner_bottom_left, 0) , (number)box->getElementY(sub_corner_bottom_left, 0) };
	std::vector<number> displacement = box->min_image(corner_old, corner_new);
	offset[0] =  int(offset[0] + LsubX - displacement[0])%LsubX; offset[1] =  int(offset[1] + LsubY - displacement[1])%LsubY;
}

void MultiPhaseField::set_positions(int offsetx, int offsety, int corner) {
	sub_corner_bottom_left=corner;
	offset[0]=offsetx;offset[1]=offsety;
}

/*transform subdomain sites into grid sites*/

int MultiPhaseField::GetSubIndex(int site, BaseBox *box){
	return box->getElementX(sub_corner_bottom_left, int((site%LsubX)+offset[0])%LsubX) + box->getElementY(sub_corner_bottom_left, int((site/LsubY)+offset[1])%LsubY) * box->getXsize();
}

int MultiPhaseField::GetSubXIndex(int site, BaseBox *box){
	return box->getElementX(sub_corner_bottom_left, int((site%LsubX)+offset[0])%LsubX);
}

int MultiPhaseField::GetSubYIndex(int site, BaseBox *box){
	return box->getElementY(sub_corner_bottom_left, int((site/LsubY)+offset[1])%LsubY);
}

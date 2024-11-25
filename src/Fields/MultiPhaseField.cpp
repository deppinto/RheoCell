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

void MultiPhaseField::resizing() {
	fieldScalar.resize(subSize);
	fieldDX.resize(subSize);
	fieldDY.resize(subSize);	
	freeEnergy.resize(subSize);
	fieldScalar_old.resize(subSize);
	dfield_old.resize(subSize);
	neighbors_sub.resize(subSize*9);
	velocityX.resize(subSize);
        velocityY.resize(subSize);
	map_sub_to_box.resize(subSize);
	map_sub_to_box_x.resize(subSize);
	map_sub_to_box_y.resize(subSize);
}

void MultiPhaseField::init(int Lx, int Ly) {
        LsubX = Lx;
	LsubY = Ly;
	subSize = Lx*Ly;
	x_sub_left=LsubX;
	y_sub_bottom=LsubY;

	resizing();
	S00=0.;
        S01=0.;

	thetaQ = thetaQ_old = PI * (1-2*drand48());
	Q00=cos(2*thetaQ);
	Q01=sin(2*thetaQ);
	nemQ.resize(2);
	nemQ_old.resize(2);
	nemQ = {0., 0.};
	nemQ_old = {0., 0.};
	Q00 = 0.5 * (nemQ[0] * nemQ[0] - nemQ[1] * nemQ[1]);
	Q01 = nemQ[0] * nemQ[1];

	Fpassive= std::vector<number> {0.,0.};
	Factive= std::vector<number> {0.,0.};
	area=0;
	sumF=0;
	offset.resize(2);
	offset[0]=0; offset[1]=0;
	setNeighborsSub();
}


void MultiPhaseField::init() {
        subSize = LsubX*LsubY;
	x_sub_left=LsubX;
	y_sub_bottom=LsubY;

	//standard field properties
	resizing();

	//shape tensor
	S00=0.;
        S01=0.;

	//nematic part
	thetaQ = thetaQ_old = PI * (1-2*drand48());
	Q00=cos(2*thetaQ);
	Q01=sin(2*thetaQ);
	nemQ.resize(2);
	nemQ_old.resize(2);
	nemQ = {0., 0.};
	nemQ_old = {0., 0.};
	Q00 = 0.5 * (nemQ[0] * nemQ[0] - nemQ[1] * nemQ[1]);
	Q01 = nemQ[0] * nemQ[1];

	//minor bookkeeping
	Fpassive = std::vector<number> {0.,0.};
        Factive = std::vector<number> {0.,0.};
	area=0;
	offset.resize(2);
	offset[0]=0; offset[1]=0;
	int x,y;
	for(int i=0; i<subSize; i++){
		y=LsubY/2-i/LsubX;
		x=LsubX/2-(i-int(i/LsubX)*LsubX);
		velocityX[i]=0;
		velocityY[i]=0;
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
		y=i/LsubX;
		x=i-y*LsubX;
		ss=0;
		for(int j=-1; j<=1; j++){
			for(int k=-1; k<=1; k++){
				xx=x+k;
				yy=y+j;
				if(x==0 && k==-1)site=-1;
				else if(x==LsubX-1 && k==1)site=-1;
				if(y==0 && j==-1)site=-1;
				else if(y==LsubY-1 && j==1)site=-1;
				site=xx+yy*LsubX;
				neighbors_sub[ss+i*9]=site;
				ss++;
			}
		}
	}
}

void MultiPhaseField::setNeighborsSubSquarePeriodic() {
	int x,y,xx,yy,site,ss;
	for(int i =0; i<subSize; i++) {
		y=i/LsubX;
		x=i-y*LsubX;
		ss=0;
		for(int j=-1; j<=1; j++){
			for(int k=-1; k<=1; k++){
				xx=x+k;
				yy=y+j;
				if(x==0 && k==-1)xx=LsubX-1;
				else if(x==LsubX-1 && k==1)xx=0;
				if(y==0 && j==-1)yy=LsubY-1;
				else if(y==LsubY-1 && j==1)yy=0;
				site=xx+yy*LsubX;
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

	if(x_sub_left>=border && y_sub_bottom>=border){
		int x=CoM[0];
		int y=CoM[1];
		int site=x+y*box->getXsize();
		std::vector<number> corner_old = std::vector<number> { (number)box->getElementX(sub_corner_bottom_left, 0) , (number)box->getElementY(sub_corner_bottom_left, 0) };
		sub_corner_bottom_left=box->getElementX(site, (int)-LsubX/2) + box->getElementY(site, (int)-LsubY/2) * box->getXsize();
	
		std::vector<number> corner_new = std::vector<number> { (number)box->getElementX(sub_corner_bottom_left, 0) , (number)box->getElementY(sub_corner_bottom_left, 0) };
		std::vector<number> displacement = box->min_image(corner_old, corner_new);
		offset[0] =  int(offset[0] + LsubX - displacement[0])%LsubX; offset[1] =  int(offset[1] + LsubY - displacement[1])%LsubY;
	}
	else{
		//create new sizes
		int new_LsubX = LsubX - (2 * x_sub_left) + (2 * border);
		int new_LsubY = LsubY - (2 * y_sub_bottom) + (2 * border);
		int new_subSize = new_LsubX * new_LsubY;
		int new_sub_corner_bottom_left=box->getElementX(sub_corner_bottom_left, -(-(2 * x_sub_left)+(2 * border))/2) + box->getElementY(sub_corner_bottom_left, -(-(2 * y_sub_bottom)+(2 * border))/2) * box->getXsize();
		std::vector<number> new_field_scalar(new_subSize, 0.);
		area=0;
	
		//do mapping between old and new fieldScalars (due to the way the code is written all other variables will be overwritten so a simple resize is enough)
		int start_site_y, start_site_x, dimension_x, dimension_y;
		if(new_LsubY<LsubY){start_site_y=(LsubY-new_LsubY)/2; dimension_y=new_LsubY;}
		else{start_site_y=0; dimension_y=LsubY;}
		if(new_LsubX<LsubX){start_site_x=(LsubX-new_LsubX)/2; dimension_x=new_LsubX;}
		else{start_site_x=0; dimension_x=LsubX;}

		for(int i=start_site_y; i<dimension_y; i++){
			int yy = box->getElementY(sub_corner_bottom_left, i);
			int off_y = i - offset[1];
			if(off_y<0)off_y+=LsubY;
			for(int j=start_site_x; j<dimension_x; j++){
				int xx = box->getElementX(sub_corner_bottom_left, j);

				int dx = (new_sub_corner_bottom_left-int(new_sub_corner_bottom_left/box->getXsize())*box->getXsize()) - xx;
				if(dx>box->getXsize()/2)dx-=box->getXsize();
				else if(dx<-box->getXsize()/2)dx+=box->getXsize();

				int dy = (new_sub_corner_bottom_left/box->getXsize()) - yy;
				if(dy>box->getYsize()/2)dy-=box->getYsize();
				else if(dy<-box->getYsize()/2)dy+=box->getYsize();

				int off_x = j - offset[0];
				if(off_x<0)off_x+=LsubX;
				new_field_scalar[abs(dx)+abs(dy)*new_LsubX]=fieldScalar[off_x+off_y*LsubX];
				area+=fieldScalar[off_x+off_y*LsubX]*fieldScalar[off_x+off_y*LsubX];
			}
		}

		LsubX = new_LsubX;
		LsubY=new_LsubY;
		sub_corner_bottom_left = new_sub_corner_bottom_left;
		offset[0] = 0; offset[1] = 0;	
		subSize = LsubX*LsubY;
		resizing();
		fieldScalar = new_field_scalar;
		setNeighborsSub();
	}
}

void MultiPhaseField::set_positions(int offsetx, int offsety, int corner) {
	sub_corner_bottom_left=corner;
	offset[0]=offsetx;offset[1]=offsety;
}


void MultiPhaseField::set_properties_to_zero() {
	S00=0;
	S01=0;
	x_sub_left = LsubX;
	y_sub_bottom = LsubY;
	area=0;
	sumF=0;
	CoM[0] = 0.; 
	CoM[1] = 0.;
}


void MultiPhaseField::check_borders(int q, int box_size_x, int box_size_y) {

	int dd = (sub_corner_bottom_left-int(sub_corner_bottom_left/box_size_x)*box_size_x) - map_sub_to_box_x[q];
	if(dd>box_size_x/2)dd-=box_size_x;
	else if(dd<-box_size_x/2)dd+=box_size_x;
	if(abs(dd) < x_sub_left)x_sub_left=abs(dd);

	dd = (sub_corner_bottom_left/box_size_x) - map_sub_to_box_y[q];
	if(dd>box_size_y/2)dd-=box_size_y;
	else if(dd<-box_size_y/2)dd+=box_size_y;
	if(abs(dd) < y_sub_bottom)y_sub_bottom=abs(dd);
}

/*transform subdomain sites (patch) into grid sites (box)*/

int MultiPhaseField::GetSubIndex(int site, BaseBox *box){
	return box->getElement(sub_corner_bottom_left, ((site - int(site/LsubX) * LsubX)+offset[0])%LsubX , ((site/LsubX)+offset[1])%LsubY );
}

int MultiPhaseField::GetSubXIndex(int site, BaseBox *box){
	return box->getElementX(sub_corner_bottom_left, ((site - int(site/LsubX) * LsubX)+offset[0])%LsubX);
}

int MultiPhaseField::GetSubYIndex(int site, BaseBox *box){
	return box->getElementY(sub_corner_bottom_left, ((site/LsubX)+offset[1])%LsubY);
}

#include "LEBcMultiPhaseField.h"


LEBcMultiPhaseField::LEBcMultiPhaseField() :  BaseField() {

}


LEBcMultiPhaseField::LEBcMultiPhaseField(int Lx, int Ly) : BaseField() {

	LsubX=Lx;
	LsubY=Ly;
	subSize=Lx*Ly;

}

LEBcMultiPhaseField::~LEBcMultiPhaseField() {

}

void LEBcMultiPhaseField::get_interaction_values(int R) {
        BaseField::get_interaction_values(R);
	init_radius=R/2;
	init_radius2=init_radius*init_radius;
}

void LEBcMultiPhaseField::resizing() {
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

	cos_x_table.resize(LsubX);
	cos_y_table.resize(LsubY);
	sin_x_table.resize(LsubX);
	sin_y_table.resize(LsubY);
	for(int i =0; i<LsubX; i++){
		cos_x_table[i]=cos(2*PI*i/LsubX);
		sin_x_table[i]=sin(2*PI*i/LsubX);
	}
        for(int i =0; i<LsubY; i++){
                cos_y_table[i]=cos(2*PI*i/LsubY);
                sin_y_table[i]=sin(2*PI*i/LsubY);
        }
}

void LEBcMultiPhaseField::init(int Lx, int Ly) {
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
	offset=std::vector<int>(2, 0);

	setNeighborsSub();
}


void LEBcMultiPhaseField::init() {
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
	offset.resize(2);
	offset=std::vector<int>(2, 0);
	area=0;
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

void LEBcMultiPhaseField::setNeighborsSub() {
	//setNeighborsSubSquareDirichlet();
	setNeighborsSubSquarePeriodic();
}


void LEBcMultiPhaseField::setNeighborsSubSquareDirichlet() {
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
				if(site!=-1)site=xx+yy*LsubX;
				neighbors_sub[ss+i*9]=site;
				ss++;
			}
		}
	}
}

void LEBcMultiPhaseField::setNeighborsSubSquarePeriodic() {
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


void LEBcMultiPhaseField::set_positions_initial(BaseBox *box) {
	int x=CoM[0];
	int y=CoM[1];
	int site=x+y*box->getXsize();
	x = box->getElementX(site, (int)-LsubX/2);
	y = box->getElementY(site, (int)-LsubY/2);
	sub_corner_bottom_left = x + y * box->getXsize();
	sub_corner_bottom_left_old = sub_corner_bottom_left;
	unrap_sub_corner_bottom_left_x = x;
	unrap_sub_corner_bottom_left_y = y;
}

void LEBcMultiPhaseField::set_positions(BaseBox *box) {

	if(index>=0){
		int siteCoM = box->getElement(sub_corner_bottom_left, int(CoM[0]), int(CoM[1]));
		CoM[1] = int(siteCoM / box->getXsize()) + (CoM[1]-int(CoM[1]));
		CoM[0] = int(siteCoM - int(siteCoM / box->getXsize()) * box->getXsize()) + (CoM[0]-int(CoM[0]));

		int new_sub_corner_bottom_left = box->getElement(siteCoM, (int)-LsubX/2, (int)-LsubY/2);
		int new_y = new_sub_corner_bottom_left / box->getXsize();
		int new_x = new_sub_corner_bottom_left - new_y * box->getXsize();
		int old_y = sub_corner_bottom_left / box->getXsize();
		int old_x = sub_corner_bottom_left - old_y * box->getXsize();
		std::vector<number> new_field_scalar(subSize, 0.);

		//do mapping between old and new fieldScalars (due to the way the code is written all other variables will be overwritten so a simple resize is enough)
		int yy, xx, site;
		std::vector<number> displacement = std::vector<number> { 0., 0.};

		for(int i=0; i<LsubY; i++){
			for(int j=0; j<LsubX; j++){
				site = j + i * LsubX;

				yy = map_sub_to_box_y[site];
				xx = map_sub_to_box_x[site];
				displacement = box->min_image(std::vector<number> {(number)new_x, (number)new_y}, std::vector<number> {(number)xx, (number)yy});
				//xx = int(LsubX + displacement[0])%LsubX;
				//yy = int(LsubY + displacement[1])%LsubY;
				if(displacement[0]>=0)
					xx = abs(int(displacement[0]));
				else
					xx = int(box->getXsize() + displacement[0])%LsubX;

				if(displacement[1]>=0)
					yy = abs(int(displacement[1]));
				else
					yy = int(box->getYsize() + displacement[1])%LsubY;


				//if(i==LsubY-1 && j==LsubX-1)std::cout<<"mid: "<<xx+yy*LsubX<<" "<<site<<" "<<displacement[0]<<" "<<displacement[1]<<" "<<xx<<" "<<yy<<" "<< map_sub_to_box_x[site]<<" "<<map_sub_to_box_y[site] <<" "<<j<<" "<<i<<std::endl;
				//if(i==LsubY-1 && j==0)std::cout<<"mid: "<<xx+yy*LsubX<<" "<<site<<" "<<displacement[0]<<" "<<displacement[1]<<" "<<xx<<" "<<yy<<" "<< map_sub_to_box_x[site]<<" "<<map_sub_to_box_y[site] <<" "<<j<<" "<<i<<std::endl;
				if(index==14)std::cout<<"mult: "<<xx+yy*LsubX<<" "<<site<<" "<<displacement[0]<<" "<<displacement[1]<<" "<<xx<<" "<<yy<<" "<< map_sub_to_box_x[site]<<" "<<map_sub_to_box_y[site] <<" "<<j<<" "<<i<<std::endl;

				new_field_scalar[xx+yy*LsubX]=fieldScalar[site];
			}
		}

		displacement = box->min_image(std::vector<number> { (number)old_x , (number)old_y } , std::vector<number> { (number)new_x , (number)new_y} );
		unrap_sub_corner_bottom_left_x += displacement[0];
		unrap_sub_corner_bottom_left_y += displacement[1];

		//std::cout<<"here-------"<<new_x<<" "<<new_y<<" "<<old_x<<" "<<old_y<<" "<<displacement[0]<<" "<<displacement[1]<<" "<<new_field_scalar[subSize-1]<<" "<<fieldScalar[subSize-1]<<" "<<new_field_scalar[subSize-1]  <<" "<<CoM[0]<<" "<<CoM[1]<<" "<<unrap_sub_corner_bottom_left_x<<" "<<unrap_sub_corner_bottom_left_y<<" "<<offset[0]<<" "<<offset[1]  <<std::endl;
		//std::cout<<"here-------"<<new_x<<" "<<new_y<<" "<<old_x<<" "<<old_y<<" "<<displacement[0]<<" "<<displacement[1]<<" "<<new_field_scalar[26]<<" "<<fieldScalar[0]<<" "<<new_field_scalar[0]<<" "<<fieldScalar[0]  <<" "<<CoM[0]<<" "<<CoM[1]<<" "<<unrap_sub_corner_bottom_left_x<<" "<<unrap_sub_corner_bottom_left_y<<" "<<offset[0]<<" "<<offset[1]  <<std::endl;

		//exit(911);

		sub_corner_bottom_left_old = sub_corner_bottom_left;
		sub_corner_bottom_left = new_sub_corner_bottom_left;
		fieldScalar = new_field_scalar;
	}
}

void LEBcMultiPhaseField::set_positions(int offsetx, int offsety, int corner, int corner_x, int corner_y, int size_x) {
	unrap_sub_corner_bottom_left_x = corner_x;
	unrap_sub_corner_bottom_left_y = corner_y;
	sub_corner_bottom_left = corner;
	sub_corner_bottom_left_old = sub_corner_bottom_left;
	offset[0]=offsetx;offset[1]=offsety;
}


void LEBcMultiPhaseField::set_properties_to_zero() {
	S00=0;
	S01=0;
	x_sub_left = LsubX;
	y_sub_bottom = LsubY;
	area=0;
	sumF=0;
	CoM[0] = 0.; 
	CoM[1] = 0.;
}


void LEBcMultiPhaseField::check_borders(int q) {
	
	int y = q / LsubX;
	int x = q - y * LsubX;

	x = (x + offset[0])%LsubX;
	y = (y + offset[1])%LsubY;

	if(x < x_sub_left)x_sub_left= x;
	if(y < y_sub_bottom)y_sub_bottom=y;
}

/*transform subdomain sites (patch) into grid sites (box)*/

int LEBcMultiPhaseField::GetSubIndex(int site, BaseBox *box){
	return box->getElement(sub_corner_bottom_left, ((site - int(site/LsubX) * LsubX)+offset[0])%LsubX , ((site/LsubX)+offset[1])%LsubY );
}

int LEBcMultiPhaseField::GetSubXIndex(int site, BaseBox *box){
	return box->getElementX(sub_corner_bottom_left, ((site - int(site/LsubX) * LsubX)+offset[0])%LsubX);
}

int LEBcMultiPhaseField::GetSubYIndex(int site, BaseBox *box){
	return box->getElementY(sub_corner_bottom_left, ((site/LsubX)+offset[1])%LsubY);
}

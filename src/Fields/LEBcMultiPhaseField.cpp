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
	laplacianPhi.resize(subSize);	
	freeEnergy.resize(subSize);
	fieldScalar_old.resize(subSize);
	dfield_old.resize(subSize);
	neighbors_sub.resize(subSize*9);
	velocityX.resize(subSize);
        velocityY.resize(subSize);

	velocityX_correction.resize(LsubY);
	phi_correction.resize(LsubY);
        shear_velocity_sign.resize(subSize);

	map_sub_to_box.resize(subSize);
	map_sub_to_box_x.resize(subSize);
	map_sub_to_box_y.resize(subSize);

	Fpassive_x.resize(subSize);
	Fpassive_y.resize(subSize);
	Factive_x.resize(subSize);
	Factive_y.resize(subSize);

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

void LEBcMultiPhaseField::update_positions(BaseBox *box) {
	int x,y,xx,yy,ss;
	for(int i =0; i<subSize; i++) {
		if(map_sub_to_box_y[i]!=0 && map_sub_to_box_y[i]!=box->getYsize()-1)continue;
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

				if(map_sub_to_box_y[i]==0 && j==-1){
					xx = xx + (int)box->get_shear_displacement();
					while(xx>=LsubX)xx-=LsubX;
				}
				else if(map_sub_to_box_y[i]==box->getYsize()-1 && j==1){
					xx = xx - (int)box->get_shear_displacement();
					while(xx<0)xx+=LsubX;
				}

				neighbors_sub[ss+i*9]=xx+yy*LsubX;
				ss++;
			}
		}
	}


	/*int x, y, xx, yy;
	for(int i =0; i<subSize; i++) {
		y = i/LsubX;
		x = i - y * LsubX;
		if(int(sub_corner_bottom_left/box->getXsize())+y<box->getYsize())continue;

		yy = sub_corner_bottom_left / box->getXsize();
		xx = sub_corner_bottom_left - yy * box->getXsize();
		xx += x;
		if(xx>=box->getXsize())xx-=box->getXsize();
		map_sub_to_box_x[i] = xx - int(box->get_shear_displacement());
		
		if(map_sub_to_box_x[i]<0)map_sub_to_box_x[i]+=box->getXsize();
		map_sub_to_box[i] = map_sub_to_box_x[i] + map_sub_to_box_y[i] * box->getXsize();
		//if(i==subSize-1){std::cout<<"here--------------"<<map_sub_to_box[i]<<" "<<map_sub_to_box_x[i]<<" "<<map_sub_to_box_x[0]<<std::endl;}
	}*/

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

	int choose_condition=0;
	if(border<=0)choose_condition=1;
	else if(x_sub_left>=border && x_sub_left<=2*border)choose_condition=1;
	else if(x_sub_left<=2*border && LsubX==box->getXsize())choose_condition=1;
	else choose_condition=2;

	//std::cout<<"here: "<<choose_condition<<" "<<fieldScalar[734]<<std::endl;

	switch(choose_condition){
	case 0:
		std::cout<<"ERROR FIELD SET POSITION CONDITION!!"<<std::endl;
		exit (911);
		break;
	case 1:
	//if((x_sub_left>=border && x_sub_left<=2*border) || (x_sub_left<=2*border && LsubX==box->getXsize())){
		{
		sub_corner_bottom_left_old = sub_corner_bottom_left;
		std::vector<int> displacement = std::vector<int> {int(CoM[0])-int(LsubX/2) , int(CoM[1])-int(LsubY/2)};
		//int new_sub_corner_bottom_left = box->getElement(sub_corner_bottom_left, displacement[0], displacement[1]);
		int new_sub_corner_bottom_left = box->getElementX(sub_corner_bottom_left, displacement[0]) + box->getElementY(sub_corner_bottom_left, displacement[1]) * box->getXsize();

		//if(index==10)std::cout<<"Before: "<< CoM_old[0]<<" "<<CoM_old[1]<<" "<<CoM[0]<<" "<<CoM[1]<<" "<< LsubX/2<<" "<<LsubY/2<<" "<<(int)CoM[0]-(int)LsubX/2  <<" "<< (int)CoM[1]-(int)LsubY/2 << " "<<new_sub_corner_bottom_left <<" "<< sub_corner_bottom_left<<" "<<sub_corner_bottom_left/box->getXsize()  <<std::endl;

		if(new_sub_corner_bottom_left == sub_corner_bottom_left){
			CoM[0]=CoM_old[0];
			CoM[1]=CoM_old[1];
			return;
		}
		else{
			//int siteCoM = box->getElement(((int)CoM_old[0] + (int)CoM_old[1] * box->getXsize()), (int)CoM[0]-(int)LsubX/2, (int)CoM[1]-(int)LsubY/2);
			int siteCoM = box->getElementX(((int)CoM_old[0] + (int)CoM_old[1] * box->getXsize()) , displacement[0]) + box->getElementY(((int)CoM_old[0] + (int)CoM_old[1] * box->getXsize()) , displacement[1]) * box->getXsize();
			CoM[1] = int(siteCoM / box->getXsize()) + (CoM[1]-int(CoM[1]));
			CoM[0] = int(siteCoM - int(siteCoM / box->getXsize()) * box->getXsize()) + (CoM[0]-int(CoM[0]));
		}

		//if(index==0)std::cout<<"First: "<< CoM_old[0]<<" "<<CoM_old[1]<<" "<<CoM[0]<<" "<<CoM[1]<<" "<< LsubX/2<<" "<<LsubY/2<<" "<<(int)CoM[0]-(int)LsubX/2  <<" "<< (int)CoM[1]-(int)LsubY/2 << " "<<new_sub_corner_bottom_left <<" "<< sub_corner_bottom_left<<" "<<sub_corner_bottom_left/box->getXsize()<<" "<< new_sub_corner_bottom_left - new_sub_corner_bottom_left/box->getXsize() * box->getXsize() <<" "<< sub_corner_bottom_left - sub_corner_bottom_left/box->getXsize() * box->getXsize()<<" "<<displacement[0]<<" "<<displacement[1]  <<std::endl;

		//if(index==36 && sub_corner_bottom_left - sub_corner_bottom_left/box->getXsize() * box->getXsize()==107)exit (911);

		int new_y = new_sub_corner_bottom_left / box->getXsize();
		int new_x = new_sub_corner_bottom_left - new_y * box->getXsize();
		int old_y = sub_corner_bottom_left / box->getXsize();
		int old_x = sub_corner_bottom_left - old_y * box->getXsize();
		std::vector<number> new_field_scalar(subSize, 0.);

		//do mapping between old and new fieldScalars (due to the way the code is written all other variables will be overwritten so a simple resize is enough)
		int yy, xx, site, new_box_site;

		std::vector<number> displacement_corner = box->min_image(std::vector<number> { (number)old_x , (number)old_y } , std::vector<number> { (number)new_x , (number)new_y} );
		unrap_sub_corner_bottom_left_x += displacement_corner[0];
		unrap_sub_corner_bottom_left_y += displacement_corner[1];

		for(int i=0; i<LsubY; i++){
			for(int j=0; j<LsubX; j++){
				site = j + i * LsubX;

				//if(i==LsubY-1 && j==LsubX-1)std::cout<<"mid: "<<xx+yy*LsubX<<" "<<site<<" "<<displacement[0]<<" "<<displacement[1]<<" "<<xx<<" "<<yy<<" "<< map_sub_to_box_x[site]<<" "<<map_sub_to_box_y[site] <<" "<<j<<" "<<i<<std::endl;

				//new_box_site = box->getElement(map_sub_to_box[site], displacement[0], displacement[1]);
				new_box_site = box->getElementX(map_sub_to_box[site], displacement[0]) + box->getElementY(map_sub_to_box[site], displacement[1]) * box->getXsize();

				xx = (j + displacement[0]);
				yy = (i + displacement[1]);
				while(yy<0){yy+=LsubY;}
				while(yy>=LsubY){yy-=LsubY;}
				while(xx<0){xx+=LsubX;}
				while(xx>=LsubX){xx-=LsubX;}

				/*if(box->getElementY(new_box_site, 0)-map_sub_to_box_y[site]>=box->getYsize()/2){
					xx = int(xx + int(box->get_shear_displacement()));
					while(xx>=LsubX)xx-=LsubX;
					shear_velocity_sign[xx+yy*LsubX] = 1;
				}
				else if(box->getElementY(new_box_site, 0)-map_sub_to_box_y[site]<=-box->getYsize()/2){
					xx = int(xx - int(box->get_shear_displacement()));
					while(xx<0)xx+=LsubX;
					shear_velocity_sign[xx+yy*LsubX] = -1;
				}
				else{
					shear_velocity_sign[xx+yy*LsubX] = 0;
				}*/

				//if(box->getElementX(new_box_site, 0) == 47 && box->getElementY(new_box_site, 0) == 0)std::cout<<"maybe "<<site<<" "<<xx+yy*LsubX<<" "<<xx<<" "<<yy<<" "<<  shear_velocity_sign[xx+yy*LsubX]<<" "<<displacement[0]<<" "<<displacement[1]<<" "<< j <<" "<<i <<std::endl;
				//if(box->getElementX(new_box_site, 0) == 58 && box->getElementY(new_box_site, 0) == 99)std::cout<<"maybe "<<site<<" "<<xx+yy*LsubX<<" "<<xx<<" "<<yy<<" "<<  shear_velocity_sign[xx+yy*LsubX]<<" "<<displacement[0]<<" "<<displacement[1]<<" "<< j <<" "<<i <<std::endl;
				//new_field_scalar[xx + yy * LsubX] = fieldScalar[site];
				new_field_scalar[site] = fieldScalar[xx + yy * LsubX];
				map_sub_to_box[site] = new_box_site; 
				map_sub_to_box_x[site] = box->getElementX(new_box_site, 0);
				map_sub_to_box_y[site] = box->getElementY(new_box_site, 0);
			}
		}

		//std::cout<<"here-------"<<new_x<<" "<<new_y<<" "<<old_x<<" "<<old_y<<" "<<displacement[0]<<" "<<displacement[1]<<" "<<new_field_scalar[subSize-1]<<" "<<fieldScalar[subSize-1]<<" "<<new_field_scalar[subSize-1]  <<" "<<CoM[0]<<" "<<CoM[1]<<" "<<unrap_sub_corner_bottom_left_x<<" "<<unrap_sub_corner_bottom_left_y<<" "<<offset[0]<<" "<<offset[1]  <<std::endl;

		sub_corner_bottom_left = new_sub_corner_bottom_left;
		fieldScalar = new_field_scalar;
		setNeighborsSub();
		break;
		}
	case 2:
	//else{
		{
		//create new sizes
		int new_LsubX = LsubX - (2 * x_sub_left) + (2 * border);
		if(LsubX<box->getXsize() && new_LsubX>=box->getXsize())new_LsubX=box->getXsize();
		else if(LsubX==box->getXsize() && new_LsubX>box->getXsize())return;
		if(new_LsubX == LsubX)return;

		int new_subSize = new_LsubX * LsubY;
		int new_sub_corner_bottom_left=box->getElement(sub_corner_bottom_left, -(-(2 * x_sub_left)+(2 * border))/2, 0);
		std::vector<number> new_field_scalar(new_subSize, 0.);
		std::vector<int> new_map_sub_to_box(new_subSize, 0);
		std::vector<int> new_map_sub_to_box_x(new_subSize, 0);
		std::vector<int> new_map_sub_to_box_y(new_subSize, 0);
		area=0;
	
		//do mapping between old and new fieldScalars (due to the way the code is written all other variables will be overwritten so a simple resize is enough)
		int start_site_x = (LsubX - new_LsubX)/2;

		std::vector<number> displacement;
		displacement = box->min_image(std::vector<number> { (number)box->getElementX(sub_corner_bottom_left, 0) , (number)box->getElementY(sub_corner_bottom_left, 0) } , std::vector<number> { (number)box->getElementX(new_sub_corner_bottom_left, 0) , (number)box->getElementY(new_sub_corner_bottom_left, 0)} );
		unrap_sub_corner_bottom_left_x += displacement[0];
		unrap_sub_corner_bottom_left_y += displacement[1];

		//if(index==36)std::cout<<"Second: "<< CoM_old[0]<<" "<<CoM_old[1]<<" "<<CoM[0]<<" "<<CoM[1]<<" "<< LsubX<<" "<< new_LsubX<<" "<< x_sub_left <<" "<< border <<" "<<displacement[0] <<" "<< displacement[1] << " "<<new_sub_corner_bottom_left <<" "<< sub_corner_bottom_left<<" "<< new_sub_corner_bottom_left - new_sub_corner_bottom_left/box->getXsize() * box->getXsize() <<" "<< sub_corner_bottom_left - sub_corner_bottom_left/box->getXsize() * box->getXsize()  <<std::endl;

		//if(index==36 && x_sub_left==12)exit (911);

		int row, site, site_old, new_box_site;
		for(int i=0; i<LsubY; i++){
			row = map_sub_to_box[0+i*LsubX];
			for(int j=0; j<new_LsubX; j++){
				site = j + i * new_LsubX;

				//if(i==LsubY-1 && j==LsubX-1)std::cout<<"mid: "<<xx+yy*LsubX<<" "<<site<<" "<<displacement[0]<<" "<<displacement[1]<<" "<<xx<<" "<<yy<<" "<< map_sub_to_box_x[site]<<" "<<map_sub_to_box_y[site] <<" "<<j<<" "<<i<<std::endl;

				if(new_LsubX>LsubX){
					new_box_site = box->getElement(row, j+start_site_x, 0);
					new_map_sub_to_box[site] = new_box_site; 
					new_map_sub_to_box_x[site] = box->getElementX(new_box_site, 0);
					new_map_sub_to_box_y[site] = box->getElementY(new_box_site, 0);
					if(j+start_site_x>=0 && j+start_site_x<LsubX){
						site_old = (j+start_site_x) + i * LsubX;
						new_field_scalar[site] = fieldScalar[site_old];
						area += new_field_scalar[site] * new_field_scalar[site];
					}
				}
				else{
					new_box_site = box->getElement(row, j+start_site_x, 0);
					new_map_sub_to_box[site] = new_box_site; 
					new_map_sub_to_box_x[site] = box->getElementX(new_box_site, 0);
					new_map_sub_to_box_y[site] = box->getElementY(new_box_site, 0);
					site_old = (j+start_site_x) + i * LsubX;
					new_field_scalar[site] = fieldScalar[site_old];
					area += new_field_scalar[site] * new_field_scalar[site];
				}
			}
		}

		LsubX = new_LsubX;
		sub_corner_bottom_left_old = sub_corner_bottom_left;
		sub_corner_bottom_left = new_sub_corner_bottom_left;
		offset[0] = 0; offset[1] = 0;	
		subSize = LsubX*LsubY;
		resizing();
		map_sub_to_box = new_map_sub_to_box;
		map_sub_to_box_x = new_map_sub_to_box_x;
		map_sub_to_box_y = new_map_sub_to_box_y;
		fieldScalar = new_field_scalar;
		setNeighborsSub();
		CoM[0]=CoM_old[0];
		CoM[1]=CoM_old[1];
		//int siteCoM = box->getElement(((int)CoM_old[0] + (int)CoM_old[1] * box->getXsize()), -(int)displacement[0], -(int)displacement[1]);
		//CoM[1] = int(siteCoM / box->getXsize()) + (CoM[1]-int(CoM[1]));
		//CoM[0] = int(siteCoM - int(siteCoM / box->getXsize()) * box->getXsize()) + (CoM[0]-int(CoM[0]));
		break;
		}
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

	if(x < x_sub_left)x_sub_left = x;
	if(LsubX - x < x_sub_left)x_sub_left = LsubX - x;
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

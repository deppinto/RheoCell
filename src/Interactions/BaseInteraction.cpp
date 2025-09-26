#include "BaseInteraction.h"

BaseInteraction::BaseInteraction() {
	energy_threshold = (number) 100.f;
	is_infinite = false;
	box = NULL;
	R=8;
	rcut = 7; //care initial cell size
	sqr_rcut = SQR(rcut);
	wall_slip = 0.5;
}

BaseInteraction::~BaseInteraction() {

}


void BaseInteraction::get_settings(input_file &inp) {
	getInputString(&inp, "topology", topology_filename, 1);
	getInputNumber(&inp, "energy_threshold", &energy_threshold, 0);
	getInputInt(&inp, "omp_thread_num", &omp_thread_num, 0);
	getInputInt(&inp, "R", &R, 1);
	getInputNumber(&inp, "wall_slip", &wall_slip, 0);
}


int BaseInteraction::get_N_from_topology() {
        std::ifstream topology;
        topology.open(topology_filename, std::ios::in);
        if(!topology.good()) throw RCexception("Can't read topology file '%s'. Aborting", topology_filename);
        int ret;
        topology >> ret;
        topology.close();
        return ret;
}

void BaseInteraction::read_topology(std::vector<BaseField *> &fields) {
        allocate_fields(fields);
        int idx = 0;
        for(auto p : fields) {
                p->index = idx;
                idx++;
        }
}

void BaseInteraction::apply_changes_after_equilibration(){

}

void BaseInteraction::begin_energy_computation() {

}

void BaseInteraction::begin_energy_computation(std::vector<BaseField *> &fields) {

}

void BaseInteraction::resetSums(int k){

}

void BaseInteraction::updateFieldProperties(BaseField *p, int q, int k){
	p->total_force_x = 0.;
	p->total_force_y = 0.;
	p->total_force_repulsion_x = 0.;
	p->total_force_repulsion_y = 0.;
	p->freeEnergy[q] = p->set_F_ext(q, box->getWalls(k), box->getLaplacianWalls(k));
}

void BaseInteraction::update_sub_to_box_map(BaseField *p, int q, int sub_site, int sub_site_x, int sub_site_y){
	p->map_sub_to_box[q]=sub_site;
	p->map_sub_to_box_x[q]=sub_site_x;
	p->map_sub_to_box_y[q]=sub_site_y;
}

number BaseInteraction::get_system_energy(std::vector<BaseField *> &fields) {

	/*number energy = 0.;
	for(auto &field : fields) {
		BaseField *p = field;
		for(int q=0; q<p->subSize; q++)energy += (number) field_interaction(p,q,false);
		if(get_is_infinite()) {
			return energy;
		}
	}

	return (number) energy;*/

	return U + K;
}


bool BaseInteraction::generate_random_configuration_overlap(BaseField *p, BaseField *q) {
	computed_r = box->sqr_min_image_distance(p->CoM, q->CoM);
	//printf("values: %f, %f\n", computed_r,rcut);
	if(computed_r >= 2.*sqr_rcut)return false;
	else return true;
}


void BaseInteraction::generate_random_configuration(std::vector<BaseField *> &fields) {

	for(auto p : fields) {
		p->CoM = std::vector<number> {drand48() * box->getXsize(), drand48() * box->getYsize()};
	}

	int N = fields.size();
	number totalNodes=box->getXsize()*box->getYsize();
	rcut= sqrt( (totalNodes/N)/PI );
	rcut=(double)R-3;
	sqr_rcut = SQR(rcut);
	std::cout<<"This is rcut: " << rcut <<std::endl;
	for(int i = 0; i < N; i++) {
		BaseField *p = fields[i];

		bool inserted = false;
		do {
			p->CoM = std::vector<number> {drand48() * box->getXsize(), drand48() * box->getYsize()};
			p->set_positions_initial(box);

			inserted = true;
			for(int n = 0; n < i; n++) {
				BaseField *q = fields[n];
				// particles with an index larger than p->index have not been inserted yet
				if(generate_random_configuration_overlap(p, q)) inserted = false;
			}

			//we take into account the external potential
			number ext_value=0.;
			int start=int(p->CoM[0])+int(p->CoM[1])*box->getXsize();
			int startX=box->getElementX(start, (int)rcut*(-1));
			int startY=box->getElementY(start, (int)rcut*(-1));
			start=startX+startY*box->getXsize();
			for(int y=0; y<2*rcut; y++){
				for(int x=0; x<2*rcut; x++){
					int f = box->getElementX(start, x) + box->getElementY(start, y) * box->getXsize();
					p->set_ext_potential(f, box->getWalls(f));	
					ext_value += p->ext_potential;
				}	
			}

			//std::cout<<"generator: "<<start<<" "<<ext_value<<std::endl;
			if(std::isnan(ext_value) || ext_value > 20) {
				inserted = false;
			}

		} while(!inserted);

		if(i > 0 && N > 5 && i % (N / 4) == 0) printf("Inserted %d%% of the particles (%d/%d)\n", i*100/N, i, N);
	}
	printf("Inserted %d%% of the particles (%d/%d)\n", 100, N, N);
}


void BaseInteraction::generate_lattice_configuration(std::vector<BaseField *> &fields, int n_rows, int n_columns) {


	// Choose your oblique lattice parameters
	number a = 2*R - 2;         // spacing along columns
	number b = 2*R - 2;         // spacing along rows (can be different from a)
	number theta = 60.0 * M_PI/180.0;  // angle between a and b (radians)

	// Lattice basis vectors (x right, y up; use -sin if you want rows to go downward)
	number ax = a;              number ay = 0;
	number bx = b * cos(theta); number by = b * sin(theta);

	// Center the lattice in the box
	int N = fields.size();
	number cx = box->getXsize() / 2.0;
	number cy = box->getYsize() / 2.0;

	// If you want the whole n_columns x n_rows block centered:
	number start_x = cx - 0.5 * ((n_columns - 1) * ax + (n_rows - 1) * bx);
	number start_y = cy - 0.5 * ((n_columns - 1) * ay + (n_rows - 1) * by);


	// Honeycomb hex: add a second point per site
	bool make_honeycomb = true; // set false for triangular
	std::vector<std::vector<number>> basis;
	basis.push_back({0.0, 0.0});                         // A sublattice
	if (make_honeycomb) basis.push_back({(ax+bx)/3.0, (ay+by)/3.0}); // B sublattice

	int placed = 0;
	for (int r = 0; r < n_rows && placed < N; ++r) {
		for (int c = 0; c < n_columns && placed < N; ++c) {
			number x0 = start_x + c*ax + r*bx;
			number y0 = start_y + c*ay + r*by;
			
			for (auto &d : basis) {
				if (placed >= N) break;
				BaseField* p = fields[placed++];
				p->CoM = std::vector<number>{ x0 + d[0], y0 + d[1] };
				p->set_positions_initial(box);
				if(placed > 0 && N > 5 && placed % (N / 4) == 0) printf("Inserted %d%% of the particles (%d/%d)\n", placed*100/N, placed, N);
			}
		}
	}
	printf("Inserted %d%% of the particles (%d/%d)\n", 100, placed, N);


	/*int N = fields.size();
	//number totalNodes = box->getXsize() * box->getYsize();
	number lattice_const = 2 * R - 2;
	number lattice_factor = sqrt(3.) / 2.; 
	//int n_rows = 20;
	//int n_columns = 5;
	number start_x = box->getXsize()/2. - (n_columns/2.) * lattice_const;
	number start_y = box->getYsize()/2. + (n_rows/2.) * lattice_const * lattice_factor;
	int lattice_row = 0;
	int lattice_column = 0;
	for(int i = 0; i < N; i++) {
		BaseField *p = fields[i];

		bool inserted = false;
		do {
			if(lattice_row % 2 == 0){p->CoM = std::vector<number> {start_x + lattice_column * lattice_const, start_y - (lattice_row+0.5) * lattice_const * lattice_factor};}
			else {p->CoM = std::vector<number> {start_x + lattice_column * lattice_const + lattice_const/2, start_y - (lattice_row+0.5) * lattice_const * lattice_factor};}
			p->set_positions_initial(box);

			lattice_column += 1;
			if(i > 0 && (i+1) % n_columns == 0){
				lattice_column = 0;
				lattice_row += 1;
			}

			inserted = true;
		} while(!inserted);

		if(i > 0 && N > 5 && i % (N / 4) == 0) printf("Inserted %d%% of the particles (%d/%d)\n", i*100/N, i, N);
	}
	printf("Inserted %d%% of the particles (%d/%d)\n", 100, N, N);*/
}


void BaseInteraction::generate_cluster_configuration(std::vector<BaseField *> &fields) {

	for(auto p : fields) {
		p->CoM = std::vector<number> {drand48() * box->getXsize(), drand48() * box->getYsize()};
	}

	int N = fields.size();
	number totalNodes=box->getXsize()*box->getYsize();
	number cluster_radius = R * R * N;
	std::vector<number> reference_point = std::vector<number> {box->getXsize()/2., box->getYsize()/2.};
	rcut= sqrt( (totalNodes/N)/PI );
	rcut=(double)R-1;
	sqr_rcut = SQR(rcut);
	std::cout<<"This is rcut: " << rcut <<std::endl;
	for(int i = 0; i < N; i++) {
		BaseField *p = fields[i];

		bool inserted = false;
		do {
			p->CoM = std::vector<number> {drand48() * box->getXsize(), drand48() * box->getYsize()};
			p->set_positions_initial(box);

			inserted = true;


			if(box->sqr_min_image_distance(p->CoM, reference_point) > cluster_radius)inserted = false;


			for(int n = 0; n < i; n++) {
				BaseField *q = fields[n];
				// particles with an index larger than p->index have not been inserted yet
				if(generate_random_configuration_overlap(p, q)) inserted = false;
			}

			//we take into account the external potential
			number ext_value=0.;
			int start=int(p->CoM[0])+int(p->CoM[1])*box->getXsize();
			int startX=box->getElementX(start, (int)rcut*(-1));
			int startY=box->getElementY(start, (int)rcut*(-1));
			start=startX+startY*box->getXsize();
			for(int y=0; y<2*rcut; y++){
				for(int x=0; x<2*rcut; x++){
					int f = box->getElementX(start, x) + box->getElementY(start, y) * box->getXsize();
					p->set_ext_potential(f, box->getWalls(f));	
					ext_value += p->ext_potential;
				}	
			}

			//std::cout<<"generator: "<<start<<" "<<ext_value<<std::endl;
			if(std::isnan(ext_value) || ext_value > 20) {
				inserted = false;
			}

		} while(!inserted);

		if(i > 0 && N > 5 && i % (N / 4) == 0) printf("Inserted %d%% of the particles (%d/%d)\n", i*100/N, i, N);
	}
	printf("Inserted %d%% of the particles (%d/%d)\n", 100, N, N);
}

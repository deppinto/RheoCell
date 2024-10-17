#include "BaseInteraction.h"

BaseInteraction::BaseInteraction() {
	energy_threshold = (number) 100.f;
	is_infinite = false;
	box = NULL;
	rcut = 5*5; //care initial cell size
	sqr_rcut = SQR(rcut);
}

BaseInteraction::~BaseInteraction() {

}


void BaseInteraction::get_settings(input_file &inp) {
	getInputString(&inp, "topology", topology_filename, 1);
	getInputNumber(&inp, "energy_threshold", &energy_threshold, 0);
}


int BaseInteraction::get_N_from_topology() {
        std::ifstream topology;
        topology.open(topology_filename, std::ios::in);
        if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", topology_filename);
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

void BaseInteraction::begin_energy_computation() {

}

void BaseInteraction::begin_energy_computation(std::vector<BaseField *> &fields) {

}

void BaseInteraction::resetSums(int k){

}

void BaseInteraction::updateFieldProperties(BaseField *p, int q){

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
	sqr_rcut = SQR(rcut);
	for(int i = 0; i < N; i++) {
		BaseField *p = fields[i];

		bool inserted = false;
		do {
			p->CoM = std::vector<number> {drand48() * box->getXsize(), drand48() * box->getYsize()};

			p->set_positions_initial(box);
			p->set_ext_potential(0, box);

			inserted = true;


			for(int n = 0; n < i; n++) {
				BaseField *q = fields[n];
				// particles with an index larger than p->index have not been inserted yet
				if(generate_random_configuration_overlap(p, q)) inserted = false;
			}

			// we take into account the external potential
			//number boltzmann_factor = exp(-p->ext_potential / temperature);
			//if(std::isnan(p->ext_potential) || drand48() > boltzmann_factor) {
			//	inserted = false;
			//}

		} while(!inserted);

		if(i > 0 && N > 5 && i % (N / 2) == 0) printf("Inserted %d%% of the particles (%d/%d)\n", i*100/N, i, N);
	}
}

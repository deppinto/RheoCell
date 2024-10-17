#include "BackendFactory.h"

#include "FD_CPUBackend.h"
#ifndef NOCUDA
#include "../CUDA/Backends/FD_CUDABackend.h"
#ifndef CUDA_DOUBLE_PRECISION
#include "../CUDA/Backends/FD_CUDAMixedBackend.h"
#endif
#endif

std::shared_ptr<SimBackend> BackendFactory::make_backend(input_file &inp) {
	std::string backend_opt, backend_prec, sim_type;
	getInputString(&inp, "backend", backend_opt, 1);

	int precision_state = getInputString(&inp, "backend_precision", backend_prec, 0);

	if(precision_state == KEY_FOUND && backend_opt == "CPU") {
		OX_LOG(Logger::LOG_NOTHING, "");
		OX_LOG(Logger::LOG_WARNING, "The 'backend_precision' option cannot be set by input file when running on CPU\n");
	}

	if(getInputString(&inp, "sim_type", sim_type, 0) == KEY_NOT_FOUND) {
		OX_LOG(Logger::LOG_INFO, "Simulation type not specified, using FD");
		sim_type = "FD";
	}
	else {
		OX_LOG(Logger::LOG_INFO, "Simulation type: %s", sim_type.c_str());
	}


	SimBackend *new_backend = nullptr;
	if(sim_type == "FD") {
		if(backend_opt == "CPU") {
			new_backend = new FD_CPUBackend();
		}
#ifndef NOCUDA
		else if(backend_opt == "CUDA") {
#ifndef CUDA_DOUBLE_PRECISION
			if(precision_state == KEY_NOT_FOUND) {
				backend_prec = "mixed";
			}
			if(backend_prec == "mixed") {
				new_backend = new CUDAMixedBackend();
			}
			else if(backend_prec == "float") {
				new_backend = new FD_CUDABackend();
			}
			else {
				throw oxDNAException("Backend precision '%s' is not allowed, as the code has been compiled with 'float' and 'mixed' support only", backend_prec.c_str());
			}
#else
			if(precision_state == KEY_NOT_FOUND) {
				backend_prec = "double";
			}
			if(backend_prec == "double") {
				new_backend = new FD_CUDABackend();
			}
			else {
				throw oxDNAException("Backend precision '%s' is not allowed, as the code has been compiled with 'double' support only", backend_prec.c_str());
			}
#endif
			OX_LOG(Logger::LOG_INFO, "CUDA backend precision: %s", backend_prec.c_str());
		}
#endif
		else {
			throw oxDNAException("Backend '%s' not supported", backend_opt.c_str());
		}
	}
	else {
		throw oxDNAException("Simulation type '%s' not supported", sim_type.c_str());
	}

	return std::shared_ptr<SimBackend>(new_backend);
}

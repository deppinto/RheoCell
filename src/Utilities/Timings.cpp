#include "Timings.h"
#include "RCexception.h"
#include "Logger.h"

#include <algorithm>

#ifdef NOCUDA
#define SYNCHRONIZE()
#else
#include <cuda_runtime_api.h>
#define SYNCHRONIZE() cudaDeviceSynchronize()
#endif

#ifdef MOSIX
#define OXDNA_CLOCK() 0
#else
#define OXDNA_CLOCK() clock()
#endif

Timer::Timer(bool ssync) {
	time = (clock_t) 0;
	last = (clock_t) 0;
	active = false;
	sync = ssync;
	desc = std::string("Uninitialized timer");
}

Timer::Timer(std::string arg, bool ssync) :
	Timer(ssync) {
	desc = std::string(arg);
}

Timer::~Timer() {
	OX_DEBUG("Timer with desc %s deleted", desc.c_str());
}

void Timer::resume() {
	if(active) {
		throw RCexception("resuming already active timer %s", desc.c_str());
	}
	last = OXDNA_CLOCK();
	active = true;
}

void Timer::pause() {
	if(sync) {
		SYNCHRONIZE();
	}
	if(!active) {
		throw RCexception("pausing resuming already inactive timer %s", desc.c_str());
	}
	time += (OXDNA_CLOCK() - last);
	active = false;
}

// this should work regardless of the timers being active
long long int Timer::get_time() {
	if(active) {
		return (long long int) (time + (OXDNA_CLOCK() - last));
	}
	else {
		return (long long int) time;
	}
}

/***************** END OF TIMER CLASS *********************/

// singleton
TimingManager *TimingManager::timingManager = nullptr;

// time manager class
TimingManager::TimingManager() {

}

TimingManager::~TimingManager() {

}

void TimingManager::enable_sync() {
	sync = true;
	for(auto t: timers) {
		t->set_sync(sync);
	}
}

void TimingManager::disable_sync() {
	sync = false;
	for(auto t: timers) {
		t->set_sync(sync);
	}
}

void TimingManager::init() {
	if(timingManager != nullptr) {
		throw RCexception("initializing an already initialized TimingManager");
	}
	timingManager = new TimingManager();
}

void TimingManager::clear() {
	if(timingManager != nullptr) {
		delete timingManager;
	}
	timingManager = nullptr;
}

TimingManager *TimingManager::instance() {
	if(timingManager == nullptr) {
		throw RCexception("accessing uninitialized TimingManager");
	}
	return timingManager;
}

void TimingManager::add_timer(TimerPtr arg) {
	timers.push_back(arg);
	parents.insert(std::make_pair(arg, nullptr));
	desc_map.insert(std::make_pair(arg->get_desc(), arg));
}

TimerPtr TimingManager::new_timer(std::string desc) {
	if(desc_map.count(desc) != 0) {
		throw RCexception("timer %s already used! Aborting", desc.c_str());
	}

	TimerPtr timer = std::make_shared<Timer>(desc, sync);

	timers.push_back(timer);
	parents[timer] = nullptr;
	desc_map[desc] = timer;

	OX_DEBUG("Adding new timer with description %s and no parent", desc.c_str());

	return timer;
}

TimerPtr TimingManager::new_timer(std::string desc, std::string parent_desc) {
	if(desc_map.count(desc) != 0) {
		throw RCexception("timer %s already used! Aborting", desc.c_str());
	}
	if(desc_map.count(parent_desc) == 0) {
		throw RCexception("Cannot add timer %s because parent timer %s does not exist", desc.c_str(), parent_desc.c_str());
	}

	TimerPtr timer = std::make_shared<Timer>(desc, sync);

	timers.push_back(timer);
	parents[timer] = get_timer_by_desc(parent_desc);
	desc_map[desc] = timer;

	OX_DEBUG("Adding new timer with description %s and parent %s", desc.c_str(), parent_desc.c_str());

	return timer;
}

void TimingManager::add_timer(TimerPtr arg, std::string parent_desc) {
	std::string my_parent_desc;
	TimerPtr my_parent_ptr;
	if(desc_map.count(parent_desc) > 0) {
		my_parent_desc = std::string(parent_desc);
		my_parent_ptr = desc_map[parent_desc];
	}
	else {
		OX_LOG(Logger::LOG_WARNING, "Trying to add timer \"%s\" with an unknown parent \"%s\". Setting parent to \"None\"", arg->get_desc().c_str(), parent_desc.c_str());
		my_parent_desc = std::string("None");
		my_parent_ptr = nullptr;
	}

	timers.push_back(arg);
	parents.insert(std::make_pair(arg, my_parent_ptr));
	desc_map.insert(std::make_pair(arg->get_desc(), arg));
}

void TimingManager::print(long long int total_steps) {
	// times (including children) 
	std::map<TimerPtr, long long int> totaltimes;
	for(unsigned int i = 0; i < timers.size(); i++)
		totaltimes[timers[i]] = timers[i]->get_time();

	// times in children 
	std::map<TimerPtr, long long int> sum_of_children;
	for(unsigned int i = 0; i < timers.size(); i++) {
		sum_of_children[timers[i]] = 0;
	}
	for(unsigned int i = 0; i < timers.size(); i++) {
		TimerPtr t = timers[i];
		TimerPtr p = parents[t];
		while(p != nullptr) {
			sum_of_children[p] += totaltimes[t];
			p = parents[p];
		}
	}

	// own time (not in children)
	std::map<TimerPtr, long long int> own_time;
	for(unsigned int i = 0; i < timers.size(); i++) {
		TimerPtr t = timers[i];
		own_time[t] = totaltimes[t] - sum_of_children[t];
	}

	// mylist will be ordered as a tree
	std::vector<std::string> mylist;
	while(mylist.size() < timers.size()) {
		for(unsigned int i = 0; i < timers.size(); i++) {
			TimerPtr t = timers[i];
			TimerPtr p = parents[t];

			if(p == nullptr) {
				mylist.push_back(t->get_desc());
			}
			else {
				// troviamo il nome del parente
				std::vector<std::string>::iterator it = std::find(mylist.begin(), mylist.end(), p->get_desc());
				if(it != mylist.end()) {
					it++;
					mylist.insert(it, t->get_desc());
				}
			}
		}
	}

	// now the list is ordered in the order we want to print it
	double tot = (double) get_timer_by_desc("SimBackend")->get_time() / CPSF;
	if(tot < 1e-10) {
		OX_LOG(Logger::LOG_INFO, "No timings available (either oxDNA was compiled with MOSIX=1 or no simulation steps were performed)");
		return;
	}

	OX_LOG(Logger::LOG_NOTHING, "");
	OX_LOG(Logger::LOG_INFO, "Total Running Time: %g s, per step: %g ms", tot, tot / total_steps * 1000.);
	OX_LOG(Logger::LOG_INFO, "Timings, in seconds, by Timer (total, own, spent in children)");
	for(unsigned int i = 0; i < mylist.size(); i++) {
		char mystr[512] = "";
		TimerPtr t = get_timer_by_desc(mylist[i]);
		TimerPtr p = parents[t];
		int generations = 0;
		while(p != nullptr) {
			generations++;
			p = parents[p];
		}
		for(int j = 0; j < generations; j++) {
			strcat(mystr, "***");
		}
		strcat(mystr, "> ");
		strcat(mystr, t->get_desc().c_str());
		OX_LOG(Logger::LOG_NOTHING, "%-30s %12.3lf (%5.1lf\%) %12.3lf (%5.1f\%) %12.3lf (%5.1f\%)",
		(char *) mystr,
		totaltimes[t] / CPSF, totaltimes[t] / CPSF / tot * 100.,
		own_time[t] / CPSF, own_time[t] / CPSF / tot * 100.,
		sum_of_children[t] / CPSF, sum_of_children[t] / CPSF / tot * 100.);
	}
	OX_LOG(Logger::LOG_NOTHING, "");

	return;
}


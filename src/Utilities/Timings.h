#ifndef TIMINGS_H_
#define TIMINGS_H_

#include <ctime>
#include <vector>
#include <map>
#include <string>
#include <memory>

#ifndef CPS
#define CPS (CLOCKS_PER_SEC)
#define CPSF ((double)(CLOCKS_PER_SEC))
#endif

class Timer {
private:
	/// string describing the timer
	std::string desc;

	/// time elapsed between pauses
	clock_t time;

	/// last time that the timer was started
	clock_t last;

	/// whether the timer is running or paused
	bool active;

	/// whether we should call SYNCHRONIZE() or not
	bool sync;

public:
	Timer(bool sync);
	Timer(std::string desc, bool ssync);
	~Timer();

	/// resumes (or starts) the timer
	void resume();

	/// pauses the timer
	void pause();

	/// resets the timer
	void reset() {
		time = (clock_t) 0;
	}

	/// returns the time passed as a llint; works for active and inactive timers
	long long int get_time();

	/// returns the time passed, converted in seconds, as a double
	double get_seconds() {
		return (double) time / (double) CPS;
	}

	/// returns the description of the clock
	std::string get_desc() {
		return desc;
	}

	/// returns whether the timer is active
	bool is_active() {
		return active;
	}

	void set_sync(bool state) {
		sync = state;
	}
};

using TimerPtr = std::shared_ptr<Timer>;

/// new attempt at a singleton to be able to add from anywhere within the code a timer and have it do what expected
class TimingManager {
private:
	std::vector<TimerPtr> timers;
	std::map<TimerPtr, TimerPtr> parents;
	std::map<std::string, TimerPtr> desc_map;
	bool sync = false;

	static TimingManager *timingManager;

	/**
	 * @brief Default constructor. It is kept private to enforce the singleton pattern.
	 *
	 * @return
	 */
	TimingManager();

public:
	TimingManager(TimingManager const&) = delete;

	/// creates a new orphan timer
	TimerPtr new_timer(std::string desc);

	/// creates a new orphan timer
	TimerPtr new_timer(std::string desc, std::string parent_desc);

	/// adds the timer to the manager, using the description to set the parent
	void add_timer(TimerPtr arg, std::string parent_desc);

	/// adds the timer to the manager, setting no parent
	void add_timer(TimerPtr arg);

	/// return the Timer pointer associated to a given description
	TimerPtr get_timer_by_desc(std::string desc) {
		return desc_map[desc];
	}

	void enable_sync();
	void disable_sync();

	/// singleton
	static TimingManager *instance();

	/// init function
	static void init();

	/// clear function
	static void clear();

	/// prints 
	void print(long long int total_steps);

	virtual ~TimingManager();
};
#endif


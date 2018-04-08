#include <queue>
#include <vector>
#include "vnEmbed.h"

using namespace std;
#define EVENT_ARRIVE  0
#define EVENT_DEPART  1
#define EVENT_LOG   3
#define EVENT_END_SIM 2
#define EVENT_PREDICT_PRICE 4
#define EVENT_MIGRATION 5
#define TIME_INTERVAL 100
#define EVENT_MIGCOM 6
class Event;
class Simulator;

class Event {
public:
	int type, index;
	int time;
	int from, to;
	int reqindex;
	Event(int _type, int _time, int _index, int _from, int _to, int _reqindex);
	virtual ~Event();
};

class EventComparer {
public:
	bool operator()(const Event& e1, const Event& e2) const {
		if (e1.time > e2.time)
			return 1;
		else
			return 0;
	}
};

class Simulator {
public:
	priority_queue<Event, vector<Event>, EventComparer> PQ;

	const Event &top();
	void pop();
	bool empty();
	void push(const Event &ev);

	Simulator();
	Simulator(vector<Request> &req);
	virtual ~Simulator();
};

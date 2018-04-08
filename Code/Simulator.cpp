
#include"stdafx.h"
#include "Simulator.h"
using namespace std;
Simulator::Simulator()
{}
Simulator::Simulator(vector<Request> &req){
	for (int i = 0; i<req.size(); i++)
	{
		this->push(Event(EVENT_ARRIVE, req[i].time, i, -1, -1, -1));
	}

	for (int i = 1; i < 1000; i++)
		this->push(Event(EVENT_LOG, i * 100, 0, -1, -1, -1));

	for (int i = 1; i < 1000; i++)
		this->push(Event(EVENT_MIGRATION, i * 100, 0, -1, -1, -1));

	this->push(Event(EVENT_END_SIM, 40000, 0, -1, -1, -1));
}

Simulator::~Simulator() {

}

const Event &Simulator::top() {
	return PQ.top();
}

void Simulator::pop() {
	PQ.pop();
}

bool Simulator::empty() {
	return PQ.empty();
}

void Simulator::push(const Event &ev)
{
	this->PQ.push(ev);
}
Event::Event(int _type, int _time, int _index, int _from, int _to, int _reqindex)
	:type(_type), time(_time), index(_index), from(_from), to(_to), reqindex(_reqindex)  {

}

Event::~Event() {

}

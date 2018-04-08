#include <string.h>
#include <iostream>
#include <ctime>
#include <fstream>
#include <math.h>
#include <string>
#include <algorithm> 
#include <stdio.h>
#include <stdlib.h>
#include <vector>

using namespace std;

#define MAX_PRICE_NUM 50000/60+1
#define MAX_REQ_LINKS 100
#define MAX_REQ_NODES 100
#define MAX_REQ_PER_NODE 1000
#define NOTEXIST -999999
#define INTER_NODE_CPU 0.01
#define BODY_NUM 100
#define DELTA 20
#define Init_NUM 6
#define P1 0.1
#define P2 0.2
#define P3 0.7
#define GROUP 5
#define rest_group -1
class Link;
class Node;
class SubInfo;
class Request;
class SubstrateNetwork;
class CurrentNode;
class CurrentLink;
class CurrentSn;
class Antibody;
class virNode;
class phyNode;

class virNode{
public:
	int req_index;
	int vir_index;
	int jindex;
	int snode_index;
	float cpu;
	float memory;
	float totalcpu;
	vector<int> visitY;
	bool devided;
	int devided_group;
};
class phyNode{
public:
	static const int NODE_ON = 1;
	static const int NODE_OFF = 0;
public:
	int physicalNodeIndex;
	int req_count;//the number of requests this current node is serving
	int req[MAX_REQ_PER_NODE];//vector that notes down each req NO.
	int vnode[MAX_REQ_PER_NODE];//notes down the set of virtual nodes deployed on this node
	float cpu[MAX_REQ_PER_NODE];//the set of cpus allocated to each virtual node
	float memory[MAX_REQ_PER_NODE];
	float rest_cpu;//the rest cpu on this node
	float rest_memory;
	int orireq_count;
	int virqueueindex[MAX_REQ_PER_NODE];
	int hot;
	int state;//indicate the operation condition of this node////////////////////////////////////////
	float energy;
	bool devided = false;

};







class Link{
public:
	int from, to;
	float bw;
	int devided_group;
	Link(int _from = -1, int _to = -1, float _bw = 0);
};

class Node{
public:
	int x;
	int y;
	float cpu;
	float memory;
};

class SubInfo{
public:
	int x;
	int y;

	int nodes;
	int links;

	float price[MAX_PRICE_NUM];
};


class SubstrateNetwork {
public:
	SubInfo subInfo;
	vector<float> cpu;
	vector<float> memory;
	double totalCPU;
	double totalBW;
	vector<Link> link;
	vector<vector<int> > link_matrix;
	vector<vector<int> > matrix;
	SubstrateNetwork();
	int GetLinkIndex(int from, int to) const;

	float cost;

	struct shortest_path *spath;

};


class Request{
public:
	int req_index, duration;//***   
	int nodes;
	int links;

	vector<Node> node;
	vector<Link> link;

	int time;
	int state;
	int split;//£¿
	int maxD;//?
	vector<int> node_mapping_method;
	vector<vector<int> > link_mapping_method;
	vector<int> temp_node_mapping_method;
	vector<vector<int> > temp_link_mapping_method;
	double GetRevenue();

};


struct shortest_path {
	int length;
	int next;
	float min_bw;
};


class CurrentNode{
public:
	static const int NODE_ON = 1;
	static const int NODE_OFF = 0;
public:
	int req_count;//the number of requests this current node is serving
	int req[MAX_REQ_PER_NODE];//vector that notes down each req NO.
	int vnode[MAX_REQ_PER_NODE];//notes down the set of virtual nodes deployed on this node
	float cpu[MAX_REQ_PER_NODE];//the set of cpus allocated to each virtual node
	float memory[MAX_REQ_PER_NODE];
	float rest_cpu;//the rest cpu on this node
	float rest_memory;
	int interNodeConstrain = 0;
	int hot = 0;
	int timeduration = 0;
	int state = 0;//indicate the operation condition of this node
	float energy;

	CurrentNode(float _rest_cpu = 0, float _rest_memory = 0);
};

class CurrentLink{
public:
	int link_count;
	int req[MAX_REQ_LINKS];
	int vlink[MAX_REQ_LINKS];
	float bw[MAX_REQ_LINKS];
	float rest_bw;
	float pre_used;
	CurrentLink(float _rest_bw = 0);

	void printUsing() const;
};
class CurrentSn{
public:
	int time_recent;//last time this sn is updated
	vector<CurrentNode> currentNode;
	vector<CurrentLink> currentLink;

	CurrentSn(SubstrateNetwork &sn);
};

//a specific solution
class PSO
{
public:
	vector<int> node_mapping_method;
	vector<int> pbest;
	vector<int> speed;

	float energy;//indicate the energy(calculated by the energy cost) of current Antibody
	float resource;//indicate the resource cost(calculated by adding all the bw used together)
	float best_resource;
	float fitness;
	float prob;
	double x;
	double y;
	double comp;
	int position;
	PSO();
	void set(const PSO& ant);
	void init();
	void print();
	PSO(const PSO& ant);
	void clear();
};
class candidate_node{
public:
	float total_bw;
	float rest_cpu;
	int node;
	int cpu_rank;
	int bw_rank;
	int price;
	int state;
	candidate_node();
	candidate_node(float _total_bw, float _rest_cpu, int _node, int _cpu_rank, int _bw_rank);
	candidate_node(float _total_bw, float _rest_cpu, int _node, int _cpu_rank, int _bw_rank, int _price);
	candidate_node(float _total_bw, float _rest_cpu, int _node, int _cpu_rank, int _bw_rank, int _price, int _state);
};
class totalW
{
public:
	vector<int> tW;
	vector<int> index;
	vector<int>TW;
};
bool read_sn_info(string sn_dir, SubstrateNetwork &sn);
bool read_req_info(string req_dir, int req_count, vector<Request> &req);
void updateNodes(const SubstrateNetwork &sn, CurrentSn &curSn, int curTime, double &done_cost, int iNode);
void costUpdate(const SubstrateNetwork &sn, CurrentSn &curSn, int curTime, double &done_cost, int &onCount);
float getCost(const SubstrateNetwork &sn, CurrentSn &curSn, int iNode, int base, int duration);
void calc_shortest_path(SubstrateNetwork &sn, struct shortest_path *spath);
bool select_candidate_node(const SubstrateNetwork &sn,
	const CurrentSn &curSn, Request &curReq,
	vector<int> &node_mapping_method, vector<vector<int> > &nodepool);

bool node_mapping_first(const SubstrateNetwork &sn, const CurrentSn &curSn,
	Request &curReq, vector<int> &node_mapping_method, int times);

bool node_mapping(const SubstrateNetwork &sn, const CurrentSn &curSn,
	Request &curReq, vector<int> &node_mapping_method, vector<vector<int> > &nodepool);
bool link_mapping(SubstrateNetwork &sn, CurrentSn &curSn, const Request &curReq,
	vector<int> &node_mapping_method, vector<vector<int> > &link_mapping_method,
struct shortest_path  *spath);


int find_shortest_path(const SubstrateNetwork &sn, const Link &req_link,
	CurrentSn &curSn,
struct shortest_path  *spath,
	int start, int end, int len, vector<int> &sp);
double GetEnergy(const vector<int> &node_mapping_method, const vector<vector<int> > &link_mapping_method,
	const Request &request, const SubstrateNetwork &sn, const CurrentSn &curSn);
bool do_node_mapping(CurrentSn &curSn, Request &curReq, vector<int> &node_mapping_method);
bool do_internode_mapping(CurrentSn &curSn, Request &curReq,
	vector<vector<int> > &link_mapping_method);
bool do_link_mapping(SubstrateNetwork &sn, CurrentSn &curSn, Request &curReq,
	vector<vector<int> > &link_mapping_method);
bool do_mapping(SubstrateNetwork &sn, CurrentSn &curSn, Request &curReq, vector<int> &node_mapping_method,
	vector<vector<int> > &link_mapping_method);
bool release_link(const SubstrateNetwork &sn, CurrentSn &curSn, Request &curReq, int total_nodes);
bool release_cpu(const SubstrateNetwork &sn, CurrentSn &curSn, Request &curReq, int total_nodes);
int comp_cpu(class candidate_node a, class candidate_node b);
int comp_cpu_bestfit(class candidate_node a, class candidate_node b);
int comp_cpu_worstfit(class candidate_node a, class candidate_node b);
int comp_bw(class candidate_node a, class candidate_node b);
int comp_bw_bestfit(class candidate_node a, class candidate_node b);
int comp_bw_worstfit(class candidate_node a, class candidate_node b);
int comp(class candidate_node a, class candidate_node b);
int comp_pr(class candidate_node a, class candidate_node b);
double updateAlpha(SubstrateNetwork &sn, CurrentSn &curSn);
void exchange(Antibody& a, Antibody& b, SubstrateNetwork &sn);
void mutate(const SubstrateNetwork &sn, const CurrentSn &curSn, Request &curReq,
	vector<int> &node_mapping_method, double odd);
template<typename T>
void MemSet(T * array, const T & init_value, int size)
{
	for (int i = 0; i<size; i++)
		array[i] = init_value;
}

int comparevirCPU(class virNode a, class virNode b);
int comparephyCPU(class phyNode a, class phyNode b);
int comparePhyIndex(class phyNode a, class phyNode b);
int  comparephyCPU_reverse(class phyNode a, class phyNode b);
int compare(vector<phyNode> a, vector<phyNode> b);
bool link_mapping_migration(SubstrateNetwork &sn, CurrentSn &curSn, const Request &curReq, vector<int> &node_mapping_method, vector<vector<int> > &link_mapping_method, struct shortest_path  *spath, int vir_index);
bool find_augument_path(int u, vector<virNode> &virNodequeue, vector<int> &visitX, vector<int> &visitY, vector<vector<int> > &edge, vector<int> &matchXY, vector<vector<int> > &matchYX, vector<phyNode> &physicalqueue, SubstrateNetwork &sn, CurrentSn &copy_currentSn, vector<Request> &requests, struct shortest_path *spath);
bool do_link_mapping_migration(const SubstrateNetwork &sn, CurrentSn &copy_curSn, Request &curReq, vector<vector<int> > &temp_link_mapping_method, int &vir_index);
bool checksameVnode(phyNode &phynode, int &reqIndex, int &temp);
bool release_link_migration(const SubstrateNetwork &sn, CurrentSn &copy_curSn, Request &curReq, virNode &virnode);
void do_node_mapping_migration(virNode &virnode, phyNode &phynode, int &u);
void release_node_migration(virNode &virnode, phyNode &phynode, vector<virNode> &virNodequeue);
int find_shortestpath_with_min_offnode(int start, int end, struct shortest_path  *spath, const SubstrateNetwork &sn, vector<int> &shortestpath, CurrentSn &curSn);
bool tpathmapping(SubstrateNetwork &sn, CurrentSn &curSn, int index, vector<int> &shortestpath, float maxbw);
bool releaselinkMIC(const SubstrateNetwork &sn, CurrentSn &curSn, int index);
bool unique_finder(phyNode a, phyNode b);
bool sort_finder(phyNode a, phyNode b);
virNode find_min(vector<virNode> a, vector<virNode> b, vector<vector<int>> c);//find min weight sum
bool check(vector<virNode> a, vector<vector<float>> b);
int get_max(int a, vector<virNode>b, vector<vector<float>> w);//get max weight
int get_min(int a, vector<virNode>b, vector<vector<float>> w);
int cmp(int a, int b);
class Label
{
public:
	vector<float> weight;
	vector<int> label;
};
int get_groupNum(vector<virNode>a);

bool sortfind(virNode a, virNode b);
bool uniquefind(virNode a, virNode b);
int findmax(vector<virNode>a, int b, vector<vector<float>>w_label);
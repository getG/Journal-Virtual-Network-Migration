#include"stdafx.h"
#include "vnEmbed.h"
#include "Stack.h"
#include<numeric>
using namespace std;
extern int iterate_time;

Link::Link(int _from, int _to, float _bw)
	:from(_from), to(_to), bw(_bw)
{
}

SubstrateNetwork::SubstrateNetwork()
{
	this->subInfo = SubInfo();
}

int SubstrateNetwork::GetLinkIndex(int from, int to) const
{
	return this->link_matrix[from][to];
}

candidate_node::candidate_node(float _total_bw, float _rest_cpu, int _node, int _cpu_rank, int _bw_rank)
	:total_bw(_total_bw), rest_cpu(_rest_cpu), node(_node), cpu_rank(_cpu_rank), bw_rank(_bw_rank)
{

}
candidate_node::candidate_node(float _total_bw, float _rest_cpu, int _node, int _cpu_rank, int _bw_rank, int _price)
	: total_bw(_total_bw), rest_cpu(_rest_cpu), node(_node), cpu_rank(_cpu_rank), bw_rank(_bw_rank), price(_price)
{

}
candidate_node::candidate_node(float _total_bw, float _rest_cpu, int _node, int _cpu_rank, int _bw_rank, int _price, int _state)
	: total_bw(_total_bw), rest_cpu(_rest_cpu), node(_node), cpu_rank(_cpu_rank), bw_rank(_bw_rank), price(_price), state(_state)
{

}


bool read_sn_info(string sn_dir, SubstrateNetwork &sn)
{
	FILE * fp_substrate_info;
	string fileName = sn_dir + "/sub.txt";
	fp_substrate_info = fopen(fileName.c_str(), "r");
	if (!fp_substrate_info)
	{
		cout << "open sub.txt failed" << endl;
		return false;
	}
	sn.subInfo.x = 25;
	sn.subInfo.y = 25;

	fscanf(fp_substrate_info, "%d", &(sn.subInfo.nodes));


	fscanf(fp_substrate_info, "%d", &(sn.subInfo.links));


	vector<int> tempLinkMatrixRow(sn.subInfo.nodes, -1);
	sn.link_matrix.resize(sn.subInfo.nodes, tempLinkMatrixRow);
	tempLinkMatrixRow.clear();
	sn.matrix.clear();
	sn.matrix.resize(sn.subInfo.nodes);

	sn.link.resize(sn.subInfo.links);

	int dummy1, dummy2;
	sn.cpu.resize(sn.subInfo.nodes);
	sn.memory.resize(sn.subInfo.nodes);
	sn.totalCPU = 0;
	sn.totalBW = 0;
	for (int j = 0; j<sn.subInfo.nodes; j++)
	{
		for (int i = 0; i<sn.subInfo.nodes; i++)
		{
			sn.matrix[j].push_back(0);
		}
	}

	for (int j = 0; j<sn.subInfo.nodes; j++)
	{
		fscanf(fp_substrate_info, "%d %d %f %f", &dummy1, &dummy2, &(sn.cpu[j]), &(sn.memory[j]));

		sn.totalCPU += sn.cpu[j];

	}

	for (int k = 0; k<sn.subInfo.links; k++)
	{
		int from_node, to_node;
		float rest_bw;
		float dummy3;
		fscanf(fp_substrate_info, "%d %d %f %f", &from_node, &to_node, &rest_bw, &dummy3);
		sn.totalBW += rest_bw;
		sn.link[k].from = from_node;
		sn.link[k].to = to_node;
		sn.link[k].bw = rest_bw;
		int from, to;
		from = sn.link[k].from;
		to = sn.link[k].to;
		sn.link_matrix[from][to] = k;
		sn.link_matrix[to][from] = k;
		sn.matrix[to_node][from_node] = 1;
		sn.matrix[from_node][to_node] = 1;
	}
	fclose(fp_substrate_info);

	//read price info of this sn
	FILE * fp_price_info;
	string price_file;
	price_file.clear();
	price_file = sn_dir + "/price.txt";

	fp_price_info = fopen(price_file.c_str(), "r");
	if (!fp_price_info)
		return false;

	for (int j = 0; j < 50000 / 60 + 1; j++)
	{
		fscanf(fp_price_info, "%f", &(sn.subInfo.price[j]));
		sn.subInfo.price[j] *= 10;
	}
	fclose(fp_price_info);
	//cout<<"read_sn_info success"<<endl;
	return true;
}


bool read_req_info(string req_dir, int req_count, vector<Request> &req)
{

	FILE * fp_req_info = NULL;
	string fileName;
	int no_use_topo;
	for (int i = 0; i<req_count; i++)
	{
		char temp[10];
		fileName.clear();
		sprintf(temp, "%d", i);
		fileName = req_dir + "/req" + temp + ".txt";
		fp_req_info = fopen(fileName.c_str(), "r");
		if (!fp_req_info)
			return false;

		fscanf(fp_req_info, "%d %d %d %d %d %d %d\n",
			&req[i].nodes, &req[i].links, &req[i].split, &req[i].time,
			&req[i].duration, &no_use_topo, &req[i].maxD);
		req[i].duration = req[i].duration;
		//req[i].time=int(req[i].time/100)*100;
		req[i].time = req[i].time;
		req[i].maxD = 10000;
		req[i].node.clear();
		req[i].node.resize(req[i].nodes);
		req[i].req_index = i;

		for (int j = 0; j<req[i].node.size(); j++)
		{
			fscanf(fp_req_info, "%d %d %f %f\n",
				&(req[i].node[j].x), &(req[i].node[j].y), &(req[i].node[j].cpu), &(req[i].node[j].memory));
			//cout<<req[i].node[j].memory<<endl;
			//getchar();
			//req[i].node[j].cpu=req[i].node[j].cpu*2;
		}

		req[i].link.clear();
		req[i].link.resize(req[i].links);
		float no_use_tmp;
		for (int j = 0; j<req[i].link.size(); j++)
		{
			fscanf(fp_req_info, "%d %d %f %f\n",
				&(req[i].link[j].from), &(req[i].link[j].to), &(req[i].link[j].bw), &no_use_tmp);
			//req[i].link[j].bw=req[i].link[j].bw*1.5;
			//new here
			//int from = req[i].link[j].from;
			//int to = req[i].link[j].to;
			//req[i].link[j].bw = (req[i].node[from].cpu + req[i].node[to].cpu) * 50 / 40;
		}

		fclose(fp_req_info);
	}
	return true;
}

// Before function called: adj_metrix & spath were allocated but not initialized
void calc_shortest_path(SubstrateNetwork &sn, struct shortest_path *spath)
{
	// Floyd's All pairs shortest path algorithm (O (n^3) ) 
	// input is adjacency matrix output is matrix of shortest paths
	// adj_matrix is the adjacency matrix
	// n is the order of the square matrix 
	// spath is the all pairs shortest paths matrix 
	// we assume that adj_matrix & spath is allocated by the caller

	int i, j, k;
	int n = sn.subInfo.nodes;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (i == j) {
				(*(spath + i*n + j)).length = 0;//spath[i][i].length=0;
				(*(spath + i*n + j)).next = i;//spath[i][j].next=i;
			}
			else {
				(*(spath + i*n + j)).length = 999999999; //infinity
				(*(spath + i*n + j)).next = -1;//no way back
			}
			(*(spath + i*n + j)).min_bw = 0.0f;
		}
	}

	for (k = 0; k<sn.subInfo.links; k++) {
		if (sn.link[k].from == -1 && sn.link[k].to == -1)
			continue;
		i = sn.link[k].from;
		j = sn.link[k].to;
		(*(spath + i*n + j)).length = 1;//spath[i][j].length=1;
		(*(spath + i*n + j)).next = j;//spath[i][j].next=j;
		(*(spath + j*n + i)).length = 1;
		(*(spath + j*n + i)).next = i;
	}


	// for each route via k from i to j pick 
	// any better routes and replace A[i][j]
	// path with sum of paths i-k and j-k

	for (k = 0; k<n; k++)
	{
		for (i = 0; i<n; i++)
		{
			for (j = 0; j<n; j++)
			{
				if ((*(spath + i*n + k)).length + (*(spath + k*n + j)).length <  (*(spath + i*n + j)).length)
					(*(spath + i*n + j)).length = (*(spath + i*n + k)).length + (*(spath + k*n + j)).length;
				(*(spath + i*n + j)).next = (*(spath + i*n + k)).next;//if the condition is satisfied, then spath[i][j].next will be k for sure
			}
		}
	}

} // Floyd's algorithm

CurrentNode::CurrentNode(float _rest_cpu, float _rest_memory)
{
	this->req_count = 0;
	memset(req, -1, MAX_REQ_PER_NODE*sizeof(int));
	memset(this->vnode, -1, MAX_REQ_PER_NODE*sizeof(int));
	memset(this->cpu, 0, MAX_REQ_PER_NODE*sizeof(float));
	this->rest_cpu = _rest_cpu;
	this->rest_memory = _rest_memory;

}

CurrentLink::CurrentLink(float _rest_bw)
{
	this->link_count = 0;
	this->pre_used = 0;
	memset(this->req, -1, MAX_REQ_LINKS*sizeof(int));
	memset(this->vlink, -1, MAX_REQ_LINKS*sizeof(int));
	memset(this->bw, 0, MAX_REQ_LINKS*sizeof(float));
	this->rest_bw = _rest_bw;
}

CurrentSn::CurrentSn(SubstrateNetwork &sn)
{
	this->time_recent = 0;
	this->currentNode.clear();
	this->currentNode.resize(sn.subInfo.nodes);
	this->currentLink.clear();
	this->currentLink.resize(sn.subInfo.links);

	for (int i = 0; i<currentNode.size(); i++)
	{
		this->currentNode[i].rest_cpu = sn.cpu[i];
		this->currentNode[i].rest_memory = sn.memory[i];

	}

	for (int i = 0; i<currentLink.size(); i++)
	{
		this->currentLink[i].rest_bw = sn.link[i].bw;
	}
}

void PSO::clear()
{
	this->energy = 0;
	this->resource = 0;
	this->x = 0;
	this->y = 0;
	this->node_mapping_method.clear();
	this->comp = 0;
	this->position = 0;
}

void PSO::init()
{
	this->energy = 0;
	this->resource = 0;
	this->comp = 0;
	this->position = 0;
}
void PSO::print()
{
	for (int i = 0; i<this->node_mapping_method.size(); i++)
		cout << this->node_mapping_method[i] << " ";
	cout << endl;
}
void PSO::set(const PSO& ant)
{
	this->energy = ant.energy;
	this->resource = ant.resource;
	this->x = ant.x;
	this->y = ant.y;
	this->comp = ant.comp;
	this->position = ant.position;
	this->node_mapping_method.resize(ant.node_mapping_method.size());
	for (int i = 0; i<this->node_mapping_method.size(); i++)
	{
		this->node_mapping_method[i] = ant.node_mapping_method[i];
	}
}
PSO::PSO()
{
	this->energy = 0;
	this->resource = 0;
	this->x = 0;
	this->y = 0;
	this->node_mapping_method.clear();
	this->comp = 0;
	this->position = 0;

}

PSO::PSO(const PSO& ant)
{
	this->energy = ant.energy;
	this->resource = ant.resource;
	this->x = ant.x;
	this->y = ant.y;
	this->node_mapping_method.clear();
	this->node_mapping_method.resize(ant.node_mapping_method.size());
	this->comp = ant.comp;
	this->position = ant.position;
	for (int i = 0; i<this->node_mapping_method.size(); i++)
	{
		this->node_mapping_method[i] = ant.node_mapping_method[i];
	}
}
float getCost(const SubstrateNetwork &sn, CurrentSn &curSn,
	int iNode, int base, int duration)
{
	float total = sn.cpu[iNode];
	float rest = curSn.currentNode[iNode].rest_cpu;
	float ret;
	if (curSn.currentNode[iNode].state == CurrentNode::NODE_OFF)
		return 0.0;

	float p = sn.subInfo.price[base];
	ret = 1 * (((total - rest) / total) * 150.0f + 165)*duration / 60.0f;//
	return ret;

}

void costUpdate(const SubstrateNetwork &sn, CurrentSn &curSn, int curTime,
	double &done_cost, int &onCount)
{
	if (curTime == curSn.time_recent)
		return;
	onCount = 0;

	for (int i = 0; i<sn.subInfo.nodes; i++)
	{
		if (curSn.currentNode[i].state == CurrentNode::NODE_ON)
		{
			onCount++;
			updateNodes(sn, curSn, curTime, done_cost, i);
		}
	}
	curSn.time_recent = curTime;
}

void updateNodes(const SubstrateNetwork &sn, CurrentSn &curSn, int curTime,
	double &done_cost, int iNode)
{
	float cost = 0;
	int timeRecent = curSn.time_recent;
	cost = getCost(sn, curSn, iNode, timeRecent / 60, curTime - timeRecent);
	done_cost += cost;
	return;


}
void mutate_cpu(vector<candidate_node> &rank)
{
	long *odds = new long[rank.size()];
	long *odds2 = new long[rank.size()];
	for (int i = 0; i<rank.size(); i++)
	{
		odds[i] = odds2[i] = pow(i + 1, 2);
	}
	for (int i = 0; i<rank.size(); i++)
	{
		for (int j = 0; j<rank.size(); j++)
		{
			odds2[j] = odds[j];
		}
		for (int j = 0; j<rank.size() - 1; j++)
		{
			odds2[j + 1] += odds2[j];
		}
		int num = rand() % odds2[rank.size() - 1];
		for (int j = 0; j<rank.size(); j++)
		{
			if (num<odds2[j])
			{
				rank[i].cpu_rank = rank.size() - j;
				odds[j] = 0;
				break;
			}
		}
	}
	delete odds;
	delete odds2;
}

void mutate_bw(vector<candidate_node> &rank)
{
	long *odds = new long[rank.size()];
	long *odds2 = new long[rank.size()];
	for (int i = 0; i<rank.size(); i++)
	{
		odds[i] = odds2[i] = pow(i + 1, 2);
	}
	for (int i = 0; i<rank.size(); i++)
	{
		for (int j = 0; j<rank.size(); j++)
		{
			odds2[j] = odds[j];
		}
		for (int j = 0; j<rank.size() - 1; j++)
		{
			odds2[j + 1] += odds2[j];
		}
		int num = rand() % odds2[rank.size() - 1];
		for (int j = 0; j<rank.size(); j++)
		{
			if (num<odds2[j])
			{
				rank[i].bw_rank = rank.size() - j;
				odds[j] = 0;
				break;
			}
		}
	}
	delete odds;
	delete odds2;
}

bool select_candidate_node(const SubstrateNetwork &sn, const CurrentSn &curSn, Request &curReq,
	vector<int> &node_mapping_method, vector<vector<int> > &nodepool)
{
	vector<int> tmp;
	tmp.clear();

	//rank_cpu by need - rest
	for (int i = 0; i<curReq.nodes; i++)
	{
		float need_cpu = curReq.node[i].cpu;
		float need_memory = curReq.node[i].memory;
		vector<candidate_node> rank_cpu;

		for (int j = 0; j<sn.subInfo.nodes; j++)
		{
			float rest_cpu = curSn.currentNode[j].rest_cpu;
			float rest_memory = curSn.currentNode[j].rest_memory;
			if ((rest_cpu >= need_cpu) && (rest_memory >= need_memory))
			{
				float total_bw = 0.0f;
				for (int k = 0; k<sn.subInfo.links; k++)
				{
					if (sn.link[k].from == j || sn.link[k].to == j)
					{
						total_bw += curSn.currentLink[k].rest_bw;
					}
				}

				rank_cpu.push_back(candidate_node(total_bw, rest_cpu - need_cpu, j, -1, -1));
				nodepool[i].push_back(j);
			}
		}
		if (rank_cpu.size() == 0)
			return false;
		sort(rank_cpu.begin(), rank_cpu.end(), comp_cpu_bestfit);

		mutate_cpu(rank_cpu);

		sort(rank_cpu.begin(), rank_cpu.end(), comp_bw_worstfit);
		mutate_bw(rank_cpu);
		sort(rank_cpu.begin(), rank_cpu.end(), comp);

		//check if it is used before
		int flag = 0;
		for (int m = 0; m<tmp.size(); m++)
		{
			if (rank_cpu[flag].node == tmp[m])
			{
				flag++;
				if (flag >= rank_cpu.size())
					return false;
				m = -1;
			}
		}


		tmp.push_back(rank_cpu[flag].node);
		node_mapping_method[i] = rank_cpu[flag].node;

		rank_cpu.clear();
	}//end of  i=0; i<curReq.nodes; i++
	return true;
}

bool node_mapping(const SubstrateNetwork &sn, const CurrentSn &curSn,
	Request &curReq, vector<int> &node_mapping_method, vector<vector<int> > &nodepool)
{

	if (!(select_candidate_node(sn, curSn, curReq, node_mapping_method, nodepool)))
	{
		return false;
	}
	return true;
}

bool link_mapping(SubstrateNetwork &sn, CurrentSn &curSn, const Request &curReq, vector<int> &node_mapping_method, vector<vector<int> > &link_mapping_method, struct shortest_path  *spath)
{
	int from, to, len;

	for (int i = 0; i<curReq.links; i++)
	{
		int req_from = curReq.link[i].from;
		from = node_mapping_method[req_from];//from=the node at the SN 

		int req_to = curReq.link[i].to;
		to = node_mapping_method[req_to];//to=the node at the SN
		int dis = 1;
		int ret;

		len = spath[from*sn.subInfo.nodes + to].length;
		vector<int> sp(len + 9);
		ret = find_shortest_path(sn, curReq.link[i], curSn, spath, from, to, len, sp);
		while (!ret && dis <= 3)
		{
			sp.clear();
			sp.resize(len + 9);
			ret = find_shortest_path(sn, curReq.link[i], curSn, spath, from, to, len + dis, sp);
			dis++;
		}

		if (ret)//dis >3 ||dis =1
		{
			int sp_len = sp[0];
			sp[0] = from;
			sp[sp_len] = to;
			//check sp
			for (int j = 0; j <= sp_len; j++)
			{
				for (int k = 0; k <= sp_len; k++)
					if (j != k && sp[j] == sp[k])//loop, error
					{
						printf("fatal errorx,exit\n");
						exit(0);
					}
			}
			//pre used
			for (int j = 0; j <= sp_len - 1; j++)
			{
				int from = sp[j];
				int to = sp[j + 1];
				int link_id = sn.GetLinkIndex(from, to);

				curSn.currentLink[link_id].pre_used += curReq.link[i].bw;
			}
			for (int j = 0; j <= sp_len; j++)
				link_mapping_method[i].push_back(sp[j]);
			//here: link_mapping_method[i][j]is the jth step of the ith link

		}//end of if(ret)
		else
		{
			for (int j = 0; j<i; j++)
				link_mapping_method[j].clear();//if false, clear all the paths recorded so far and return a false signal
			return false;
		}

	}
	return true;
}


int find_shortest_path(const SubstrateNetwork &sn, const Link &req_link,
	CurrentSn &curSn, struct shortest_path  *spath,
	int start, int end, int len, vector<int> &sp)
{
	time_t start_t = 0, end_t = 0;
	if (len > 10)
		return 0;
	if (start == end)
	{
		printf("fatal error1(from equals end. I am in find_shortest_path.exit\n"); exit(0);
	}

	if (len == 1)
	{
		sp[0] = 1;//will be collected as the len of the sp
		sp[1] = end;//next?
		int k = sn.GetLinkIndex(start, end);
		if (k == -1)
			return 0;
		if (req_link.bw > (curSn.currentLink[k].rest_bw - curSn.currentLink[k].pre_used))
			return 0;
		return 1;
	}
	int node_off = 0;
	int from;
	int to;

	vector<vector<int> > node(len + 1);
	node[0].push_back(start);
	node[len].push_back(end);

	vector<bool> finded(sn.subInfo.nodes, false);
	finded[start] = finded[end] = true;

	for (int iLen = 1; iLen < len; iLen++)//ilen from 1 to len-1
	{
		for (int iNode = 0; iNode < sn.subInfo.nodes; iNode++)//for all nodes within this sn
		{
			//if iLen equals the length of the SP from start to iNode
			if (spath[start*sn.subInfo.nodes + iNode].length == iLen
				&& !finded[iNode])
			{
				if (iLen >= 2)
				{
					for (int i = 0; i<node[iLen - 1].size(); i++)
					{
						int from = node[iLen - 1][i];
						int iLink = sn.GetLinkIndex(from, iNode);

						if ((iLink != -1) && (curSn.currentLink[iLink].rest_bw - curSn.currentLink[iLink].pre_used >= req_link.bw))
						{
							if (iLen == len - 1)
							{
								int iendlink = -1;
								iendlink = sn.GetLinkIndex(iNode, end);
								if (iendlink != -1 && (curSn.currentLink[iendlink].rest_bw - curSn.currentLink[iendlink].pre_used >= req_link.bw))
								{
									finded[iNode] = true;
									node[iLen].push_back(iNode);
									break;
								}
							}
							else
							{
								finded[iNode] = true;
								node[iLen].push_back(iNode);
								break;
							}

						}
					}
				}
				else	// iLen==1
				{
					int iLink = sn.GetLinkIndex(start, iNode);
					if ((iLink != -1) && (curSn.currentLink[iLink].rest_bw - curSn.currentLink[iLink].pre_used >= req_link.bw))
					{
						finded[iNode] = true;
						node[iLen].push_back(iNode);
					}
				}
			}
		}
	}


	int path_num = 1;
	for (int i = 1; i <= len; i++)
		path_num *= node[i].size();
	if (path_num <= 0)
		return 0;

	Stack<int> s;
	s.push(node[0][0]);

	vector<int> tmp(sn.subInfo.nodes, 0);
	vector<vector<int> > flag(len + 1, tmp);
	tmp.clear();
	flag[0][0] = 1;

	const int MAX_PATHS = 1000000;
	if (path_num >= MAX_PATHS)	// path_num must not exceed MAX_PATH
		path_num = MAX_PATHS;

	vector<int> tmp2(path_num);
	vector<vector<int> > path(len + 1, tmp2);
	tmp2.clear();

	int slen = 1;
	int paths = 0;
	vector<int> cur_num(len + 1, 0);

	while (!s.empty())//
	{

		int from, to;
		int iLink;
		int k;
		if (slen != len)
		{
			for (k = 0; k < node[slen].size(); k++)
			{
				if (flag[slen][k] == 0)
				{
					from = s.top();
					to = node[slen][k];
					iLink = sn.GetLinkIndex(from, to);
					if (iLink != -1 && (curSn.currentLink[iLink].rest_bw - curSn.currentLink[iLink].pre_used >= req_link.bw))//这个K可以用
						break;
					else
						flag[slen][k] = 1;//
				}
			}
			if (k == node[slen].size()) // if not found
			{
				for (int t = 0; t < node[slen].size(); t++)
					flag[slen][t] = 0;
				s.pop();
				slen--;
			}
			else
			{
				s.push(node[slen][k]);
				flag[slen][k] = 1;
				slen++;
			}

		}
		else
		{
			// assure that the last point is adjacent to the end point
			from = s.top();
			to = end;
			iLink = sn.GetLinkIndex(from, to);

			if (iLink != -1 && (curSn.currentLink[iLink].rest_bw - curSn.currentLink[iLink].pre_used >= req_link.bw))
			{
				for (int k = slen - 1, top = s.size() - 1; k >= 0; k--, top--)
				{
					path[k][paths] = s[top];
				}
				paths++;
				if (paths >= MAX_PATHS)
					break;
			}
			s.pop();
			slen--;
		}

	}

	for (int i = 0; i < paths; i++)
		path[0][i] = start;
	for (int i = 0; i < paths; i++)
		path[len][i] = end;
	///////////////////////////////////////////////////////////////////////
	start_t = clock();
	for (int i = 0; i < paths; i++)
	{
		node_off = 0;
		int j, k;
		int Flag;
		for (j = 0; j < len; j++)
		{
			from = path[j][i];
			to = path[j + 1][i];

			if (from != start &&curSn.currentNode[from].rest_cpu < INTER_NODE_CPU)
			{
				path[0][i] = -1;
				break;
			}
			if (to != end && curSn.currentNode[to].rest_cpu < INTER_NODE_CPU)
			{
				path[0][i] = -1;
				break;
			}
			Flag = 0;


			k = sn.GetLinkIndex(from, to);
			if (k != -1)
			{
				if ((curSn.currentLink[k].rest_bw - curSn.currentLink[k].pre_used) < req_link.bw)
				{
					printf("Inside:%d->%d not enough bw\n", from, to);
					path[0][i] = -1;
					Flag = 0;
					break;
				}
				else
				{
					Flag = 1;
					if (curSn.currentNode[to].state == CurrentNode::NODE_OFF)
						node_off++;
				}
			}
			else
			{
				printf("path %d link from %d to %d not exist\n", i, from, to);
				path[0][i] = -1;
				Flag = 0;
				break;
			}
		}
		if (j == len)
		{

			path[0][i] = node_off;
		}
	}
	end_t = clock();

	int min = -1;
	for (int j = 0; j <paths; j++)
	{
		if (path[0][j] != -1)
		{
			min = j;
			break;
		}
	}
	for (int j = 1; j < paths; j++)
	{
		if (path[0][j] != -1 && path[0][j] < path[0][min])
			min = j;
	}
	sp[0] = len;

	if (min == -1)
		return 0;

	for (int i = 1; i < sp[0]; i++)
	{
		sp[i] = path[i][min];
	}

	return 1;
}
double GetNodeEnergy(const SubstrateNetwork &sn, const CurrentSn &curSn,
	int iNode, int start_time, float need_cpu)
{
	float total = sn.cpu[iNode];
	float p = sn.subInfo.price[start_time];
	float ret = 0.0f;
	if (curSn.currentNode[iNode].state == CurrentNode::NODE_ON)
		ret = p * ((need_cpu / total) * 150.0f);
	else
		ret = p * ((need_cpu / total) * 150.0f + 165);

	return ret;
}


double GetEnergy(const vector<int> &node_mapping_method, const vector<vector<int> > &link_mapping_method,
	const Request &request, const SubstrateNetwork &sn, const CurrentSn &curSn)
{
	double cost = 0;
	vector<int> inter_node_count(sn.subInfo.nodes, 0);

	for (int i = 0; i<link_mapping_method.size(); i++)
	{
		for (int j = 1; j<link_mapping_method[i].size() - 1; j++)
		{
			int iNode = link_mapping_method[i][j];
			inter_node_count[iNode]++;
		}
	}

	//get node energy
	for (int i = 0; i<link_mapping_method.size(); i++)//there are i links
	{
		for (int j = 1; j<link_mapping_method[i].size() - 1; j++)//there are j inter nodes
		{
			int timeRecent = request.time;
			int flag = 0;
			int iNode = link_mapping_method[i][j];
			if (inter_node_count[iNode] > 0)//iNode is a internode more than once 
			{
				float need_cpu = INTER_NODE_CPU*inter_node_count[iNode];

				for (int k = 0; k<node_mapping_method.size(); k++)
				{
					if (iNode == node_mapping_method[k])
					{
						need_cpu += request.node[k].cpu;
						flag = -1;
						break;
					}
				}

				cost += GetNodeEnergy(sn, curSn, iNode, timeRecent / 60, need_cpu);
				if (flag == -1)
					inter_node_count[iNode] = -1;
				else
					inter_node_count[iNode] = 0;
				continue;
			}
		}
	}


	//get internode energy
	for (int i = 0; i<node_mapping_method.size(); i++)
	{
		int timeRecent = request.time;
		int iNode = node_mapping_method[i];
		if (inter_node_count[iNode] == -1)//-1 means also a allocated node
			continue;

		float need_cpu = request.node[i].cpu;
		cost += GetNodeEnergy(sn, curSn, iNode, timeRecent / 60, need_cpu);
	}

	return cost;
}





bool do_mapping(SubstrateNetwork &sn, CurrentSn &curSn, Request &curReq, vector<int> &node_mapping_method,
	vector<vector<int> > &link_mapping_method)
{
	if (!do_node_mapping(curSn, curReq, node_mapping_method))
	{
		cout << "do node mapping failed" << endl;
		return false;
	}

	if (!do_internode_mapping(curSn, curReq, link_mapping_method))
	{
		cout << "do internode mapping failed" << endl;
		return false;
	}
	if (!do_link_mapping(sn, curSn, curReq, link_mapping_method))
	{
		cout << "do link mapping failed" << endl;
		return false;
	}
	curReq.node_mapping_method.clear();
	curReq.node_mapping_method.resize(node_mapping_method.size());
	curReq.link_mapping_method.clear();
	curReq.link_mapping_method.resize(link_mapping_method.size());

	curReq.node_mapping_method = node_mapping_method;
	curReq.link_mapping_method = link_mapping_method;
	curReq.temp_node_mapping_method.clear();
	curReq.temp_node_mapping_method.resize(node_mapping_method.size());
	curReq.temp_link_mapping_method.clear();
	curReq.temp_link_mapping_method.resize(link_mapping_method.size());

	curReq.temp_node_mapping_method = node_mapping_method;
	curReq.temp_link_mapping_method = link_mapping_method;


	return true;
}


bool do_node_mapping(CurrentSn &curSn, Request &curReq, vector<int> &node_mapping_method)
{
	int reqNodes = curReq.nodes;

	for (int i = 0; i<reqNodes; i++)
	{
		float need_cpu = curReq.node[i].cpu;
		float need_memory = curReq.node[i].memory;
		int snode = node_mapping_method[i];

		if (curSn.currentNode[snode].rest_cpu < need_cpu)
			return false;

		curSn.currentNode[snode].state = CurrentNode::NODE_ON;

		int req_count = curSn.currentNode[snode].req_count;
		curSn.currentNode[snode].req[req_count] = curReq.req_index;
		curSn.currentNode[snode].vnode[req_count] = i;
		curSn.currentNode[snode].cpu[req_count] = need_cpu;
		curSn.currentNode[snode].memory[req_count] = need_memory;
		curSn.currentNode[snode].rest_cpu -= need_cpu;
		curSn.currentNode[snode].rest_memory -= need_memory;
		curSn.currentNode[snode].req_count++;
	}
	return true;
}
bool do_internode_mapping(CurrentSn &curSn, Request &curReq,
	vector<vector<int> > &link_mapping_method)
{
	int links = curReq.links;
	for (int i = 0; i<links; i++)
	{
		int len = link_mapping_method[i].size();
		if (len >2)
		{
			for (int j = 1; j<len - 1; j++)
			{
				int snode = link_mapping_method[i][j];

				int req_count = curSn.currentNode[snode].req_count;
				curSn.currentNode[snode].req[req_count] = curReq.req_index;
				curSn.currentNode[snode].vnode[req_count] = -1;
				curSn.currentNode[snode].cpu[req_count] = INTER_NODE_CPU;
				curSn.currentNode[snode].rest_cpu -= INTER_NODE_CPU;
				curSn.currentNode[snode].req_count++;


				if (curSn.currentNode[snode].state == CurrentNode::NODE_OFF)
					curSn.currentNode[snode].state = CurrentNode::NODE_ON;
			}
		}
	}

	return true;
}
bool do_link_mapping(SubstrateNetwork &sn, CurrentSn &curSn, Request &curReq,
	vector<vector<int> > &link_mapping_method)
{
	int links = curReq.links;

	for (int i = 0; i<links; i++)
	{
		float need_bw = curReq.link[i].bw;
		int len = link_mapping_method[i].size();
		for (int j = 0; j<len - 1; j++)
		{
			int from = link_mapping_method[i][j];
			int to = link_mapping_method[i][j + 1];

			int link_index = sn.GetLinkIndex(from, to);
			if (link_index == -1)
			{
				cout << "no such link" << endl;
				return false;
			}
			if (curSn.currentLink[link_index].rest_bw < need_bw)
			{
				cout << i << endl;
				cout << "There is not enough bw " << endl;
				return false;
			}
			int req_count = curSn.currentLink[link_index].link_count;
			curSn.currentLink[link_index].req[req_count] = curReq.req_index;
			curSn.currentLink[link_index].vlink[req_count] = i;
			curSn.currentLink[link_index].bw[req_count] = need_bw;
			curSn.currentLink[link_index].rest_bw -= need_bw;
			curSn.currentLink[link_index].link_count++;
		}
	}
	return true;
}

int comp_cpu_bestfit(class candidate_node a, class candidate_node b)
{
	if (a.state == CurrentNode::NODE_OFF && b.state == CurrentNode::NODE_ON)
		return 0;
	if (a.state == CurrentNode::NODE_ON && b.state == CurrentNode::NODE_OFF)
		return 1;


	return a.rest_cpu < b.rest_cpu;
}

int comp_bw_worstfit(class candidate_node a, class candidate_node b)
{
	return a.total_bw > b.total_bw;
}

int comp(class candidate_node a, class candidate_node b)
{
	return a.bw_rank + a.cpu_rank < b.bw_rank + b.cpu_rank;
}



bool release_link(const SubstrateNetwork &sn, CurrentSn &curSn, Request &curReq, int total_nodes)
{

	int curIndex = curReq.req_index;

	for (int i = 0; i<sn.subInfo.links; i++)
	{
		for (int j = 0; j<MAX_REQ_LINKS; j++)
		{
			if (curSn.currentLink[i].req[j] == curIndex)
			{
				int link_count = curSn.currentLink[i].link_count;
				if (link_count == 1)
				{
					curSn.currentLink[i].link_count = 0;
					curSn.currentLink[i].rest_bw += curSn.currentLink[i].bw[j];
					curSn.currentLink[i].bw[j] = 0;
					curSn.currentLink[i].vlink[j] = -1;
					curSn.currentLink[i].req[j] = -1;
				}
				else
				{
					curSn.currentLink[i].link_count--;
					curSn.currentLink[i].rest_bw += curSn.currentLink[i].bw[j];

					curSn.currentLink[i].req[j] = curSn.currentLink[i].req[link_count - 1];
					curSn.currentLink[i].vlink[j] = curSn.currentLink[i].vlink[link_count - 1];
					curSn.currentLink[i].bw[j] = curSn.currentLink[i].bw[link_count - 1];

					curSn.currentLink[i].req[link_count - 1] = -1;
					curSn.currentLink[i].vlink[link_count - 1] = -1;
					curSn.currentLink[i].bw[link_count - 1] = 0;
					j--;

				}
			}
		}
	}//end of int i=0; i<sn.subInfo.links

	//just check
	for (int i = 0; i<sn.subInfo.links; i++)
	{
		float rest_bw = curSn.currentLink[i].rest_bw;
		for (int j = 0; j<MAX_REQ_LINKS; j++)
		{
			if (curSn.currentLink[i].req[j] != -1)
			{
				rest_bw += curSn.currentLink[i].bw[j];
			}
		}

		float tmp = sn.link[i].bw;
		if (fabs(rest_bw - tmp) >= 0.01)
		{
			cout << "there is something wrong with releasing link" << endl;
			return false;
		}
	}
	return true;
}
bool release_cpu(const SubstrateNetwork &sn, CurrentSn &curSn, Request &curReq, int total_nodes)
{
	int curIndex = curReq.req_index;

	for (int i = 0; i<total_nodes; i++)
	{
		for (int j = 0; j<MAX_REQ_PER_NODE; j++)
		{
			if (curSn.currentNode[i].req[j] == curIndex)
			{
				int req_count = curSn.currentNode[i].req_count;
				if (req_count == 1)
				{
					curSn.currentNode[i].req_count = 0;
					curSn.currentNode[i].rest_cpu += curSn.currentNode[i].cpu[j];
					curSn.currentNode[i].rest_memory += curSn.currentNode[i].memory[j];
					curSn.currentNode[i].cpu[j] = 0;
					curSn.currentNode[i].memory[j] = 0;
					curSn.currentNode[i].vnode[j] = -1;
					curSn.currentNode[i].req[j] = -1;
					curSn.currentNode[i].state = CurrentNode::NODE_OFF;
				}
				else
				{
					curSn.currentNode[i].req_count--;
					curSn.currentNode[i].rest_cpu += curSn.currentNode[i].cpu[j];
					curSn.currentNode[i].rest_memory += curSn.currentNode[i].memory[j];

					curSn.currentNode[i].req[j] = curSn.currentNode[i].req[req_count - 1];
					curSn.currentNode[i].vnode[j] = curSn.currentNode[i].vnode[req_count - 1];
					curSn.currentNode[i].cpu[j] = curSn.currentNode[i].cpu[req_count - 1];
					curSn.currentNode[i].memory[j] = curSn.currentNode[i].memory[req_count - 1];
					curSn.currentNode[i].req[req_count - 1] = -1;
					curSn.currentNode[i].vnode[req_count - 1] = -1;
					curSn.currentNode[i].cpu[req_count - 1] = 0;
					j--;
				}
			}
		}
	}

	//just check

	for (int i = 0; i<total_nodes; i++)
	{
		float rest_cpu = curSn.currentNode[i].rest_cpu;
		for (int j = 0; j<MAX_REQ_PER_NODE; j++)
		{
			if (curSn.currentNode[i].req[j] != -1)
			{
				rest_cpu += curSn.currentNode[i].cpu[j];
			}
		}

		float tmp = sn.cpu[i];
		if (fabs(rest_cpu - tmp) >= 0.01)
		{
			cout << "There is something wrong with release_cpu (just check) " << endl;
			return false;
		}
	}

	return true;
}

double updateAlpha(SubstrateNetwork &sn, CurrentSn &curSn)
{
	double restCPU = 0;
	double restBW = 0;

	for (int i = 0; i<sn.subInfo.nodes; i++)
	{
		restCPU += curSn.currentNode[i].rest_cpu;
	}
	for (int i = 0; i<sn.subInfo.links; i++)
	{
		restBW += curSn.currentLink[i].rest_bw;
	}
	return (sn.totalCPU + sn.totalBW - restCPU - restBW) / (sn.totalCPU + sn.totalBW);
}
void random_array(int a[], int size)
{
	srand(time(0));
	int r;
	int i, j;
	int flag;
	for (i = 0; i < size; i++)
	{
		flag = 0;
		while (flag == 0)
		{
			r = rand() % size;
			for (j = 0; j < i; j++)
				if (r == a[j])
					break;
			if (j == i)
			{
				flag = 1;
				a[i] = r;
			}

		}
	}
}



void mutate(const SubstrateNetwork &sn, const CurrentSn &curSn, Request &curReq,
	vector<int> &node_mapping_method, double odd)
{
	for (int i = 0; i<node_mapping_method.size(); i++)
	{
		double odd2 = double((rand() % 100 + 1)) / 100.0;
		if (odd2<odd)
		{
			float need_cpu = curReq.node[i].cpu;

			vector<candidate_node> rank_cpu;
			for (int j = 0; j<sn.subInfo.nodes; j++)
			{
				float rest_cpu = curSn.currentNode[j].rest_cpu;
				if (rest_cpu >= need_cpu)
				{
					float total_bw = 0.0f;
					for (int k = 0; k<sn.subInfo.links; k++)
					{
						if (sn.link[k].from == j || sn.link[k].to == j)
						{
							total_bw += curSn.currentLink[k].rest_bw;
						}
					}
					rank_cpu.push_back(candidate_node(total_bw, rest_cpu - need_cpu, j, -1, -1));
				}
			}
			if (rank_cpu.size() == 0)
				return;
			sort(rank_cpu.begin(), rank_cpu.end(), comp_cpu_bestfit);
			mutate_cpu(rank_cpu);

			sort(rank_cpu.begin(), rank_cpu.end(), comp_bw_worstfit);
			mutate_bw(rank_cpu);

			sort(rank_cpu.begin(), rank_cpu.end(), comp);

			//check if it is used before
			int flag = 0;
			for (int m = 0; m<node_mapping_method.size(); m++)
			{
				if (m == i)
					continue;
				if (rank_cpu[flag].node == node_mapping_method[m])
				{
					flag++;
					if (flag >= rank_cpu.size())
						return;
					m = -1;
				}
			}
			node_mapping_method[i] = rank_cpu[flag].node;

		}
	}
}


double Request::GetRevenue()
{
	float cpu_rev = 0;
	float bw_rev = 0;

	for (int i = 0; i<this->nodes; i++)
	{
		cpu_rev += node[i].cpu;
	}

	for (int i = 0; i<this->links; i++)
	{
		bw_rev += link[i].bw;
	}

	return (cpu_rev + bw_rev);
}



///////////////////////////
int  comparevirCPU(class virNode a, class virNode b)
{
	if (a.totalcpu != b.totalcpu)
		return a.totalcpu > b.totalcpu;
	if (a.totalcpu == a.totalcpu)
		return 	a.cpu>b.cpu;

}
int  comparephyCPU(class phyNode a, class phyNode b)
{
	if (a.hot != b.hot)
		return a.hot > b.hot;
	if (a.hot == b.hot)
		return a.rest_cpu < b.rest_cpu;
}

int  comparephyCPU_reverse(class phyNode a, class phyNode b)
{
	if (a.hot != b.hot)
		return a.hot < b.hot;
	if (a.hot == b.hot)
		return a.rest_cpu > b.rest_cpu;
}
int compare(vector<phyNode> a, vector<phyNode> b)
{
	//int w[m][n] = -1;

	int	w = -1;
	int flag = 0;
	for (int i = 0; i<b.size(); i++)
	{
		for (int j = 0; j<a.size(); j++)
		{
			//int w[m][n] = -1;
			if (a[j].physicalNodeIndex == b[i].physicalNodeIndex)
			{
				w = pow(i, 1) + pow(j, 1);
				flag = 1;
				break;
			}
		}
		if (flag == 1)
		{
			break;
		}
	}
	return w;
}
bool link_mapping_migration(SubstrateNetwork &sn, CurrentSn &curSn, const Request &curReq, vector<int> &node_mapping_method, vector<vector<int> > &link_mapping_method, struct shortest_path  *spath, int vir_index)
{
	int from, to, len;

	for (int i = 0; i<curReq.links; i++)
	{
		int req_from = curReq.link[i].from;
		from = node_mapping_method[req_from];//from,the node at the SN 

		int req_to = curReq.link[i].to;
		to = node_mapping_method[req_to];//to, the node at the SN
		if ((req_from == vir_index) || (req_to == vir_index))
		{
			int dis = 1;
			int ret;

			len = spath[from*sn.subInfo.nodes + to].length;
			vector<int> sp(len + 9);
			ret = find_shortest_path(sn, curReq.link[i], curSn, spath, from, to, len, sp);
			while (!ret && dis <= 3)
			{
				sp.clear();
				sp.resize(len + 9);
				ret = find_shortest_path(sn, curReq.link[i], curSn, spath, from, to, len + dis, sp);
				dis++;
			}

			if (ret)//dis >3 ||dis =1
			{
				int sp_len = sp[0];
				sp[0] = from;
				sp[sp_len] = to;
				//check sp
				for (int j = 0; j <= sp_len; j++)
				{
					for (int k = 0; k <= sp_len; k++)
						if (j != k && sp[j] == sp[k])//loop, error
						{
							printf("fatal errorx,exit,\n");
							exit(0);

						}
				}
				//pre used
				for (int j = 0; j <= sp_len - 1; j++)
				{
					int from = sp[j];
					int to = sp[j + 1];
					int link_id = sn.GetLinkIndex(from, to);

					curSn.currentLink[link_id].pre_used += curReq.link[i].bw;
				}
				for (int j = 0; j <= sp_len; j++)
					link_mapping_method[i].push_back(sp[j]);
				//here: link_mapping_method[i][j]is the jth step of the ith link

			}//end of if(ret)
			else
			{
				for (int j = 0; j<i; j++)
					link_mapping_method[j].clear();//if false, clear all the paths recorded so far and return a false signal
				return false;
			}

		}//if(from==objectiveNode||to==objectiveNode)

	}
	return true;
}

bool find_augument_path(int u, vector<virNode> &virNodequeue, vector<int> &visitX, vector<int> &visitY, vector<vector<int> > &edge, vector<int> &matchXY, vector<vector<int> > &matchYX,
	vector<phyNode> &physicalqueue, SubstrateNetwork &sn, CurrentSn &copy_currentSn, vector<Request> &requests, struct shortest_path *spath)
{
	int temp = -5;
	int reqindex = virNodequeue[u].req_index;
	int vir_index = virNodequeue[u].vir_index;
	int req_from;
	int req_to;
	int loop;
	virNode viru = virNodequeue[u];
	CurrentSn temp_copy_currentSn(sn);
	vector<phyNode> temp_physicalqueue;
	vector<vector<int> > temp_matchYX;
	vector<int> temp_node_mapping_method;
	vector<vector<int> > temp_link_mapping_method;
	temp_node_mapping_method.clear();
	temp_link_mapping_method.clear();
	temp_physicalqueue.clear();
	temp_physicalqueue.resize(physicalqueue.size());
	temp_matchYX.clear();
	temp_matchYX.resize(matchYX.size());
	visitX[u] = 1;
	//cout<<"haha0"<<endl;

	for (int v = 0; v<physicalqueue.size(); v++)
	{

		if ((edge[u][v] == 1) && (matchXY[u] != v))
		{
			//cout<<u<<" "<<v<<endl;
			//getchar();
			reqindex = virNodequeue[u].req_index;

			visitY[v] = 1;
			temp_node_mapping_method.resize(requests[reqindex].temp_node_mapping_method.size());
			temp_node_mapping_method = requests[reqindex].temp_node_mapping_method;//////////////////////////
			temp_node_mapping_method[vir_index] = physicalqueue[v].physicalNodeIndex;//
			temp_link_mapping_method.resize(requests[reqindex].temp_link_mapping_method.size());
			temp_link_mapping_method = requests[reqindex].temp_link_mapping_method;//link提取////////////////////////
			//cout<<"haha2"<<endl;

			loop = 0;
			for (int i = 0; i<temp_node_mapping_method.size(); i++)
			{
				if ((i != vir_index) && (temp_node_mapping_method[vir_index] == temp_node_mapping_method[i]))
				{
					loop = 1;
					break;
				}
			}
			//cout<<"loop"<<loop<<endl;

			for (int i = 0; i<requests[reqindex].links; i++)
			{
				int req_from = requests[reqindex].link[i].from;
				int req_to = requests[reqindex].link[i].to;
				if ((vir_index == req_from) || (vir_index == req_to))//vir_index req_from req_to 
				{
					temp_link_mapping_method[i].clear();
				}

			}

			//cout<<"haha3"<<endl;
			//cout<<"loop"<<loop<<endl;
			if (loop == 0)
			{

				if ((physicalqueue[v].rest_cpu>virNodequeue[u].cpu) && (physicalqueue[v].rest_memory >= virNodequeue[u].memory) && (!checksameVnode(physicalqueue[v], reqindex, temp)))//if(true)加上memory
				{
					//cout<<"success"<<endl;
					//getchar();
					for (int i = 0; i<sn.subInfo.links; i++)
					{
						copy_currentSn.currentLink[i].pre_used = 0;
					}

					do_link_mapping_migration(sn, copy_currentSn, requests[reqindex], temp_link_mapping_method, vir_index);
					do_node_mapping_migration(virNodequeue[u], physicalqueue[v], u);
					matchYX[v].push_back(u);
					matchXY[u] = v;
					requests[reqindex].temp_node_mapping_method = temp_node_mapping_method;
					requests[reqindex].temp_link_mapping_method = temp_link_mapping_method;


					return true;
				}

				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

				else{

					//getchar();
					for (int i = 0; i<sn.subInfo.links; i++)
					{
						copy_currentSn.currentLink[i].pre_used = 0;
					}
					checksameVnode(physicalqueue[v], reqindex, temp);
					//cout<<"temp"<<temp<<endl;
					//getchar();
					if (temp == -1)
					{

						//getchar();

						for (int i = 0; i<matchYX[v].size(); i++)
						{
							int release_in_phy = i + physicalqueue[v].orireq_count;
							int search1;

							search1 = matchYX[v][i];







							temp_copy_currentSn = copy_currentSn;
							temp_physicalqueue = physicalqueue;
							temp_matchYX = matchYX;

							release_link_migration(sn, copy_currentSn, requests[virNodequeue[matchYX[v][i]].req_index], virNodequeue[matchYX[v][i]]);//virNodequeue[match[v][i]].req_index
							release_node_migration(virNodequeue[matchYX[v][i]], physicalqueue[v], virNodequeue);
							int matchsize = matchYX[v].size() - 1;
							matchYX[v][i] = matchYX[v][matchsize];
							matchYX[v].resize(matchYX[v].size() - 1);








							if ((physicalqueue[v].rest_cpu >= virNodequeue[u].cpu) && (physicalqueue[v].rest_memory >= virNodequeue[u].memory))
							{
								for (int i = 0; i<sn.subInfo.links; i++)
								{
									copy_currentSn.currentLink[i].pre_used = 0;
								}

								do_link_mapping_migration(sn, copy_currentSn, requests[reqindex], temp_link_mapping_method, vir_index);
								do_node_mapping_migration(virNodequeue[u], physicalqueue[v], u);
								matchYX[v].push_back(u);
								requests[reqindex].temp_node_mapping_method = temp_node_mapping_method;
								requests[reqindex].temp_link_mapping_method = temp_link_mapping_method;
								//cout<<"match上了2"<<endl;






								/////////////////////////
								if (find_augument_path(search1, virNodequeue, visitX, visitY, edge, matchXY, matchYX, physicalqueue, sn, copy_currentSn, requests, spath))
								{
									for (int i = 0; i<sn.subInfo.links; i++)
									{
										copy_currentSn.currentLink[i].pre_used = 0;
									}

									do_link_mapping_migration(sn, copy_currentSn, requests[reqindex], temp_link_mapping_method, vir_index);
									do_node_mapping_migration(virNodequeue[u], physicalqueue[v], u);
									matchXY[u] = v;
									requests[reqindex].temp_node_mapping_method = temp_node_mapping_method;
									requests[reqindex].temp_link_mapping_method = temp_link_mapping_method;

									return true;
								}
								else
								{
									for (int i = 0; i<sn.subInfo.links; i++)
									{
										copy_currentSn.currentLink[i].pre_used = 0;
									}

									copy_currentSn = temp_copy_currentSn;
									physicalqueue = temp_physicalqueue;

									matchYX.resize(temp_matchYX.size());
									matchYX = temp_matchYX;

								}


								//////////////////////////


							}
							for (int i = 0; i<sn.subInfo.links; i++)
							{
								copy_currentSn.currentLink[i].pre_used = 0;
							}
						}

						//////////////////////////////////////////////////////////                          
					}//if(temp==-1)
					if (temp == -2)
					{

						//getchar();
						return false;

					}
					if (temp >= 0)
					{


						int release_in_match = temp - physicalqueue[v].orireq_count;
						int search;
						temp_copy_currentSn = copy_currentSn;
						temp_physicalqueue = physicalqueue;
						temp_matchYX = matchYX;




						search = matchYX[v][release_in_match];
						release_link_migration(sn, copy_currentSn, requests[virNodequeue[matchYX[v][release_in_match]].req_index], virNodequeue[matchYX[v][release_in_match]]);
						release_node_migration(virNodequeue[matchYX[v][release_in_match]], physicalqueue[v], virNodequeue);
						int matchsize = matchYX[v].size() - 1;
						matchYX[v][release_in_match] = matchYX[v][matchsize];
						matchYX[v].resize(matchYX[v].size() - 1);

						if ((physicalqueue[v].rest_cpu >= virNodequeue[u].cpu) && (physicalqueue[v].rest_memory >= virNodequeue[u].memory))
						{
							for (int i = 0; i<sn.subInfo.links; i++)
							{
								copy_currentSn.currentLink[i].pre_used = 0;
							}


							do_link_mapping_migration(sn, copy_currentSn, requests[reqindex], temp_link_mapping_method, vir_index);
							do_node_mapping_migration(virNodequeue[u], physicalqueue[v], u);
							matchYX[v].push_back(u);
							requests[reqindex].temp_node_mapping_method = temp_node_mapping_method;
							requests[reqindex].temp_link_mapping_method = temp_link_mapping_method;


							if (find_augument_path(search, virNodequeue, visitX, visitY, edge, matchXY, matchYX, physicalqueue, sn, copy_currentSn, requests, spath))
							{

								for (int i = 0; i<sn.subInfo.links; i++)
								{
									copy_currentSn.currentLink[i].pre_used = 0;
								}

								do_link_mapping_migration(sn, copy_currentSn, requests[reqindex], temp_link_mapping_method, vir_index);
								do_node_mapping_migration(virNodequeue[u], physicalqueue[v], u);

								matchXY[u] = v;
								//cout<<"haha"<<endl;
								requests[reqindex].temp_node_mapping_method = temp_node_mapping_method;
								requests[reqindex].temp_link_mapping_method = temp_link_mapping_method;

								return true;
							}
							else
							{
								for (int i = 0; i<sn.subInfo.links; i++)
								{
									copy_currentSn.currentLink[i].pre_used = 0;
								}

								copy_currentSn = temp_copy_currentSn;
								physicalqueue = temp_physicalqueue;

								matchYX.resize(temp_matchYX.size());
								matchYX = temp_matchYX;

							}


						}






						for (int i = 0; i<sn.subInfo.links; i++)
						{
							copy_currentSn.currentLink[i].pre_used = 0;
						}


					}






				}
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////








			}//if(loop==0)
			else
			{
				;
			}




		} 

	}
	return false;
}


bool checksameVnode(phyNode &phynode, int &reqindex, int &temp)
{

	for (int i = 0; i<phynode.req_count; i++)
	{
		if (reqindex == phynode.req[i])
		{
			if (i>(phynode.orireq_count - 1))
			{
				temp = i;

				return true;
			}
			else
			{
				temp = -2;
				return true;
			}
		}
	}
	temp = -1;
	return false;
}


bool do_link_mapping_migration(const SubstrateNetwork &sn, CurrentSn &copy_curSn, Request &curReq, vector<vector<int> > &temp_link_mapping_method, int &vir_index)//这个函数只映射上刚刚补上的缺失链路
{
	int links = curReq.links;
	int req_from;
	int req_to;
	float need_bw;
	int len;
	for (int i = 0; i<links; i++)
	{
		need_bw = curReq.link[i].bw;
		len = temp_link_mapping_method[i].size();
		req_from = curReq.link[i].from;
		req_to = curReq.link[i].to;
		if ((vir_index == req_from) || (vir_index == req_to))
		{
			for (int j = 0; j<len - 1; j++)
			{
				int from = temp_link_mapping_method[i][j];
				int to = temp_link_mapping_method[i][j + 1];

				int link_index = sn.GetLinkIndex(from, to);
				if (link_index == -1)
				{
					//cout<<"no such link"<<endl;
					//getchar();
					return false;
				}
				if (copy_curSn.currentLink[link_index].rest_bw < need_bw)
				{
					//cout<<i<<endl;
					//cout<<"There is not enough bw "<<endl;
					//getchar();
					return false;
				}
				int req_count = copy_curSn.currentLink[link_index].link_count;
				copy_curSn.currentLink[link_index].req[req_count] = curReq.req_index;
				copy_curSn.currentLink[link_index].vlink[req_count] = i;
				copy_curSn.currentLink[link_index].bw[req_count] = need_bw;
				copy_curSn.currentLink[link_index].rest_bw -= need_bw;
				copy_curSn.currentLink[link_index].link_count++;
			}
		}



	}

	return true;
}

bool release_link_migration(const SubstrateNetwork &sn, CurrentSn &copy_curSn, Request &curReq, virNode &virnode)											 	{

	int curIndex = curReq.req_index;
	int req_from;
	int req_to;
	int temp = 0;
	int from;
	int to;
	int len;
	int link_index;
	for (int i = 0; i<curReq.links; i++)
	{
		req_from = curReq.link[i].from;
		req_to = curReq.link[i].to;
		if ((virnode.vir_index == req_from) || (virnode.vir_index == req_to))
		{
			len = curReq.temp_link_mapping_method[i].size();
			for (int j = 0; j<len - 1; j++)
			{
				from = curReq.temp_link_mapping_method[i][j];
				to = curReq.temp_link_mapping_method[i][j + 1];
				link_index = sn.GetLinkIndex(from, to);

				for (int k = 0; k<MAX_REQ_LINKS; k++)
				{

					if (copy_curSn.currentLink[link_index].req[k] == curIndex)
					{


						int link_count = copy_curSn.currentLink[link_index].link_count;
						if (link_count == 1)
						{
							copy_curSn.currentLink[link_index].link_count = 0;
							copy_curSn.currentLink[link_index].rest_bw += copy_curSn.currentLink[link_index].bw[k];
							copy_curSn.currentLink[link_index].bw[k] = 0;
							copy_curSn.currentLink[link_index].vlink[k] = -1;
							copy_curSn.currentLink[link_index].req[k] = -1;
						}
						else
						{
							copy_curSn.currentLink[link_index].link_count--;
							copy_curSn.currentLink[link_index].rest_bw += copy_curSn.currentLink[link_index].bw[k];

							copy_curSn.currentLink[link_index].req[k] = copy_curSn.currentLink[link_index].req[link_count - 1];
							copy_curSn.currentLink[link_index].vlink[k] = copy_curSn.currentLink[link_index].vlink[link_count - 1];
							copy_curSn.currentLink[link_index].bw[k] = copy_curSn.currentLink[link_index].bw[link_count - 1];

							copy_curSn.currentLink[link_index].req[link_count - 1] = -1;
							copy_curSn.currentLink[link_index].vlink[link_count - 1] = -1;
							copy_curSn.currentLink[link_index].bw[link_count - 1] = 0;
							k--;
						}

					}
				}
			}


		}


	}








	//just check
	for (int i = 0; i<sn.subInfo.links; i++)
	{
		float rest_bw = copy_curSn.currentLink[i].rest_bw;
		for (int j = 0; j<MAX_REQ_LINKS; j++)
		{
			if (copy_curSn.currentLink[i].req[j] != -1)
			{
				rest_bw += copy_curSn.currentLink[i].bw[j];
			}
		}

		float tmp = sn.link[i].bw;
		if (fabs(rest_bw - tmp) >= 0.01)
		{
			cout << "there is something wrong with releasing link" << endl;
			getchar();
			return false;
		}
	}
	return true;
}

void do_node_mapping_migration(virNode &virnode, phyNode &phynode, int &u)
{
	int req_count = phynode.req_count;


	virnode.jindex = req_count;
	phynode.rest_cpu -= virnode.cpu;
	phynode.rest_memory -= virnode.memory;
	phynode.req[req_count] = virnode.req_index;
	phynode.vnode[req_count] = virnode.vir_index;
	phynode.cpu[req_count] = virnode.cpu;
	phynode.memory[req_count] = virnode.memory;
	phynode.virqueueindex[req_count] = u;
	phynode.req_count++;
}
void release_node_migration(virNode &virnode, phyNode &phynode, vector<virNode> &virNodequeue)
{

	int req_count = phynode.req_count;
	int virnode_index_inqueue = phynode.virqueueindex[req_count - 1];

	phynode.rest_cpu += virnode.cpu;
	phynode.rest_memory += virnode.memory;
	phynode.req[virnode.jindex] = phynode.req[req_count - 1];
	phynode.vnode[virnode.jindex] = phynode.vnode[req_count - 1];
	phynode.cpu[virnode.jindex] = phynode.cpu[req_count - 1];
	phynode.memory[virnode.jindex] = phynode.memory[req_count - 1];
	phynode.virqueueindex[virnode.jindex] = phynode.virqueueindex[req_count - 1];
	virNodequeue[virnode_index_inqueue].jindex = virnode.jindex;
	phynode.req_count--;
	phynode.req[req_count - 1] = -1;
	phynode.vnode[req_count - 1] = -1;
	phynode.virqueueindex[req_count - 1] = -1;
	phynode.cpu[req_count - 1] = 0;
	phynode.memory[req_count - 1] = 0;



}

///////////////////////////////////////////////////////						
int find_shortestpath_with_min_offnode(int start, int end, struct shortest_path  *spath, const SubstrateNetwork &sn, vector<int> &shortestpath, CurrentSn &curSn)
{
	int len = spath[start*sn.subInfo.nodes + end].length;

	time_t start_t = 0, end_t = 0;
	shortestpath.resize(len + 1);
	cout << "len" << len << endl;
	//getchar();

	if (len > 10)
		return 0;
	if (start == end)
	{
		printf("fatal error1(from equals end. I am in find_shortestpath_with_min_offnode.exit\n"); exit(0);
	}

	if (len == 1)
	{

		int k = sn.GetLinkIndex(start, end);
		if (k == -1)
			return 0;
		else
		{
			shortestpath[0] = start;
			shortestpath[1] = end;

			return 1;
		}

	}
	int node_off = 0;
	int from;
	int to;

	vector<vector<int> > node(len + 1);
	node[0].push_back(start);
	node[len].push_back(end);

	vector<bool> finded(sn.subInfo.nodes, false);
	finded[start] = finded[end] = true;

	for (int iLen = 1; iLen < len; iLen++)//ilen from 1 to len-1
	{
		for (int iNode = 0; iNode < sn.subInfo.nodes; iNode++)//for all nodes within this sn
		{
			//if iLen equals the length of the SP from start to iNode
			if ((spath[start*sn.subInfo.nodes + iNode].length == iLen) && (!finded[iNode]))
			{
				if (iLen >= 2)
				{
					for (int i = 0; i<node[iLen - 1].size(); i++)
					{
						int from = node[iLen - 1][i];
						int iLink = sn.GetLinkIndex(from, iNode);

						if (iLink != -1)
						{
							if (iLen == len - 1)
							{
								int iendlink = -1;
								iendlink = sn.GetLinkIndex(iNode, end);
								if (iendlink != -1)
								{
									finded[iNode] = true;
									node[iLen].push_back(iNode);//node[i][j] means the jth node of the ith sp
									break;
								}
							}
							else
							{
								finded[iNode] = true;
								node[iLen].push_back(iNode);
								break;
							}

						}
					}
				}
				else	// iLen==1
				{
					int iLink = sn.GetLinkIndex(start, iNode);
					if (iLink != -1)
					{
						finded[iNode] = true;
						node[iLen].push_back(iNode);
					}

				}
			}
		}

	}



	int path_num = 1;
	for (int i = 1; i <= len; i++)
		path_num *= node[i].size();
	if (path_num <= 0)
		return 0;

	Stack<int> s;
	s.push(node[0][0]);

	vector<int> tmp(sn.subInfo.nodes, 0);
	vector<vector<int> > flag(len + 1, tmp);
	tmp.clear();
	flag[0][0] = 1;

	const int MAX_PATHS = 1000000;
	if (path_num >= MAX_PATHS)	// path_num must not exceed MAX_PATH
		path_num = MAX_PATHS;

	vector<int> tmp2(path_num);
	vector<vector<int> > path(len + 1, tmp2);
	tmp2.clear();

	int slen = 1;
	int paths = 0;
	vector<int> cur_num(len + 1, 0);

	while (!s.empty())
	{

		int from, to;
		int iLink;
		int k;
		if (slen != len)
		{
			for (k = 0; k < node[slen].size(); k++)
			{
				if (flag[slen][k] == 0)
				{
					from = s.top();
					to = node[slen][k];
					iLink = sn.GetLinkIndex(from, to);
					if (iLink != -1)
						break;
					else
						flag[slen][k] = 1;
				}
			}
			if (k == node[slen].size()) // if not found
			{
				for (int t = 0; t < node[slen].size(); t++)
					flag[slen][t] = 0;
				s.pop();
				slen--;
			}
			else
			{
				s.push(node[slen][k]);
				flag[slen][k] = 1;
				slen++;
			}

		}
		else
		{
			// assure that the last point is adjacent to the end point
			from = s.top();
			to = end;
			iLink = sn.GetLinkIndex(from, to);

			if (iLink != -1)
			{
				for (int k = slen - 1, top = s.size() - 1; k >= 0; k--, top--)
				{
					path[k][paths] = s[top];
				}
				paths++;
				if (paths >= MAX_PATHS)
					break;
			}
			s.pop();
			slen--;
		}

	}



	for (int i = 0; i < paths; i++)
		path[0][i] = start;
	for (int i = 0; i < paths; i++)
		path[len][i] = end;

	for (int i = 0; i<path.size(); i++)
	{
		for (int j = 0; j<paths; j++)
		{
			//cout<<path[i][j]<<" ";
		}
		//cout<<""<<endl;
	}
	//getchar();







	for (int i = 0; i < paths; i++)
	{
		node_off = 0;
		int j, k;
		int Flag;
		float maxbw = 999;
		for (j = 0; j < len; j++)
		{
			from = path[j][i];
			to = path[j + 1][i];

			if (from != start &&curSn.currentNode[from].rest_cpu < INTER_NODE_CPU)
			{
				path[0][i] = -1;
				break;
			}
			if (to != end && curSn.currentNode[to].rest_cpu < INTER_NODE_CPU)
			{
				path[0][i] = -1;
				break;
			}
			Flag = 0;


			k = sn.GetLinkIndex(from, to);
			if (k != -1)
			{

				Flag = 1;
				//if(curSn.currentNode[to].state == CurrentNode::NODE_OFF)
				//	node_off ++;
				if (curSn.currentLink[k].rest_bw<maxbw)
					maxbw = curSn.currentLink[k].rest_bw;
				//cout<<"maxbw"<<maxbw<<endl;
			}

		}
		if (j == len)
		{
			path[0][i] = maxbw;
		}
	}

	int min = -1;
	for (int j = 0; j <paths; j++)
	{
		if (path[0][j] != -1)
		{
			min = j;
			break;
		}
	}
	for (int j = 1; j < paths; j++)
	{
		if (path[0][j] != -1 && path[0][j] > path[0][min])
			min = j;
	}


	if (min == -1)
		return 0;

	for (int i = 0; i <= len; i++)
	{
		//cout<<path[i][min]<<endl;
		if (i == 0)
		{
			shortestpath[0] = start;
		}
		else
		{
			shortestpath[i] = path[i][min];
		}
	}

	return 1;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool tpathmapping(SubstrateNetwork &sn, CurrentSn &curSn, int index, vector<int> &shortestpath, float maxbw)
{
	float need_bw = maxbw;
	int len = shortestpath.size();
	for (int j = 0; j<len - 1; j++)
	{
		int from = shortestpath[j];
		int to = shortestpath[j + 1];

		int link_index = sn.GetLinkIndex(from, to);
		if (link_index == -1)
		{
			cout << "no such link" << endl;
			return false;
		}
		if (curSn.currentLink[link_index].rest_bw < need_bw)
		{

			cout << "There is not enough bw " << endl;
			return false;
		}
		int req_count = curSn.currentLink[link_index].link_count;
		curSn.currentLink[link_index].req[req_count] = index;
		curSn.currentLink[link_index].vlink[req_count] = -1;
		curSn.currentLink[link_index].bw[req_count] = need_bw;
		curSn.currentLink[link_index].rest_bw -= need_bw;
		curSn.currentLink[link_index].link_count++;
	}

	return true;
}




bool releaselinkMIC(const SubstrateNetwork &sn, CurrentSn &curSn, int index)
{

	int curIndex = index;

	for (int i = 0; i<sn.subInfo.links; i++)
	{
		for (int j = 0; j<MAX_REQ_LINKS; j++)
		{
			if (curSn.currentLink[i].req[j] == curIndex)
			{
				int link_count = curSn.currentLink[i].link_count;
				if (link_count == 1)
				{
					curSn.currentLink[i].link_count = 0;
					curSn.currentLink[i].rest_bw += curSn.currentLink[i].bw[j];
					curSn.currentLink[i].bw[j] = 0;
					curSn.currentLink[i].vlink[j] = -1;
					curSn.currentLink[i].req[j] = -1;
				}
				else
				{
					curSn.currentLink[i].link_count--;
					curSn.currentLink[i].rest_bw += curSn.currentLink[i].bw[j];

					curSn.currentLink[i].req[j] = curSn.currentLink[i].req[link_count - 1];
					curSn.currentLink[i].vlink[j] = curSn.currentLink[i].vlink[link_count - 1];
					curSn.currentLink[i].bw[j] = curSn.currentLink[i].bw[link_count - 1];

					curSn.currentLink[i].req[link_count - 1] = -1;
					curSn.currentLink[i].vlink[link_count - 1] = -1;
					curSn.currentLink[i].bw[link_count - 1] = 0;
					j--;

				}
			}
		}
	}//end of int i=0; i<sn.subInfo.links

	//just check
	for (int i = 0; i<sn.subInfo.links; i++)
	{
		float rest_bw = curSn.currentLink[i].rest_bw;
		for (int j = 0; j<MAX_REQ_LINKS; j++)
		{
			if (curSn.currentLink[i].req[j] != -1)
			{
				rest_bw += curSn.currentLink[i].bw[j];
			}
		}

		float tmp = sn.link[i].bw;
		if (fabs(rest_bw - tmp) >= 0.01)
		{
			cout << "there is something wrong with releasing link" << endl;
			return false;
		}
	}
	return true;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool unique_finder(class phyNode a, class phyNode b)
{
	return a.physicalNodeIndex == b.physicalNodeIndex;
}

bool sort_finder(phyNode a, phyNode b)
{
	return a.physicalNodeIndex < b.physicalNodeIndex;
}
int comparePhyIndex(class phyNode a, class phyNode b)
{

	return a.physicalNodeIndex > b.physicalNodeIndex;

}
int cmp(int a, int b)
{
	return(a > b) ? 1 : 0;
}
virNode find_min(vector<virNode> a, vector<virNode> b, vector<vector<int>> c)
{
	int p = 0;
	totalW totalW1;
	for (int i = 0; i<b.size(); i++)
	{
		if (b[i].devided == false)
		{
			for (int j = 0; j<a.size(); j++)
			{
				totalW1.TW.push_back(c[i][j]);
			}
		}
		totalW1.tW.push_back(accumulate(totalW1.TW.begin(), totalW1.TW.end(), 0));
		totalW1.index.push_back(i);
	}
	for (int i = 0; i < totalW1.tW.size() - 1; i++)
	{
		if (totalW1.tW[i] <= totalW1.tW[i + 1])
		{
			int totalW_final = totalW1.tW[i];
			p = totalW1.index[i];
		}
	}
	return b[p];
}

bool check(vector<virNode> a, vector<vector<float>> b)
{


	for (int i = 0; i < a.size(); i++)
	{

		int label_destination = a[i].devided_group;

		if (label_destination == findmax(a, i, b))
			continue;
		else
			return 0;
	}
	return 1;

}

int get_max(int a, vector<virNode>b, vector<vector<float>> w)//
{
	vector<vector<float>>label;
	Label label_sum;
	label.clear();
	label.resize(b.size());

	int index;
	for (int i = 0; i < b.size(); i++)
	{
		for (int j = 0; j < w[a].size(); j++)
		{
			if (a == j)continue;
			if ((b[j].devided_group == i) && (w[a][j] != 0))
			{
				label[i].push_back(w[a][j]);

			}
		}

		label_sum.weight.push_back(accumulate(label[i].begin(), label[i].end(), 0));

	}
	//sort(label_sum.weight.begin(), label_sum.weight.end(), cmp);
	vector<float> ::iterator maxweight = max_element(label_sum.weight.begin(), label_sum.weight.end());
	index = distance(label_sum.weight.begin(), maxweight);
	return index;
}
int get_min(int a, vector<virNode>b, vector<vector<float>> w)
{
	int min = 0;
	for (int i = 0; i < b.size(); i++)
	{
		//if (a == i)continue;
		if (w[a][min] < w[a][i])
			min = i;
	}
	return b[min].devided_group;
}
int get_groupNum(vector<virNode>a)
{
	int count;
	//sort(a.begin(), a.end(), sortfind);
	//a.erase(unique(a.begin(), a.end(), uniquefind),a.end());
	//for (int i = 0; i < a.size(); i++)
	//{
	//	if (a[i].devided_group ==i )
	//		 count++;
	//}
	count = a.size();
	return count;
}
bool sortfind(virNode a, virNode b)
{
	return a.devided_group < b.devided_group;
}

bool uniquefind(virNode a, virNode b)
{
	return a.devided_group == b.devided_group;
}
int findmax(vector<virNode>a, int b, vector<vector<float>>w_label)
{
	int index = 0;
	vector<int>c;
	c.clear();
	int group_destination;

	for (int i = 0; i < a.size(); i++)
	{
		if (w_label[b][i] != 0)
		{
			c.push_back(a[i].devided_group);
		}
	}

	for (int i = 0; i < c.size(); i++)
	{
		int temp = 0; int max = 0;
		temp = count(c.begin(), c.end(), c[i]);
		if (temp >= max)
		{
			max = temp;
			index = c[i];
		}
	}
	return index;
}
// TestVN.cpp
//

#include "stdafx.h"
#include "Simulator.h"
#include <queue>  
#include <iostream>
#include"windows.h"
#include"stdlib.h"
#include<numeric>
using namespace std;
int iterate_time;
//the struct used in the priority queue for comparison of Antibodies
struct cmp{
	bool operator()(PSO a, PSO b){
		return a.resource<b.resource;
	}
};

int main(int argc, char* argv[])
{
	time_t start_t, end_t;
	//record the start time of the program
	start_t = clock();
	srand(time(NULL));
	//an instance of SubstrateNetwork
	//it consists the topology of the network, 
	//including the number of nodes and links,
	//its Bandwidth and CPU resource distribution
	//and its path information. For more detail, 
	//please refer to "vnEmbed.h"
	SubstrateNetwork sn;
	//a vector of Class Requests including the 
	//arrival time, duration, and most importantly
	//the topology this request is requesting
	vector<Request> requests;
	//the standard input format includes the number of requests,
	//the directory containing the substrate network info file 
	//the directory containing the request files
	//the time of iteration
	//NOTE: there are other versions of this program that asks for 
	//six parameters, the extra one is the parameter that alters either
	//the duration of each request or their arraival rate. Of course, this
	//alteration can also be achieved by modifying the request info files.
	if (argc != 5)
	{
		cout << "InputFormat:" << endl;
		cout << "./vnEmbed <req_count> <sn_dir> <req_dir> <iterate_time>" << endl;
		getchar();
		system("pause");
		return 0;
	}
	//read substrate network information
	if (!read_sn_info(argv[2], sn))//sn_dir=subtest
	{
		cout << "Read SN info failed " << endl;
		return 0;
	}
	int REQsize = atoi(argv[1]);
	requests.resize(REQsize);
	if (!read_req_info(argv[3], REQsize, requests))
	{
		cout << "Read REQ info failed " << endl;
		return 0;
	}

	//floyd shortest path calculation algorithm
	struct shortest_path *spath;
	spath = (shortest_path*)malloc(sizeof(shortest_path)*sn.subInfo.nodes*sn.subInfo.nodes);
	calc_shortest_path(sn, spath);

	//spath[i][j]= k, k is the next step from i
	//sn represents the sn that we have (static information), currentsn represents the sn that we currently have (dynamic information)

	CurrentSn currentSn(sn);
	CurrentSn copy_currentSn(sn);
	//We simulate the process of requests arriving randomly
	//by this Simulator class. We take out an event each time from it
	Simulator sim(requests);
	//The ordered queue of Antibodies
	//two variables that record the revenue and cost produced so far
	double done_rev = 0;
	double done_cost = 0;
	double switchcost = 0;
	double cost = 0;
	//the variable that records the machines that are turned on
	int onCount = 0;
	//init[] records the first Antibodies we produce before the multiply and mutation 
	PSO init[GROUP];
	int matchflag[2000];

	//int before[50] = { 0 };
	//int after[50] = { 0 };
	//int state[50] = { 0 };
	//int stateNO = 0;
	//antibody2[][] records the Antibodies from "set" that do NOT go through the multiply and mutation
	//which means antibody2[][] preserves the original copy of antibody[][]
	//both of them will be put back to "set" again, from which we will get the best Antibodies
	//either from the modified Antibodies or their original copies
	//This enables us to preserve the "history best"
	//The time for which we allow the initialization of each Antibody to fail
	int find_antibody_time = 0;
	//counting the number of Antibodies in init[]
	int PNumber = 0;
	iterate_time = atoi(argv[4]);
	//output all the information into file
	//every 100 time unit (not real time)
	char logname[100];
	sprintf(logname, "TIMELOG_%s_REQ%d_ITERATE%d.txt", argv[3], REQsize, iterate_time);
	FILE *fp_log = fopen(logname, "w");
	FILE *fplog = fopen("time", "w");
	if (!fplog)
	{
		cout << "time.txt failed" << endl;
		return 0;
	}
	if (!fp_log)
	{
		cout << "fp_log.txt failed" << endl;
		return 0;
	}
	FILE *fp_log1 = fopen("node.txt", "w");
	if (!fp_log1)
	{
		cout << "node.txt failed" << endl;
		return 0;
	}
	fprintf(fp_log, "%s", "   time  revenue           cost              CPU%    BW%  ratio   oncount\n");
	//output the revenue and energy information into file
	//every 4000 time unit (not real time)
	char graphname[100];
	sprintf(graphname, "result_%s_REQ%d_ITERATE%d.txt", argv[3], REQsize, iterate_time);
	FILE *graph = fopen(graphname, "w");
	double restCPU;
	double restBW;
	int finalindex;
	float Maxfitness;
	vector<int> gbest;
	//vector<Widget> vWidgets(500, Widget(0));
	float gbest_resource;
	int gindex;
	int accept = 0;
	int count = 0;
	int reqindex;
	virNode virnode;
	phyNode phynode;
	vector<virNode> virNodequeue;
	virNodequeue.clear();
	vector<virNode>virNodequeue_temp;
	vector<phyNode> physicalqueue;
	vector<phyNode>physicalqueue_temp;
	physicalqueue.clear();
	int vir_index;
	int loop;
	vector< vector<float> > Weight;
	vector< vector<int> > edge;
	vector<float> xn;
	vector<float> yn;

	vector< vector<int> > copy_Weight;
	vector< vector<int> > copy_edge;
	vector<float> copy_xn;
	vector<float> copy_yn;
	vector<int> temp_node_mapping_method;
	vector< vector<int> > temp_link_mapping_method;
	float fcost;
	float bcost;
	float maxW;
	int temp1;
	int alert = 0;
	int error = 0;
	int number = 0;
	int numberon = 0;
	int stateNoOut = 0;
	float totaltime;
	vector<int> visitY;
	vector<int> visitX;
	vector<int> matchXY;
	vector< vector<int> > matchYX;
	int temp;
	int migration = 0;
	vector<int> shortestpath;
	vector<vector<int>> w;
	vector<vector<float>>w_label;
	//w.resize(virNodequeue.size());
	vector<vector<phyNode>>	devide_nodepool;
	vector<vector<int>>virNodequeue_devided;
	totalW totalW1;//
	long int timer = 0;
	int  constent = 5;
	int K;
	vector<vector<int>>physical_groupnodeindex;
	vector<vector<phyNode>>physical_devided;
	vector<vector<virNode>>virNode_devided;
	/*
	for(int i=0;i<50;i++)
	cout<<currentSn.currentNode[i].rest_memory<<endl;
	getchar();
	for(int i=0;i<10;i++)
	for(int j=0;j<requests[i].node.size();j++)
	cout<<requests[i].node[j].memory<<endl;
	getchar();
	*/
	fprintf(graph, "%s", "   time  revenue           cost              CPU%    BW%     ratio		oncount\n");

	int lasttime = 0;
	while (!sim.empty())
	{

		const Event &curEvent = sim.top();//get the top event
		int curIndex = curEvent.index;//get the index of current event 
		vector<int> node_mapping_method;
		vector< vector<int> > link_mapping_method;
		vector< vector<int> > nodepool;
		int nodesOfCurReq = 0;
		int linksOfCurReq = 0;
		int timeduration = 0;
		//update the variable "oncount", and calculate done_cost
		//cout<<"curEvent.time"<<curEvent.time<<endl;
		costUpdate(sn, currentSn, curEvent.time, done_cost, onCount);

		/*
		timeduration=curEvent.time-lasttime;
		for(int i=0;i<50;i++)
		{
		if(currentSn.currentNode[i].state == CurrentNode::NODE_ON)
		currentSn.currentNode[i].timeduration+=timeduration;
		}
		lasttime=curEvent.time;
		*/
		//cout<<"curEvent.type"<<curEvent.type<<endl;
		switch (curEvent.type){
			//Logging event
		case EVENT_LOG:
			if (curEvent.time <= 48000)
			{
				restCPU = 0;
				restBW = 0;
				//calculate the remaining resource of the substrate network
				for (int i = 0; i<sn.subInfo.nodes; i++)
				{
					restCPU += currentSn.currentNode[i].rest_cpu;
				}
				for (int i = 0; i<sn.subInfo.links; i++)
				{
					restBW += currentSn.currentLink[i].rest_bw;
				}
				fprintf(fp_log, "%6d   %-16lf   %-16lf   %-5.2lf   %-5.2lf  %-2lf  %-2d\n",
					curEvent.time, done_rev / curEvent.time, (done_cost + switchcost) / curEvent.time, (sn.totalCPU - restCPU + 0.0001) / sn.totalCPU * 100,
					(sn.totalBW - restBW + 0.0001) / sn.totalBW * 100, (accept*1.0) / (count*1.0), onCount);
				fflush(fp_log);
				if ((curEvent.time % 4000 == 0) && (curEvent.time>0))
				{
					fprintf(graph, "%6d   %-16lf   %-16lf %-5.2lf %-5.2lf  %-2lf  %-2d\n", curEvent.time, done_rev / curEvent.time, done_cost / curEvent.time, (sn.totalCPU - restCPU + 0.0001) / sn.totalCPU * 100, (sn.totalBW - restBW + 0.0001) / sn.totalBW * 100, (accept*1.0) / (count*1.0), onCount);
				}
				fflush(graph);
			}
			break;
			///////////////////here comes a request, get started!////////////////////
		case  EVENT_ARRIVE:
			count++;
			find_antibody_time = 0;
			nodesOfCurReq = requests[curIndex].nodes;//the number of nodes of current request
			cout << "REQ: " << curIndex << endl;
			/////////////here starts the initialization process!////////////////
			for (PNumber = 0; find_antibody_time<20 && PNumber<GROUP; PNumber++)
			{
				node_mapping_method.clear();
				node_mapping_method.resize(nodesOfCurReq);
				nodepool.clear();
				nodepool.resize(nodesOfCurReq);
				//node_mapping_method[i] is the number of node in sn for node i in request
				if (!(node_mapping(sn, currentSn, requests[curIndex], node_mapping_method, nodepool)))
				{
					nodepool.clear();
					//cout<<"Node_mapping_first failed"<<endl;
					find_antibody_time++;
					PNumber--;
					continue;
				}
				linksOfCurReq = requests[curIndex].links;
				link_mapping_method.clear();
				link_mapping_method.resize(linksOfCurReq);
				//link_mapping_method[i][j] is the jth step of the ith sp
				if (!(link_mapping(sn, currentSn, requests[curIndex], node_mapping_method, link_mapping_method, spath)))
				{
					//cout<<"link mapping first failed"<<endl;
					find_antibody_time++;
					for (int i = 0; i<sn.subInfo.links; i++)
					{
						currentSn.currentLink[i].pre_used = 0;
					}
					PNumber--;
					continue;
				}
				//clear pre_used
				for (int i = 0; i<sn.subInfo.links; i++)
				{
					currentSn.currentLink[i].pre_used = 0;
				}
				//save this mapping method
				init[PNumber].node_mapping_method = node_mapping_method;
				init[PNumber].pbest = node_mapping_method;

				for (int i = 0; i<linksOfCurReq; i++)
				{
					init[PNumber].resource += requests[curIndex].link[i].bw*(link_mapping_method[i].size() - 1);
				}
				//if(init[FNumber].resource<MIN_Resource)MIN_Resource=init[FNumber].resource;
				//if(init[FNumber].resource>MAX_Resource)MAX_Resource=init[FNumber].resource;
				init[PNumber].speed.clear();
				init[PNumber].speed.resize(nodesOfCurReq);
				for (int i = 0; i<nodesOfCurReq; i++)
				{
					if (((double)rand() / ((double)(RAND_MAX)+(double)(1)))<0.5)
						init[PNumber].speed[i] = 0;
					else
						init[PNumber].speed[i] = 1;
				}
				for (int j = 0; j<linksOfCurReq; j++)
				{
					init[PNumber].resource += requests[curIndex].link[j].bw*(link_mapping_method[j].size() - 1);
				}
				init[PNumber].best_resource = init[PNumber].resource;
			}//end of for(PNumber = 0）
			if (PNumber<5)
			{
				cout << "Error! Initialization failed!" << endl;
				break;
			}
			/*

			for(int j=0;j<PNumber;j++)
			{
			for(int k=0;k<init[j].node_mapping_method.size();k++)
			{
			cout<<init[j].node_mapping_method[k]<<" ";
			if(k==(init[j].node_mapping_method.size()-1))
			cout<<"------------------------------"<<endl;
			}
			}

			for(int j=0;j<PNumber;j++)
			{
			for(int k=0;k<init[j].speed.size();k++)
			{
			cout<<init[j].speed[k]<<" ";
			if(k==(init[j].speed.size()-1))
			cout<<"------------------------------"<<endl;
			}
			}
			*/
			//getchar();
			gbest_resource = init[0].best_resource;
			for (int i = 0; i<PNumber; i++)
			{
				if (init[i].best_resource <= gbest_resource)
				{
					gbest_resource = init[i].best_resource;
					finalindex = i;
				}
			}
			gbest = init[finalindex].pbest;
			////////////////here comes the iteration!!!//////////////// 
			//////////////////////////////////////////////////////////////////////////////////////////////////////////
			for (int iter = 0; iter<iterate_time; iter++)
			{

				double r = 0;
				vector<int> minus_pNMM;
				vector<int> minus_gNMM;
				vector<int> mutiple_NMM;
				minus_pNMM.clear();
				minus_gNMM.clear();
				mutiple_NMM.clear();
				minus_pNMM.resize(nodesOfCurReq);
				minus_gNMM.resize(nodesOfCurReq);
				mutiple_NMM.resize(nodesOfCurReq);
				for (int i = 0; i<PNumber; i++)
				{
					for (int j = 0; j<init[i].speed.size(); j++)
					{
						if (init[i].pbest[j] == init[i].node_mapping_method[j])
							minus_pNMM[j] = 1;
						else
							minus_pNMM[j] = 0;
						if (gbest[j] == init[i].node_mapping_method[j])
							minus_gNMM[j] = 1;
						else
							minus_gNMM[j] = 0;
						r = ((double)rand() / ((double)(RAND_MAX)+(double)(1)));
						if (r<P1) init[i].speed[j] = init[i].speed[j];
						if (P1<r<(P1 + P2)) init[i].speed[j] = minus_pNMM[j];
						if ((P1 + P2)<r) init[i].speed[j] = minus_gNMM[j];
					}
					bool temp;
					do
					{
						temp = true;
						for (int j = 0; j<init[i].speed.size(); j++)
						{
							if (init[i].speed[j] == 0)
							{
								int dim;
								r = ((double)rand() / ((double)(RAND_MAX)+(double)(1)));
								dim = (int)(r*nodepool[j].size());
								init[i].node_mapping_method[j] = nodepool[j][dim];
							}
						}
						for (int j = 0; j <init[i].node_mapping_method.size(); j++)
						{
							for (int k = 0; k < init[i].node_mapping_method.size(); k++)
							{
								if (j != k && init[i].node_mapping_method[j] == init[i].node_mapping_method[k])//loop, error
								{
									temp = false;
									// cout<<"there is a loop, restart"<<endl;         
								}
							}
						}
					} while (temp == false);
					int temp1 = true;
					link_mapping_method.clear();
					link_mapping_method.resize(linksOfCurReq);
					if (!(link_mapping(sn, currentSn, requests[curIndex], init[i].node_mapping_method, link_mapping_method, spath)))
					{
						for (int n = 0; n<sn.subInfo.links; n++)
						{
							currentSn.currentLink[n].pre_used = 0;
						}
						init[i].resource = 999999;
						temp1 = false;
					}
					for (int n = 0; n<sn.subInfo.links; n++)
					{
						currentSn.currentLink[n].pre_used = 0;
					}
					if (temp1 == true)
					{
						init[i].resource = 0;
						for (int j = 0; j<linksOfCurReq; j++)
						{
							init[i].resource += requests[curIndex].link[j].bw*(link_mapping_method[j].size() - 1);
						}
					}
					if (init[i].resource<init[i].best_resource)
					{
						init[i].best_resource = init[i].resource;
						init[i].pbest = init[i].node_mapping_method;
					}
					gbest_resource = init[0].best_resource;
					for (int j = 0; j<PNumber; j++)
					{
						if (init[j].best_resource <= gbest_resource)
						{
							gbest_resource = init[j].best_resource;
							gindex = j;
						}
					}
					gbest = init[gindex].pbest;
				}//end of for(int i=0;i<PNumber;i++)
			}
			link_mapping_method.clear();
			link_mapping_method.resize(linksOfCurReq);
			if (!link_mapping(sn, currentSn, requests[curIndex], gbest, link_mapping_method, spath))
			{
				//cout<<"link_mapping final failed!"<<endl;
				//exit(0);
				break;
			}
			if (!(do_mapping(sn, currentSn, requests[curIndex], gbest, link_mapping_method)))
			{
				//cout<<"do mapping final failed"<<endl;
				break;
			}
			cost += gbest_resource;
			for (int j = 0; j<requests[curIndex].node.size(); j++)
			{
				cost += requests[curIndex].node[j].cpu;
			}
			sim.push(Event(EVENT_DEPART, requests[curIndex].time + requests[curIndex].duration, curIndex, -1, -1, -1));
			accept++;
			matchflag[curIndex] = 1;
			node_mapping_method.clear();
			link_mapping_method.clear();
			break;
		case EVENT_DEPART:
			done_rev += requests[curIndex].GetRevenue();
			//release the resources
			release_cpu(sn, currentSn, requests[curIndex], sn.subInfo.nodes);
			release_link(sn, currentSn, requests[curIndex], sn.subInfo.links);
			break;
		case EVENT_MIGRATION:
			time_t StartClock, EndClock;
			int before[50];
			int after[50];
			int state[50];
			int stateNO;
			//before[50] = { 0 };
			memset(before, 0, sizeof(before));
			memset(after, 0, sizeof(after));
			memset(state, 0, sizeof(state));

			stateNO = 0;
			StartClock = clock();
			physicalqueue.clear();
			virNodequeue.clear();
			physicalqueue_temp.clear();
			//cout<<"migration"<<endl;
			if (curEvent.time <= 48000)
			{

				//getchar();
				// for(int i=0; i<sn.subInfo.nodes; i++)
				//  {
				// if(currentSn.currentNode[i].state == CurrentNode::NODE_ON)
				//  {
				//   onMichine++;
				//  } 0
				//   if(currentSn.currentNode[i].state == CurrentNode::NODE_OFF)
				//   {
				//  offMichine++;    
				// cout<<i<<endl;
				//   }
				//   }

				for (int i = 0; i < sn.subInfo.nodes; i++)
				{
					currentSn.currentNode[i].interNodeConstrain = 0;
					for (int j = 0; j < currentSn.currentNode[i].req_count; j++)
					{
						if (currentSn.currentNode[i].vnode[j] == -1)
						{
							currentSn.currentNode[i].interNodeConstrain = 1;

							break;
						}
					}
					if (currentSn.currentNode[i].state == CurrentNode::NODE_ON)
					{
						numberon++;
						before[i] = 1;
						cout << i << endl;
					}

					if ((currentSn.currentNode[i].state == CurrentNode::NODE_ON) && ((sn.cpu[i]) - currentSn.currentNode[i].rest_cpu <= (0.2*sn.cpu[i])) && (currentSn.currentNode[i].interNodeConstrain == 0) && (currentSn.currentNode[i].hot != 1) || ((currentSn.currentNode[i].hot == -1) && (currentSn.currentNode[i].state == CurrentNode::NODE_ON) && (currentSn.currentNode[i].interNodeConstrain == 0)))//节点i是待迁移节点
					{
						number++;
						for (int j = 0; j < currentSn.currentNode[i].req_count; j++)
						{
							virnode.totalcpu += currentSn.currentNode[i].cpu[j];
						}
						//cout<<"totalCPU"<<virnode.totalcpu<<endl;
						for (int j = 0; j<currentSn.currentNode[i].req_count; j++)
						{
							virnode.snode_index = i;
							virnode.jindex = j;
							virnode.req_index = currentSn.currentNode[i].req[j];
							virnode.vir_index = currentSn.currentNode[i].vnode[j];
							virnode.cpu = currentSn.currentNode[i].cpu[j];
							virnode.memory = currentSn.currentNode[i].memory[j];

							virNodequeue.push_back(virnode);
						}
						virnode.totalcpu = 0;
					}
					if ((currentSn.currentNode[i].state == CurrentNode::NODE_ON) && ((sn.cpu[i]) - currentSn.currentNode[i].rest_cpu>(0.2*sn.cpu[i])) && (currentSn.currentNode[i].hot != -1) || currentSn.currentNode[i].hot == 1)//目标物理节点
					{
						phynode.physicalNodeIndex = i;
						phynode.orireq_count = phynode.req_count = currentSn.currentNode[i].req_count;
						phynode.rest_cpu = currentSn.currentNode[i].rest_cpu;
						phynode.rest_memory = currentSn.currentNode[i].rest_memory;
						phynode.hot = currentSn.currentNode[i].hot;
						for (int j = 0; j < currentSn.currentNode[i].req_count; j++)
						{
							phynode.req[j] = currentSn.currentNode[i].req[j];
							phynode.vnode[j] = currentSn.currentNode[i].vnode[j];
							phynode.cpu[j] = currentSn.currentNode[i].cpu[j];
							phynode.memory[j] = currentSn.currentNode[i].memory[j];
						}

						physicalqueue.push_back(phynode);

						for (int j = 0; j < currentSn.currentNode[i].req_count; j++)
						{
							phynode.req[j] = 0;
							phynode.vnode[j] = 0;
							phynode.cpu[j] = 0;
							phynode.memory[j] = 0;
						}
					}
					currentSn.currentNode[i].interNodeConstrain = 0;
				}// for(int i=0; i<sn.subInfo.nodes; i++)

				virNodequeue_temp.clear();
				physicalqueue_temp = physicalqueue;
				virNodequeue_temp = virNodequeue;

				w.clear();
				w.resize(virNodequeue.size());
				w_label.clear();
				w_label.resize(virNodequeue.size());

				virNode_devided.clear();
				//virNode_devided.resize(constent);
				physical_devided.clear();
				//physical_devided.resize(constent);
				physical_groupnodeindex.clear();
				//physical_groupnodeindex.resize(constent);
				virNodequeue_devided.clear();
				//virNodequeue_devided.resize(constent);

				totalW1.TW.clear();
				totalW1.tW.clear();//jia
				totalW1.index.clear();


				devide_nodepool.clear();
				devide_nodepool.resize(virNodequeue.size());

				for (int i = 0; i < virNodequeue.size(); i++)
				{
					devide_nodepool[i].clear();
					for (int j = 0; j < physicalqueue_temp.size(); j++)
					{

						if ((physicalqueue_temp[j].rest_cpu >= virNodequeue[i].cpu) && (physicalqueue_temp[j].rest_memory >= virNodequeue[i].memory))
							//	candidate_node a[i] = physicalqueue[j];
						{
							devide_nodepool[i].push_back(physicalqueue_temp[j]);

						}
					}


					sort(devide_nodepool[i].begin(), devide_nodepool[i].end(), comparephyCPU_reverse);//从大到小排序

				}



				for (int i = 0; i < virNodequeue.size(); i++)
				{
					w[i].resize(virNodequeue.size());
					for (int j = 0; j < virNodequeue.size(); j++)
					{
						//cout << currentSn.currentNode[j].cpu << endl;
						w[i][j] = compare(devide_nodepool[i], devide_nodepool[j]);
						//cout<<devide_nodepool[i][0].devided<<endl;
						//cout<<devide_nodepool[j]<<endl;
					}
				}

				for (int i = 0; i < virNodequeue.size(); i++)
				{
					w_label[i].resize(virNodequeue.size());
					for (int j = 0; j < virNodequeue.size(); j++)
					{
						w_label[i][j] = fabs(virNodequeue[i].cpu - virNodequeue[j].cpu) + fabs(virNodequeue[i].memory - virNodequeue[j].memory);
						if (w_label[i][j] <= 16)w_label[i][j] = 0;
						//cout << i << "to" << j << ":" << w_label[i][j] << endl;
						//if (w_label[i][j] <= 4)w_label[i][j] = 0;
					}
				}

				//for (int i = 0; i < virNodequeue.size(); i++)
				//{
				//	w_label.resize(virNodequeue.size());
				//	for (int j = 0; j < virNodequeue.size(); j++)
				//	{
				//		if ((w[i][j] == 0) && (i != j))
				//		{
				//			w_label[i].push_back(9999);
				//		}
				//		else if ((!(w[i][j]>3))&&(w[i][j]!=0))
				//			w_label[i].push_back(1 / w[i][j]);
				//		else
				//			continue;
				//	}
				//}


				const int maxt = 3;
				for (int i = 0; i < virNodequeue.size(); i++)
				{
					virNodequeue[i].devided_group = i;
					//cout << "virNodequeue[" << i << "].devided_group:" <<i<< endl;
				}

				int t = 0;
				while (true)
				{
					if ((check(virNodequeue, w_label) == 1) || (t == maxt))
					{
						break;
					}
					else
					{
						t = t + 1;

						for (int i = 0; i < virNodequeue.size(); i++) 
						{
						
							virNodequeue[i].devided_group = get_max(i, virNodequeue, w_label);
							
						}
						//break;
					}
				}
				virNodequeue_temp = virNodequeue;
				sort(virNodequeue_temp.begin(), virNodequeue_temp.end(), sortfind);
				virNodequeue_temp.erase(unique(virNodequeue_temp.begin(), virNodequeue_temp.end(), uniquefind), virNodequeue_temp.end());

				constent = get_groupNum(virNodequeue_temp);
				physical_groupnodeindex.resize(constent);
				virNodequeue_devided.resize(constent);

				virNode_devided.resize(constent);

				for (int i = 0; i < constent; i++)
				{
					for (int j = 0; j < virNodequeue.size(); j++)
					{
						if (virNodequeue[j].devided_group == virNodequeue_temp[i].devided_group)
							virNode_devided[i].push_back(virNodequeue[j]);
					}
				}

				physical_devided.resize(constent);
				for (int k = 0; k < constent; k++)//
				{
					for (int j = 0; j < virNodequeue.size(); j++)
					{
						if (virNodequeue[j].devided_group == virNodequeue_temp[k].devided_group)//
						{
							for (int d = 0; d < devide_nodepool[j].size(); d++)
							{
								if (devide_nodepool[j][d].devided == false)
								{
									physical_groupnodeindex[k].push_back(devide_nodepool[j][d].physicalNodeIndex);//
									physical_devided[k].push_back(devide_nodepool[j][d]);
									devide_nodepool[j][d].devided = true;
								}
							}

						}
					}

					sort(physical_devided[k].begin(), physical_devided[k].end(), sort_finder);
					physical_devided[k].erase(unique(physical_devided[k].begin(), physical_devided[k].end(), unique_finder), physical_devided[k].end());

				}//end of for(int k=0,k<constent;k++)



				number = 0;
				numberon = 0;
				for (int group = 0; group < constent; group++)
				{
					sort(physical_devided[group].begin(), physical_devided[group].end(), comparephyCPU);
					sort(virNode_devided[group].begin(), virNode_devided[group].end(), comparevirCPU);


					for (int i = 0; i < physical_devided[group].size(); i++)
					{
						for (int j = 0; j < virNode_devided[group].size(); j++)
						{
							if (physical_devided[group][i].physicalNodeIndex == virNode_devided[group][j].snode_index)
							{
								cout << "这是很严重的错误啊" << endl;
								getchar();
							}
						}
					}
				}//end of for (int k = 0; k < constent; k++)

			
					for (int group = 0; group < constent; group++)
					{
						Weight.clear();
						copy_Weight.clear();
						Weight.resize(virNode_devided[group].size());
						copy_Weight.resize(virNode_devided[group].size());
						edge.clear();
						copy_edge.clear();
						edge.resize(virNode_devided[group].size());
						copy_edge.resize(virNode_devided[group].size());
						xn.clear();
						copy_xn.clear();
						yn.clear();
						copy_yn.clear();
						xn.resize(virNode_devided[group].size());
						yn.resize(physical_devided[group].size());
						copy_xn.resize(virNode_devided[group].size());
						copy_yn.resize(physical_devided[group].size());



						for (int u = 0; u < virNode_devided[group].size(); u++)
						{
							reqindex = virNode_devided[group][u].req_index;
							vir_index = virNode_devided[group][u].vir_index;
							for (int v = 0; v < physical_devided[group].size(); v++)
							{
								temp_node_mapping_method.clear();
								temp_link_mapping_method.clear();

								//getchar();
								temp_node_mapping_method.resize(requests[reqindex].node_mapping_method.size());
								temp_node_mapping_method = requests[reqindex].node_mapping_method;

								temp_node_mapping_method[vir_index] = physical_devided[group][v].physicalNodeIndex;
								//	for (int i = 0; i < temp_node_mapping_method.size(); i++)
								//cout<<temp_node_mapping_method[i]<<" ";
								//cout<<"一"<<endl;
								//getchar();
								temp_link_mapping_method.resize(requests[reqindex].link_mapping_method.size());
								temp_link_mapping_method = requests[reqindex].link_mapping_method;


								loop = 0;
								for (int i = 0; i < temp_node_mapping_method.size(); i++)
								{
									if ((i != vir_index) && (temp_node_mapping_method[vir_index] == temp_node_mapping_method[i]))
									{
										loop = 1;
										break;
									}
								}
								//cout<<"loop"<<loop<<endl;

								for (int i = 0; i < requests[reqindex].links; i++)
								{
									int req_from = requests[reqindex].link[i].from;
									int req_to = requests[reqindex].link[i].to;
									if ((vir_index == req_from) || (vir_index == req_to))
									{
										temp_link_mapping_method[i].clear();
									}
								}


								//getchar();
								if (loop == 0)
								{
									if ((physical_devided[group][v].rest_cpu>virNodequeue[u].cpu) && (physical_devided[group][v].rest_memory >= virNode_devided[group][u].memory) && (link_mapping_migration(sn, copy_currentSn, requests[reqindex], temp_node_mapping_method, temp_link_mapping_method, spath, vir_index)))//这里没必要使用checksameVnode(没有这个函数)这个函数吧,link_mapping_migration是重点
									{

										for (int i = 0; i < sn.subInfo.links; i++)
										{
											copy_currentSn.currentLink[i].pre_used = 0;
										}

										for (int i = 0; i < requests[reqindex].links; i++)
										{
											fcost += requests[reqindex].link[i].bw*(requests[reqindex].link_mapping_method[i].size() - 1);
											bcost += requests[reqindex].link[i].bw*(temp_link_mapping_method[i].size() - 1);
										}
										if (physical_devided[group][v].hot == 1)
										{
											Weight[u].push_back(fcost - bcost);
											copy_Weight[u].push_back(fcost - bcost);
										}
										if (physical_devided[group][v].hot == 0)
										{
											Weight[u].push_back(fcost - bcost - 1000);
											copy_Weight[u].push_back(fcost - bcost - 1000);
										}
										//cout<<u<<" "<<v<<" "<<fcost-bcost<<" "<<Weight[u][v]<<endl;
										//cout<<"matched"<<endl;
										fcost = 0;
										bcost = 0;
									}//if ((physicalqueue[v].rest_cpu>virNodequeue[u].cpu) && (physicalqueue[v].rest_memory >= virNodequeue[u].memory) && (link_mapping_migration(sn, copy_currentSn, requests[reqindex], temp_node_mapping_method, temp_link_mapping_method, spath, vir_index))
									else
									{
										for (int i = 0; i < sn.subInfo.links; i++)
										{
											copy_currentSn.currentLink[i].pre_used = 0;
										}
										Weight[u].push_back(NOTEXIST);
										copy_Weight[u].push_back(NOTEXIST);
										//cout<<u<<" "<<v<<" "<<NOTEXIST<<" "<<Weight[u][v]<<endl;

									}
								}//if(loop=0)
								else
								{
									Weight[u].push_back(NOTEXIST);
									copy_Weight[u].push_back(NOTEXIST);
									//cout<<u<<" "<<v<<" "<<NOTEXIST<<" "<<Weight[u][v]<<endl;
									//cout<<"unmatchedloop=1"<<endl;
								}
							}//for(int v=0;v<physicalqueue.size();v++)
						}//for(int u=0;u<virNode_devided[j].size();u++)
						//}//end for(int k=0;k<constent;k++)



						//{
						for (int u = 0; u < virNode_devided[group].size(); u++)
						{
							maxW = NOTEXIST;
							for (int v = 0; v < physical_devided[group].size(); v++)
							{
								if (Weight[u][v] >= maxW)
								{
									maxW = Weight[u][v];
									temp1 = v;
								}
							}
							for (int v = 0; v < physical_devided[group].size(); v++)
							{
								edge[u].push_back(0);
								copy_edge[u].push_back(0);
							}
							if (maxW > NOTEXIST)
							{
								copy_edge[u][temp1] = edge[u][temp1] = 1;
								copy_xn[u] = xn[u] = maxW;
							}
						}
						//}

						for (int v = 0; v < physical_devided[group].size(); v++)
						{
							copy_yn[v] = yn[v] = 0;
						}


						visitX.clear();
						visitX.resize(virNode_devided[group].size());
						visitY.clear();
						visitY.resize(physical_devided[group].size());

						for (int i = 0; i < virNode_devided[group].size(); i++)
						{
							visitX.push_back(0);
						}
						for (int i = 0; i < physical_devided[group].size(); i++)
						{
							visitY.push_back(0);
						}

						matchXY.clear();
						matchXY.resize(virNode_devided[group].size());
						matchYX.clear();
						matchYX.resize(physical_devided[group].size());
						for (int u = 0; u < virNode_devided[group].size(); u++)
						{
							matchXY[u] = -1;
						}

						//KM（）
						cout << "start KM" << endl;
						for (int u = 0; u < virNode_devided[group].size(); u++)
						{
							int p = 0;
							while (p < 3)
							{
								visitX.clear();
								visitY.clear();

								for (int i = 0; i < virNode_devided[group].size(); i++)
								{
									visitX.push_back(0);
								}
								for (int i = 0; i < physical_devided[group].size(); i++)
								{
									visitY.push_back(0);
								}

								cout << "start augmentpath" << endl;
								if (find_augument_path(u, virNode_devided[group], visitX, visitY, edge, matchXY, matchYX, physical_devided[group], sn, copy_currentSn, requests, spath))
								{

									for (int i = 0; i < virNode_devided[group].size(); i++)
									{
										cout << matchXY[i] << " ";
									}
									cout << "" << endl;
									//getchar();

									break;
								}
								float slack = 9999999;
								int tempi;
								int tempj;
								for (int i = 0; i < virNode_devided[group].size(); i++)
								{
									if (visitX[i])
									{
										for (int j = 0; j < physical_devided[group].size(); j++)
										{
											if ((!visitY[j]) && ((xn[i] + yn[j] - Weight[i][j]) < slack))
											{
												slack = xn[i] + yn[j] - Weight[i][j];
												tempi = i;
												tempj = j;
											}
										}
									}
								}

								for (int i = 0; i < virNode_devided[group].size(); i++)
								{
									if (visitX[i])
										xn[i] -= slack;
								}
								for (int j = 0; j < physical_devided[group].size(); j++)
								{
									if (visitY[j])
										yn[j] += slack;
								}

								edge[tempi][tempj] = 1;

								p++;
							}
						}
						//KM()


						for (int i = 0; i < virNode_devided[group].size(); i++)
						{
							cout << matchXY[i] << " ";
						}
						cout << "" << endl;
						cout << "up matchXY,outside of KM()" << endl;

						for (int u = 0; u < matchXY.size(); u++)
						{
							if (matchXY[u] != -1)
							{
								vector<int> final_node_mapping_method;
								vector<vector<int> > final_link_mapping_method;
								int reqindex = virNode_devided[group][u].req_index;
								int vir_index = virNode_devided[group][u].vir_index;
								int onstatenumber = 0;
								int ANOTHERonstatenumber = 0;
								final_node_mapping_method.clear();
								final_node_mapping_method.resize(requests[reqindex].node_mapping_method.size());
								final_node_mapping_method = requests[reqindex].node_mapping_method;
								final_node_mapping_method[vir_index] = physical_devided[group][matchXY[u]].physicalNodeIndex;
								final_link_mapping_method.clear();
								final_link_mapping_method.resize(requests[reqindex].link_mapping_method.size());
								for (int i = 0; i < 50; i++)
								{
									if (currentSn.currentNode[i].state == CurrentNode::NODE_ON)
										onstatenumber++;

								}
								release_cpu(sn, currentSn, requests[reqindex], sn.subInfo.nodes);
								release_link(sn, currentSn, requests[reqindex], sn.subInfo.links);
								error++;
								if (link_mapping(sn, currentSn, requests[reqindex], final_node_mapping_method, final_link_mapping_method, spath))
								{
									if (do_mapping(sn, currentSn, requests[reqindex], final_node_mapping_method, final_link_mapping_method))
									{
										requests[reqindex].node_mapping_method.clear();
										requests[reqindex].node_mapping_method = final_node_mapping_method;
										migration++;
									}
								}
								for (int i = 0; i < 50; i++)
								{
									if (currentSn.currentNode[i].state == CurrentNode::NODE_ON)
									{
										ANOTHERonstatenumber++;
										after[i] = 1;
										cout << i << endl;
									}
								}
								if (onstatenumber > ANOTHERonstatenumber)
								{
									alert++;
								}
								for (int i = 0; i < sn.subInfo.links; i++)
								{
									currentSn.currentLink[i].pre_used = 0;
								}
							}
						}


						for (int ii = 0; ii < 50; ii++)
						{
							//cout << before[ii] << " ";
							state[ii] = after[ii] - before[ii];
						}
						cout << endl;
						for (int jj = 0; jj < 50; jj++)
						{
							//cout << after[jj] << " ";
							//cout << endl;
							//state[jj] = after[jj] - before[jj];
							if (state[jj] == 1)
								stateNO++;
							//cout << stateNO << endl;
						}
						cout << endl;
						stateNoOut = stateNO;
						switchcost += 300 * stateNO * 100 / 60;
						vector<int> shortestpath;
						totaltime = 0;
						for (int u = 0; u < matchXY.size(); u++)
						{
							shortestpath.clear();
							if (matchXY[u] != -1)
							{
								if (find_shortestpath_with_min_offnode(virNode_devided[group][u].snode_index, physical_devided[group][matchXY[u]].physicalNodeIndex, spath, sn, shortestpath, currentSn) == 1)
								{
									//for (int i = 0; i <= (shortestpath.size() - 1); i++)
									//{
									//	cout << shortestpath[i] << endl;
									//}
									float maxbw = 999;
									int k;
									for (int j = 0; j < shortestpath.size() - 1; j++)
									{
										k = sn.GetLinkIndex(shortestpath[j], shortestpath[j + 1]);
										if (k != -1)
										{

											if (currentSn.currentLink[k].rest_bw < maxbw)
												maxbw = currentSn.currentLink[k].rest_bw;
										}
									}

									totaltime += 8 * virNode_devided[group][u].memory / maxbw;
									sim.push(Event(EVENT_MIGCOM, curEvent.time + 8 * virNode_devided[group][u].memory / maxbw, 2000 + virNode_devided[group][u].req_index + u * 10000, virNode_devided[group][u].snode_index, physicalqueue[matchXY[u]].physicalNodeIndex, virNode_devided[group][u].req_index));//修改过后有constent条记录

									tpathmapping(sn, currentSn, 2000 + virNode_devided[group][u].req_index + u * 10000, shortestpath, (maxbw - 1));

								}
							}
						}
					}
				fprintf(fplog, "%6d %lf\n", curEvent.time, totaltime);
				fflush(fplog);
				EndClock = clock();
				timer += EndClock - StartClock;
			}//if(curEvent.time<=48000)

			break;
		case EVENT_MIGCOM:
			cout << "before" << endl;
			releaselinkMIC(sn, currentSn, curEvent.index);
			cout << "after" << endl;
			break;
		default:break;
		}//end switch(curEvent.type)
		sim.pop();
	}//end while(!sim.empty())
	double all_reve = 0;
	for (unsigned int i = 0; i<requests.size(); i++)
	{
		if (matchflag[i] == 1)
			all_reve += requests[i].GetRevenue();
	}
	cout << "ITE:" << iterate_time << " All_revenue: " << all_reve << endl;
	cout << "ALL_energy_cost:" << done_cost + switchcost << endl;
	cout << "switch cost:" << switchcost << endl;
	cout << "Accepted: " << accept << endl;
	end_t = clock();
	cout << "Total time is " << end_t - start_t << " seconds" << endl << endl << endl;
	cout << "migration time is " << timer << endl;
	cout << "migration" << migration << "alert" << error << endl;
	fclose(fp_log);
	fclose(graph);
	fclose(fplog);
	getchar();
	system("pause");
	return 0;
}


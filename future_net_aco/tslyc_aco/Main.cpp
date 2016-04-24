#include "io.h"
#include "myFunc.h"


int main(int argc,char *argv[])
{
	cout<<"Begin"<<endl;
	Topo_Node* Head = Node_Init(NODE_MAX_NUM);            //初始化邻接表
	
	Info_Node info;                                   //信息点
	char *topo[5000];                                 //边的数据
	int edge_num;                                     //边的数量
	int demand_num;
	char *demand; //必经点数据
	double init_prob=0;
	int iteration = 0;

	char *topo_file = argv[1];
	
    edge_num = read_file(topo, 5000, topo_file);
    if (edge_num == 0)
    {
        printf("Please input valid topo file.\n");
        return -1;
    }
    char *demand_file = argv[2];
    demand_num = read_file(&demand, 1, demand_file);
    if (demand_num != 1)
    {
        printf("Please input valid demand file.\n");
        return -1;
    }
	srand(time(0));
	Search_Init(Head,info,topo,edge_num,demand);      //搜索初始化

	//DisplayTopo(Head,info);
	//DisplayDemand(info);       
	cout<<info.antNumber<<endl;
	A_Node Ants = Ant_Init(info);                     //蚁群初始化
	ofstream outfile;
	outfile.open("result.txt");
	while (iteration < 5){
		Build_Solution(Head,Ants,info,init_prob);
		Check_Best_Distance(Ants,info);
		Calculate_Fitness(Ants,info);
		Pheromone_Evaporates(Head,info);
		Update_Pheromone(Head,Ants,info);
		Print_Route(Ants,info,outfile,Head);
		iteration++;
	}
	DisplayBestWay(info);
	for(int i=0;i<5;i++){
		Update_Best_Route(Head,Ants,info);
	}
	outfile.close();
	DisplayBestWay(info);
	/*
	FreeAnts(Ants,info);
	FreeInfo(info);
	FreeHead(Head);
	release_buff(topo, edge_num);
    release_buff(&demand, 1);*/
    cout<<"END"<<endl;
	system("pause");
    return 0;
}

#ifndef _MY_FUNC_H_
#define _MY_FUNC_H_

#include <iostream>
#include <string>
#include <map>
#include <cstdlib> 
#include <sstream>
#include <cmath>
#include <ctime>
#include <fstream>
using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//宏定义
//bit位定义
#define BITX        0x0000
#define BIT0        0x0001          //bit位
#define BIT1        0x0002
#define BIT2        0x0004
#define BIT3        0x0008          //bit位
#define BIT4        0x0010
#define BIT5        0x0020
#define BIT6        0x0040          //bit位
#define BIT7        0x0080

#define NODE_MAX_NUM          600         //最大节点数
#define OUT_DEGREE_NUM          8          //最大出度
#define NODE_INFO_LINE          4          //一行数据包含的信息

#define NODE_MAX_ROUTS         50          //最多50个顶点

#define INVALID                -1          //负值

#define NODE_NUMBER          BITX         //0 节点数量
#define NODE_TYPE            BIT0         //1 节点类型

#define SUB_NODE_INDEX      BITX          //0
#define SUB_NODE_NEXT       BIT0          //1


#define NODE_TYPE_COMM      BITX          //通用节点
#define NODE_TYPE_START     BIT0          //起始节点
#define NODE_TYPE_END       BIT1          //终止节点
#define NODE_TYPE_MID       BIT2          //中间节点

//蚂蚁信息
#define INFO_START_NODE           0x0001 //起始节点
#define INFO_END_NODE             0x0002 //结束节点
#define INFO_ANT_NUMBER           0x0004 //蚂蚁数量
#define INFO_MUST_NUMBER          0x0008 //必经节点数量
#define INFO_NODE_NUMBER          0x0010 //顶点总数量
#define INFO_ALPHA_VAL            0x0020 //阿尔法值
#define INFO_BETA_VAL             0x0040 //贝塔值
#define INFO_BEST_DIST            0x0080 //最佳距离

#define INFO_MUST_DATA            0x0001 //必经数据
#define INFO_BEST_ROUT            0x0002 //最佳路由
//#define PHEROMONE_RATE            0.1
#define PHEROMONE_INIT            0.5     //信息素
#define MAX_WIGHT_VALUE           21     //最大权重
#define POSITIVE_CONTS            0.75
#define EVAPORATION_RATE          0.5
#define MIN_PHEROMONE             0.2
#define MAX_PHEROMONE             0.8
#define ANT_NUMBER                500      //蚂蚁数量
#define ANT_ALPHA                 1
#define ANT_BETA                  1

#define ANT_UNGO                  0
#define ANT_GOST                  1
#define ANT_TABU                  2

#define MAX_DIST_VAL              32767
#define MAX_ITERATIONS            100
//#define INIT_PROB                 0.83   //初始概率
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//存储有向图顶点节点映射
typedef struct Out_Degree{
  short data[2];                        //索引 下一节点
}*O_Node;

typedef struct Topo_Node{
    short data[2];                      //个数 类型 
    O_Node  node;                       //索引指针
}*T_Node;

typedef struct Ant_Node{

	short    dist;                    //距离
	double   fitn;                    //适应度
	short    posi;                    //记录当前位置
	short*   rout;                    //小蚂蚁走过的路
	short    node[NODE_MAX_NUM];      //保存访问信息
	short*   tabu;                    //禁忌表

}*A_Node;                             //一只萌萌哒小蚂蚁

struct Mat_Node{

	short*  data[NODE_MAX_NUM];       //600x600的矩阵

}; //用于计算信息素和保存连通图
struct Pher_Node{
	double* data[NODE_MAX_NUM];       //600x600的矩阵
};
struct Info_Node{

	short   startNode;                     //起始节点
	short   endNode;                       //结束节点

	short   mustNumber;                    //必经点数量
	short   antNumber;                     //蚂蚁数量

	short   nodeNumber;                    //节点总数量
    short   alpha;                         //alpha参数

	short   beta;                          //beta参数
	short   maxWeight;                     //最大权重

	short   bestDist;                      //最好节点
	float   bestFitn;                      //最好适应度
	short   mustData[NODE_MAX_ROUTS];      //必经节点初始化,最多50个
	short   bestLength;
	short*  bestRout;                      //最好的路由,申请节点总量大小的节点
	short*  leftRout;
	short*  rightRout;
	Pher_Node pheromone;                    //信息素矩阵
	Mat_Node distances;                    //权重矩阵

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//信息节点初始化
void Info_MustData_Init(Info_Node& info);
//初始化其他需要的信息 antNumber,maxWeight bestDist alpha beta
void Info_Some_Init(Info_Node& info,short antNumber,short maxWeight,short bestDist,short alpha,short beta);

//蚂蚁节点初始化 一个参数，返回A_Node行数据，没有分配存储空间时调用该函数
A_Node Ant_Init(Info_Node& info);
//蚂蚁节点初始化 两个参数，返回void,复位节点数据使用该函数
void Ant_Init(A_Node Ants,Info_Node& info);
//释放蚂蚁节点分配空间
void FreeAnts(A_Node Ants,Info_Node& info);
void FreeInfo(Info_Node& info);

//检测节点是否访问过,访问过返回true，否则返回false  T_Node Head,int start,int offset,A_Node Ants,int index
bool Check_Visit(T_Node Head,int start,int offset,A_Node Ants,int index);
bool Check_Visit(A_Node Ants,int index,int position);
//检查是否为必经节点
bool Check_Visit(T_Node Head,int offset);


//链表节点 蚂蚁节点 蚂蚁数量 开始节点 数量
void Build_Solution(T_Node Head,A_Node Ants,Info_Node& info,double init_prob);
void Check_Best_Distance(A_Node Ants,Info_Node& info);
void Init_Pheromone(T_Node Head,Info_Node& info);
void Update_Best_Route(T_Node Head,A_Node Ants,Info_Node& info);
void Calculate_Fitness(A_Node Ants,Info_Node& info);

//计算挥发洗漱
void Pheromone_Evaporates(T_Node Head,Info_Node& info);

void Update_Pheromone(T_Node Head,A_Node Ants,Info_Node& info);
int Get_Random(int from, int to);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool isNodeEmpty(const T_Node Head,const short offset);

//返回节点内容 NODE_NUMBER(个数) NODE_TYPE（类型）
short getNode(const T_Node Head,const short i,const short flag);
//设置节点内容 NODE_NUMBER(个数) NODE_TYPE（类型）
void setNode(T_Node Head,short i,const short data,const short flag);
//返回当前节点的第j个点 SUB_NODE_INDEX SUB_NODE_NEXT
short getSubNode(const T_Node Head,short i,short j,const short flag);
//设置当前节点的第j个点 SUB_NODE_INDEX SUB_NODE_NEXT
void setSubNode(const T_Node Head,const short i,const short j,const short data,const short flag);
//计数器加 链表 偏移值 步进值
void addSubNodeNumber(const T_Node Head,const short i,short steps);
//计数器减 链表 偏移值 步进值
void desSubNodeNumber(const T_Node Head,const short i,short steps);
//释放头节点
void FreeHead(Topo_Node* Head);

//初始化头节点
Topo_Node* Node_Init(const int num);

//初始化子节点（尾部节点）
O_Node SubNode_Init(const short num);
//存储待查找点
void Storage_Demand(T_Node Head,               //链表节点
					Info_Node& info,
					map<short,char>& demMap,   //中间节点
					const string str,const string dela,const string delb);
//读取有向图进行存储,一行数据的存储
void Storage_Topo(T_Node Head,Info_Node& info,const string str,const string del);
//搜索整体初始化
void Search_Init(T_Node Head,       //邻接表初始化
				 Info_Node& info,   //信息表初始化
				 char *topo[5000],  //边的数据
				 int edge_num,      //边的数量
				 char *demand       //必经点数据
				 );
string number_to_String(double n);
void Print_Route(A_Node Ants,Info_Node& info);
void Print_Route(A_Node Ants,Info_Node& info,ofstream& outfile,T_Node Head);
void DisplayTopo(T_Node Head,Info_Node& info);
void DisplayDemand(Info_Node& info);
void DisplayBestWay(Info_Node& info);
void ComBestWay(Info_Node& info);
#endif
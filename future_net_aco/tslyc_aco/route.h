#ifndef __ROUTE_H__
#define __ROUTE_H__
/**
    *日期 : 2016/04/03    我猜后面这个日期就不更新了
    *作者 : tslyc                 这是集体的智慧,反正我不会告诉你其他人是谁
    *说明 : 耗尽了毕生功力还是没能有个好的结局的话,先烧三炷香祭奠一下清明节逝去的青春吧.
    *            接下来就是正文了,我看起来好闲,还有心情写这么多注释,算了,刚芭蕾古达腮!!!
    */
////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <string>
#include <string.h>
#include <stdlib.h> 
#include <math.h>
 #include <time.h>
using namespace std;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//宏定义部分
//bit位定义
#define BITX        0x0000                            //0
#define BIT0        0x0001                            //1
#define BIT1        0x0002                            //2
#define BIT2        0x0004                            //4
#define BIT3        0x0008                            //8
#define BIT4        0x0010                            //16
#define BIT5        0x0020                            //32
#define BIT6        0x0040                            //64
#define BIT7        0x0080                            //128

#define INVALID                                         -1           //负值
#define RANDOM_RANGE                        100       //随机数种子范围默认值
//邻接表与矩阵相关定义
#define NODE_MAX_NUM                      600         //最大节点数
#define OUT_DEGREE_NUM                       8          //最大出度
#define NODE_INFO_LINE                         4          //顶点图,一行数据包含的信息
#define NODE_MAX_ROUTS                     50         //最多50个顶点

//用于getNode()函数获取节点类型和数量
#define NODE_NUMBER                      BITX           //0   节点数量
#define NODE_TYPE                              BIT0            //1   节点类型

//用于getSubNode()函数获取边的尾节点
#define SUB_NODE_INDEX                  BITX          //0     节点索引,目前没什么用
#define SUB_NODE_NEXT                    BIT0           //1     节点出度方向

//用于设置节点的类型
#define NODE_TYPE_COMM                 BITX          //通用节点
#define NODE_TYPE_START                   BIT0          //起始节点
#define NODE_TYPE_END                       BIT1          //终止节点
#define NODE_TYPE_MID                       BIT2         //中间节点

//蚂蚁信
#define PHEROMONE_INIT                    0.5                   //信息素初始浓度
#define MAX_WIGHT_VALUE                   21                   //最大权重,本题为21
#define POSITIVE_CONTS                      0.75                  
#define EVAPORATION_RATE                 0.5                 //挥发率
#define MIN_PHEROMONE                     0.1
#define MAX_PHEROMONE                     0.9

#define ANT_NUMBER                             300                //蚂蚁数量
#define ANT_ALPHA                                     1                 //alpha系数
#define ANT_BETA                                         2                 //beta系数

#define ANT_UNGO                                       0                //蚂蚁未经过
#define ANT_GOST                                         1                //蚂蚁经过
#define ANT_TABU                                         2                //蚂蚁在进忌表

#define MAX_DIST_VAL                               32767            //最大距离值
#define MAX_ITERATIONS                         15                  //最大迭代次数

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//存储有向图顶点节点映射
//出度节点
typedef struct Out_Degree{
  short data;                                                                             //目前仅保存下一节点
}*O_Node;
//顶点图索引
typedef struct Topo_Node{
    short data[2];                                                                      //出度总数  顶点类型 
    O_Node  node;                                                                   // 索引指针
}*T_Node;

//呆萌的小蚂蚁节点
typedef struct Ant_Node{
	short     dist;                                                     //距离和
	double   fitn;                                                    //适应度
	short    posi;                                                     //记录当前节点位置
	short*   rout;                                                   //小蚂蚁走过的路
	short    node[NODE_MAX_NUM];             //保存访问信息,记录是否在禁忌表或者已经过
	short*   tabu;                                                  //禁忌表,保存禁忌点信息
	short    counter;                                             //计数器,计数 必经点个数

}*A_Node;                                                                        //一只萌萌哒小蚂蚁

//600X600的矩阵就是这么嚣张
struct Mat_Node{
	short*  data[NODE_MAX_NUM];             //600x600的矩阵
}; //用于保存连通图和索引节点
//信息素矩阵
struct Pher_Node{
	double* data[NODE_MAX_NUM];          //600x600的矩阵
};
//一个杂乱的信息表,啥啥都在这儿了,算我机智吧
struct Info_Node{

	short   startNode;                           //起始节点
	short   endNode;                            //结束节点

	short   mustNumber;                    //必经点数量
	short   antNumber;                       //蚂蚁数量

	short   nodeNumber;                    //节点总数量
	short   maxWeight;                        //最大权重

                short   alpha;                                    //alpha参数
	short   beta;                                      //beta参数
	
	short    bestDist;                              //最好节点
	short   bestLength;                        //最佳路径长度
                short*  bestRout;                           //最好的路由,申请节点总量大小的节点

	double      bestFitn;                             //最好适应度
	double      initProb;                             //初始化概率

	short    mustData[NODE_MAX_ROUTS];      //必经节点初始化,最多50个

	Mat_Node  distances;                                         //权重矩阵
	Mat_Node  edgeindex;                                       //边索引矩阵

	Pher_Node pheromone;                                   //信息素矩阵
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//信息节点初始化
void Info_MustData_Init(Info_Node& info);
//初始化其他需要的信息 antNumber,maxWeight bestDist alpha beta
void Info_Some_Init(Info_Node& info,short antNumber,short maxWeight,short bestDist,short alpha,short beta);
//释放信息表中的堆内存
void FreeInfo(Info_Node& info);

//蚂蚁节点初始化 一个参数，返回A_Node行数据，没有分配存储空间时调用该函数
A_Node Ant_Init(Info_Node& info);
//蚂蚁节点初始化 两个参数，返回void,复位节点数据使用该函数
void Ant_Init(A_Node Ants,Info_Node& info);
//释放蚂蚁节点分配空间
void FreeAnts(A_Node Ants,Info_Node& info);


//检测节点是否访问过,访问过返回true，否则返回false  T_Node Head,int start,int offset,A_Node Ants,int index
bool Check_Visit(T_Node Head,int start,int offset,A_Node Ants,int index);
bool Check_Visit(A_Node Ants,int index,int position);
//检查是否为必经节点
bool Check_Visit(T_Node Head,int offset);
//检查禁忌表
bool Check_Tabu(A_Node Ants,int index,int position);

//链表节点 蚂蚁节点 蚂蚁数量 开始节点 数量
void Build_Solution(T_Node Head,A_Node Ants,Info_Node& info);
//检查最佳路由
void Check_Best_Distance(A_Node Ants,Info_Node& info);
//更新适应度
void Calculate_Fitness(A_Node Ants,Info_Node& info);

//计算挥发系数
void Pheromone_Evaporates(T_Node Head,Info_Node& info);
//更新信息素矩阵
void Update_Pheromone(T_Node Head,A_Node Ants,Info_Node& info);
//设置随机数初始化种子 ,默认 seed =  0 和 range = 0 
void Random_Init(const unsigned int seed = 0,const int range = 0);
//获取随机数
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

//初始化头节点
Topo_Node* Node_Init(const int num);
//初始化子节点（尾部节点）
O_Node SubNode_Init(const short num);
//释放头节点
void FreeHead(Topo_Node* Head);

//存储待查找点
void Storage_Demand(T_Node Head,               //链表节点
		             Info_Node& info,        //基本信息表
		             const string str,const string dela,const string delb); //你猜这几个参数的作用,看看代码就好了
//读取有向图进行存储,一行数据的存储
void Storage_Topo(T_Node Head,Info_Node& info,const string str,const string del);
//搜索整体初始化,先初始化Info_Node必须点,再初始化必经点,在初始化图节点,最后初始化Info_Node剩余内容,好复杂的样子 > _ >
void Search_Init(T_Node Head,                   //邻接表初始化
		Info_Node& info,             //信息表初始化
		 char *topo[5000],           //边的数据
		int edge_num,                  //边的数量
		char *demand                  //必经点数据
		);
//打印输出最好的路径,这个不是真的打印啦,是输出到华大为给的接口而已
void  PrintBestRouter(Info_Node& info);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//原始代码行
void search_route(char *graph[5000], int edge_num, char *condition);

#endif

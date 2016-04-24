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
//�궨��
//bitλ����
#define BITX        0x0000
#define BIT0        0x0001          //bitλ
#define BIT1        0x0002
#define BIT2        0x0004
#define BIT3        0x0008          //bitλ
#define BIT4        0x0010
#define BIT5        0x0020
#define BIT6        0x0040          //bitλ
#define BIT7        0x0080

#define NODE_MAX_NUM          600         //���ڵ���
#define OUT_DEGREE_NUM          8          //������
#define NODE_INFO_LINE          4          //һ�����ݰ�������Ϣ

#define NODE_MAX_ROUTS         50          //���50������

#define INVALID                -1          //��ֵ

#define NODE_NUMBER          BITX         //0 �ڵ�����
#define NODE_TYPE            BIT0         //1 �ڵ�����

#define SUB_NODE_INDEX      BITX          //0
#define SUB_NODE_NEXT       BIT0          //1


#define NODE_TYPE_COMM      BITX          //ͨ�ýڵ�
#define NODE_TYPE_START     BIT0          //��ʼ�ڵ�
#define NODE_TYPE_END       BIT1          //��ֹ�ڵ�
#define NODE_TYPE_MID       BIT2          //�м�ڵ�

//������Ϣ
#define INFO_START_NODE           0x0001 //��ʼ�ڵ�
#define INFO_END_NODE             0x0002 //�����ڵ�
#define INFO_ANT_NUMBER           0x0004 //��������
#define INFO_MUST_NUMBER          0x0008 //�ؾ��ڵ�����
#define INFO_NODE_NUMBER          0x0010 //����������
#define INFO_ALPHA_VAL            0x0020 //������ֵ
#define INFO_BETA_VAL             0x0040 //����ֵ
#define INFO_BEST_DIST            0x0080 //��Ѿ���

#define INFO_MUST_DATA            0x0001 //�ؾ�����
#define INFO_BEST_ROUT            0x0002 //���·��
//#define PHEROMONE_RATE            0.1
#define PHEROMONE_INIT            0.5     //��Ϣ��
#define MAX_WIGHT_VALUE           21     //���Ȩ��
#define POSITIVE_CONTS            0.75
#define EVAPORATION_RATE          0.5
#define MIN_PHEROMONE             0.2
#define MAX_PHEROMONE             0.8
#define ANT_NUMBER                500      //��������
#define ANT_ALPHA                 1
#define ANT_BETA                  1

#define ANT_UNGO                  0
#define ANT_GOST                  1
#define ANT_TABU                  2

#define MAX_DIST_VAL              32767
#define MAX_ITERATIONS            100
//#define INIT_PROB                 0.83   //��ʼ����
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//�洢����ͼ����ڵ�ӳ��
typedef struct Out_Degree{
  short data[2];                        //���� ��һ�ڵ�
}*O_Node;

typedef struct Topo_Node{
    short data[2];                      //���� ���� 
    O_Node  node;                       //����ָ��
}*T_Node;

typedef struct Ant_Node{

	short    dist;                    //����
	double   fitn;                    //��Ӧ��
	short    posi;                    //��¼��ǰλ��
	short*   rout;                    //С�����߹���·
	short    node[NODE_MAX_NUM];      //���������Ϣ
	short*   tabu;                    //���ɱ�

}*A_Node;                             //һֻ������С����

struct Mat_Node{

	short*  data[NODE_MAX_NUM];       //600x600�ľ���

}; //���ڼ�����Ϣ�غͱ�����ͨͼ
struct Pher_Node{
	double* data[NODE_MAX_NUM];       //600x600�ľ���
};
struct Info_Node{

	short   startNode;                     //��ʼ�ڵ�
	short   endNode;                       //�����ڵ�

	short   mustNumber;                    //�ؾ�������
	short   antNumber;                     //��������

	short   nodeNumber;                    //�ڵ�������
    short   alpha;                         //alpha����

	short   beta;                          //beta����
	short   maxWeight;                     //���Ȩ��

	short   bestDist;                      //��ýڵ�
	float   bestFitn;                      //�����Ӧ��
	short   mustData[NODE_MAX_ROUTS];      //�ؾ��ڵ��ʼ��,���50��
	short   bestLength;
	short*  bestRout;                      //��õ�·��,����ڵ�������С�Ľڵ�
	short*  leftRout;
	short*  rightRout;
	Pher_Node pheromone;                    //��Ϣ�ؾ���
	Mat_Node distances;                    //Ȩ�ؾ���

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//��Ϣ�ڵ��ʼ��
void Info_MustData_Init(Info_Node& info);
//��ʼ��������Ҫ����Ϣ antNumber,maxWeight bestDist alpha beta
void Info_Some_Init(Info_Node& info,short antNumber,short maxWeight,short bestDist,short alpha,short beta);

//���Ͻڵ��ʼ�� һ������������A_Node�����ݣ�û�з���洢�ռ�ʱ���øú���
A_Node Ant_Init(Info_Node& info);
//���Ͻڵ��ʼ�� ��������������void,��λ�ڵ�����ʹ�øú���
void Ant_Init(A_Node Ants,Info_Node& info);
//�ͷ����Ͻڵ����ռ�
void FreeAnts(A_Node Ants,Info_Node& info);
void FreeInfo(Info_Node& info);

//���ڵ��Ƿ���ʹ�,���ʹ�����true�����򷵻�false  T_Node Head,int start,int offset,A_Node Ants,int index
bool Check_Visit(T_Node Head,int start,int offset,A_Node Ants,int index);
bool Check_Visit(A_Node Ants,int index,int position);
//����Ƿ�Ϊ�ؾ��ڵ�
bool Check_Visit(T_Node Head,int offset);


//����ڵ� ���Ͻڵ� �������� ��ʼ�ڵ� ����
void Build_Solution(T_Node Head,A_Node Ants,Info_Node& info,double init_prob);
void Check_Best_Distance(A_Node Ants,Info_Node& info);
void Init_Pheromone(T_Node Head,Info_Node& info);
void Update_Best_Route(T_Node Head,A_Node Ants,Info_Node& info);
void Calculate_Fitness(A_Node Ants,Info_Node& info);

//����ӷ�ϴ��
void Pheromone_Evaporates(T_Node Head,Info_Node& info);

void Update_Pheromone(T_Node Head,A_Node Ants,Info_Node& info);
int Get_Random(int from, int to);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool isNodeEmpty(const T_Node Head,const short offset);

//���ؽڵ����� NODE_NUMBER(����) NODE_TYPE�����ͣ�
short getNode(const T_Node Head,const short i,const short flag);
//���ýڵ����� NODE_NUMBER(����) NODE_TYPE�����ͣ�
void setNode(T_Node Head,short i,const short data,const short flag);
//���ص�ǰ�ڵ�ĵ�j���� SUB_NODE_INDEX SUB_NODE_NEXT
short getSubNode(const T_Node Head,short i,short j,const short flag);
//���õ�ǰ�ڵ�ĵ�j���� SUB_NODE_INDEX SUB_NODE_NEXT
void setSubNode(const T_Node Head,const short i,const short j,const short data,const short flag);
//�������� ���� ƫ��ֵ ����ֵ
void addSubNodeNumber(const T_Node Head,const short i,short steps);
//�������� ���� ƫ��ֵ ����ֵ
void desSubNodeNumber(const T_Node Head,const short i,short steps);
//�ͷ�ͷ�ڵ�
void FreeHead(Topo_Node* Head);

//��ʼ��ͷ�ڵ�
Topo_Node* Node_Init(const int num);

//��ʼ���ӽڵ㣨β���ڵ㣩
O_Node SubNode_Init(const short num);
//�洢�����ҵ�
void Storage_Demand(T_Node Head,               //����ڵ�
					Info_Node& info,
					map<short,char>& demMap,   //�м�ڵ�
					const string str,const string dela,const string delb);
//��ȡ����ͼ���д洢,һ�����ݵĴ洢
void Storage_Topo(T_Node Head,Info_Node& info,const string str,const string del);
//���������ʼ��
void Search_Init(T_Node Head,       //�ڽӱ��ʼ��
				 Info_Node& info,   //��Ϣ���ʼ��
				 char *topo[5000],  //�ߵ�����
				 int edge_num,      //�ߵ�����
				 char *demand       //�ؾ�������
				 );
string number_to_String(double n);
void Print_Route(A_Node Ants,Info_Node& info);
void Print_Route(A_Node Ants,Info_Node& info,ofstream& outfile,T_Node Head);
void DisplayTopo(T_Node Head,Info_Node& info);
void DisplayDemand(Info_Node& info);
void DisplayBestWay(Info_Node& info);
void ComBestWay(Info_Node& info);
#endif
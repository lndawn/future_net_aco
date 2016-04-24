#include "myFunc.h"


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//初始化一次
//初始化信息素矩阵和距离矩阵,设计一个600X600的距离矩阵和信息素矩阵
void Info_MustData_Init(Info_Node& info){
	
	for(int i = 0;i<NODE_MAX_NUM;i++){   
		info.pheromone.data[i] = new double[NODE_MAX_NUM];
		info.distances.data[i] = new short[NODE_MAX_NUM];
		for(int j = 0;j<NODE_MAX_NUM;j++){
			info.pheromone.data[i][j] = INVALID;  //信息素矩阵
			info.distances.data[i][j] = INVALID;  //距离值矩阵
		}
	}
}
/*
void Info_MustData_Init(Info_Node& info){
	for(int i = 0;i<info.nodeNumber;i++){    
		for(int j = 0;j<info.nodeNumber;j++){
			info.pheromone.data[i][j] = INVALID;
			info.distances.data[i][j] = INVALID;
		}
	}
}
*/

//初始化其他需要的信息 antNumber,maxWeight bestDist alpha beta 
void Info_Some_Init(Info_Node& info,short antNumber,short maxWeight,short bestDist,short alpha,short beta){
	info.antNumber = antNumber;    //蚂蚁数量
	info.maxWeight = maxWeight;    //最大权重
	info.alpha = alpha;            //alpha
	info.beta = beta;              //beta

	info.bestDist = bestDist;           //初始化为一个最大值
	info.bestFitn = MAX_DIST_VAL;       //最好适应度
	info.bestLength = 0;                //最好矩阵长度
	if(info.nodeNumber > 0){            //如果节点数大于零
		info.bestRout = new short[info.nodeNumber];
		for(int i = 0;i<info.nodeNumber;i++){//初始化最好路由节点 
			info.bestRout[i] = INVALID;
		}
	}else{
		//输出节点数不合法信息
		cout<<"Error: nodeNumber is  Invalid!!"<<endl;
	}
}

//蚂蚁节点初始化 一个参数，返回A_Node行数据，没有分配存储空间时调用该函数
A_Node Ant_Init(Info_Node& info){
	A_Node Ants = NULL;                              //声明一个蚂蚁节点指针

	if(info.antNumber > 0){

		Ants = new Ant_Node[info.antNumber];         //申请大小为number的节点
		for(int i = 0;i < info.antNumber; i++){
			(Ants+i)->dist = 0;                      //初始化距离零
			(Ants+i)->fitn = 0;                      //初始化适应度
			(Ants+i)->posi = 0;                      //初始化位置点

			(Ants+i)->rout = new short[info.nodeNumber];  //初始化路由节点总大小
			(Ants+i)->tabu = new short[info.nodeNumber];  //初始化禁忌表

			(Ants+i)->rout[0] = info.startNode;           //获取起始点值
			(Ants+i)->tabu[0] = INVALID;                  //初始化禁忌表赋值

			for (int j = 1; j < info.nodeNumber; j++){    //初始化蚂蚁节点值
				(Ants+i)->rout[j] = INVALID;
				(Ants+i)->tabu[j] = INVALID;
			}
			for(int j = 0; j < NODE_MAX_NUM; j++){
				(Ants+i)->node[i] = ANT_UNGO;             //初始化所有节点为未经过
			}
			(Ants+i)->node[info.startNode] = ANT_GOST;  
		}
	}
	return Ants;
}

//蚂蚁节点初始化 两个参数，返回void,复位节点数据使用该函数
void Ant_Init(A_Node Ants,Info_Node& info){
	if(Ants != NULL){
		for(int i = 0;i < info.antNumber; i++){
			(Ants+i)->dist = 0;   //初始化距离零
			(Ants+i)->fitn = 0;   //初始化适应度
			(Ants+i)->posi = 0;   //初始化位置点

			(Ants+i)->rout[0] = info.startNode;         //获取起始点值
			(Ants+i)->node[0] = ANT_UNGO;               //初始化所有节点为未经过
			(Ants+i)->tabu[0] = INVALID;
			for (int j = 1; j < info.nodeNumber; j++){  //初始化蚂蚁节点值
				(Ants+i)->rout[j] = INVALID;
				(Ants+i)->tabu[j] = INVALID;            
				(Ants+i)->node[j] = ANT_UNGO;           //初始化所有节点为未经过
			}
		}
	}
}
//释放蚂蚁节点分配空间
void FreeAnts(A_Node Ants,Info_Node& info){
	if(Ants != NULL){
		for(int i = 0;i<info.antNumber;i++){
			if((Ants+i)->rout != NULL){
				delete[] ((Ants+i)->rout);
			}
			if((Ants+i)->tabu != NULL){
				delete[] ((Ants+i)->tabu);
			}
		}
		delete[] Ants;                     //释放头节点
	}
}
//释放信息素和距离矩阵
void FreeInfo(Info_Node& info){
	for(int i = 0;i<NODE_MAX_NUM;i++){   
		delete[] info.pheromone.data[i];
		delete[]info.distances.data[i];
	}
}
//检测节点是否访问过,访问过返回true，否则返回no，未知idata
bool Check_Visit(T_Node Head,int start,int offset,A_Node Ants,int index){

	int idata = getSubNode(Head,start,offset,SUB_NODE_NEXT); //获取尾节点

	if((Ants+index)->node[idata] == ANT_GOST){               //检测尾节点是否为经过
		return true;
	}
	return false;
}
//已知idata
bool Check_Visit(A_Node Ants,int index,int position){
	if((Ants+index)->node[position] == ANT_GOST){           //检测尾节点是否为经过
		return true;
	}
	return false;

}
//检查是否为必经节点
bool Check_Visit(T_Node Head,int offset)
{
	if(getNode(Head,offset,NODE_TYPE) == NODE_TYPE_MID){    //节点为必经节点
		return true;
	}
	return false;
}
//检测禁忌表
bool Check_Tabu(A_Node Ants,int index,int position){
	if((Ants+index)->node[position] == ANT_TABU){           //检测尾节点是否为经过
		return true;
	}
	return false;
}
//链表节点 蚂蚁节点 蚂蚁数量 开始节点 数量
void Build_Solution(T_Node Head,A_Node Ants,Info_Node& info,double init_prob) {

	int nodeNumber = info.nodeNumber;                      //获取节点总数
	double* transfer;                                      //转移概率分子值
	double linkRate = 0;                                   //连接度总和
	short count=0;                                         //计数器
	short tabuCount = 0;                                  //禁忌表计数器
	short tabuTemp=-1;
	if(nodeNumber>0){
	 transfer = new double[nodeNumber];                    //申请空间
	}
	Ant_Init(Ants,info);                                   //初始化蚂蚁节点，清空

	for (int i = 0; i < info.antNumber; i++) {   //根据蚂蚁只数完成一次循环
		count=0;
		short cycleCount = 0;
		tabuCount=0;
		while ((Ants+i)->rout[(Ants+i)->posi] != info.endNode){//若当前节点不是终止节点且走过的点数小于总点数
			cycleCount++;
			int position = (Ants+i)->rout[(Ants+i)->posi]; //获取当前位置节点值，相当于A->B中的A

			memset(transfer,0,nodeNumber);                 //清零概率分子值

			linkRate = 0;                                  //初始化概率分母和
			for (int j = 0; j < getNode(Head,position,NODE_NUMBER); j++){//获取当前节点的出度值
				short idata = getSubNode(Head,position,j,SUB_NODE_NEXT);// 获取出度j尾节点的值
				if (Check_Visit(Ants,i,idata) == false){
					if(info.distances.data[position][idata]>0 && Check_Visit(Ants,i,idata) == false && Check_Tabu(Ants,i,idata)==false && idata!=tabuTemp){
						if(idata == info.endNode && count != info.mustNumber){
							continue;}
						else{
							if(Check_Visit(Head,idata) && (Ants+i)->posi>1){
									info.pheromone.data[position][idata] *= 1.2;
									info.pheromone.data[(Ants+i)->rout[(Ants+i)->posi-1]][position] *= 1.2;
								}
						double current = pow(info.pheromone.data[position][idata], info.alpha) * pow(double(info.maxWeight - info.distances.data[position][idata]),info.beta);;
							linkRate += current; //连接度加和     
							transfer[idata] = current;//当前节点量值
						}
					}
				}
			}

			for (short j = 0; j < getNode(Head,position,NODE_NUMBER); j++){
				short idata = getSubNode(Head,position,j,SUB_NODE_NEXT);      //获取尾节点值
				if(info.distances.data[position][idata]>0 && Check_Visit(Ants,i,idata) == false && Check_Tabu(Ants,i,idata)==false){ 
					if(idata == info.endNode && count != info.mustNumber){
							continue;}//查看当前节点是否访问过
					else
					transfer[idata] = transfer[idata] / linkRate;
				}
			}

			double roulette = (double) Get_Random(0, 100) / 100.0;
			double minor = 0;
			double major = 0;
			short selectidata = -1;
			if(roulette<=init_prob)
			{ 
				double temp=0;
				for (int j = 0; j <getNode(Head,position,NODE_NUMBER); j++){
					short idata = getSubNode(Head,position,j,SUB_NODE_NEXT);
					if (Check_Visit(Ants,i,idata) == false && info.distances.data[position][idata]>0 && Check_Tabu(Ants,i,idata)==false)
					{
						if(idata == info.endNode && count != info.mustNumber){
							continue;}//查看当前节点是否访问过
						else
						{
							if(transfer[idata]>temp)
							{
								temp=transfer[idata];
								selectidata=idata;
							}
						}
					}
				}
			}
			else{
				for (int j = 0; j <getNode(Head,position,NODE_NUMBER); j++)
				{
					short idata = getSubNode(Head,position,j,SUB_NODE_NEXT);
					if(info.distances.data[position][idata]>0 && Check_Visit(Ants,i,idata) == false && Check_Tabu(Ants,i,idata)==false) 
					{
						if(idata == info.endNode && count != info.mustNumber){
							continue;}//查看当前节点是否访问
						else
						{
							major += transfer[idata];
							if (roulette >= minor && roulette <= major) 
							{/*
								(Ants+i)->posi += 1;                    
								(Ants+i)->rout[(Ants+i)->posi] = idata;  

								(Ants+i)->node[idata] = ANT_GOST;      
								(Ants+i)->dist += info.distances.data[position][idata];
								info.pheromone.data[position][idata] *= 0.8;

								if(Check_Visit(Head,idata))
								{
									count++;
									info.pheromone.data[position][idata] *= 1.5;
								}*/
								selectidata=idata;
								break;
							} 
						}
					}
						else
						{
							minor = major;
						}
					}
			}
			if(selectidata==-1)
			{
				if(Check_Visit(Head,(Ants+i)->rout[(Ants+i)->posi]))//查看是否为必经点
				{
					count--;
					(Ants+i)->rout[(Ants+i)->posi]=ANT_UNGO;
					(Ants+i)->tabu[tabuCount] = INVALID;
				}
				else
				{
					(Ants+i)->tabu[tabuCount] = (Ants+i)->rout[(Ants+i)->posi];
					(Ants+i)->node[position] = ANT_TABU;      //添加到已经过
				}
				(Ants+i)->rout[(Ants+i)->posi] = INVALID;
				(Ants+i)->posi -= 1;
				(Ants+i)->dist -= info.distances.data[(Ants+i)->rout[(Ants+i)->posi]][position];
				
				tabuCount++;
				//tabuTemp=(Ants+i)->rout[(Ants+i)->posi];
				if((tabuCount>nodeNumber*0.1) || (cycleCount>nodeNumber))
				{
					//cout<<cycleCount<<" "<<endl;
					break;
				}
			}
			else
			{
				(Ants+i)->posi += 1;                     //位置标记++
				(Ants+i)->rout[(Ants+i)->posi] =selectidata;  //保存下一节点值

				(Ants+i)->node[selectidata] = ANT_GOST;        //添加到已经过
				(Ants+i)->dist += info.distances.data[position][selectidata];
				info.pheromone.data[position][selectidata] *= 0.9;
				
				if(Check_Visit(Head,selectidata)){
					count++;
				}
			}
			if(count == info.mustNumber && (Ants+i)->rout[(Ants+i)->posi]==info.endNode)
					break;
		}
	}
	delete[] transfer;
}
void Check_Best_Distance(A_Node Ants,Info_Node& info) {
	int antId = INVALID;
	for (int i = 0; i < info.antNumber; i++) {//根据蚂蚁数量检测

		if ((Ants+i)->dist < info.bestDist && (Ants+i)->rout[(Ants+i)->posi] == info.endNode){
			info.bestDist = (Ants+i)->dist;
			antId = i;//一趟遍历记录最小的代价蚂蚁
		} 
	}
	if(antId != INVALID){
		info.bestLength = (Ants+antId)->posi+1;  //最佳路由长度
		for (int j = 0; j < info.nodeNumber; j++) {
			info.bestRout[j] = (Ants+antId)->rout[j];
		}
	}
}
void Init_Pheromone(T_Node Head,Info_Node& info){
	for(int i=0;i<info.nodeNumber;i++){
		for(int j=0;j<getNode(Head,i,NODE_NUMBER);j++){
			int idata=getSubNode(Head,i,j,SUB_NODE_NEXT);
			if(idata>=0){
				info.pheromone.data[i][idata]=PHEROMONE_INIT;
			}
		}
	}
}
void Update_Best_Route(T_Node Head,A_Node Ants,Info_Node& info){
	short visitCount=0;
	short left_dist=0;
	short right_dist=0;
	short left_endNode;
	short leftRout[600];
	short rightRout[600];
	short orig_dist=0;
	short start_Position=0;
	double* transfer; 
	memset(leftRout,-1,600);
	memset(rightRout,-1,600);
	Ant_Init(Ants,info); 
	//int roulette0 = Get_Random(0, info.bestLength-1);
	int roulette0=Get_Random(0, info.bestLength-1);
	int i=0;
	int k=0;
	for(int j=roulette0;j>0;j--){
		if(Check_Visit(Head,info.bestRout[j])){//将路由表进行分割，左边路由表的position为roulette
			start_Position=j;
			visitCount++;
			break;
		}
	}
	for(i=0;i<start_Position;i++)
	{
		leftRout[i]=info.bestRout[i];
		left_dist+=info.distances.data[info.bestRout[i]][info.bestRout[i+1]];
	}
	leftRout[i]=info.bestRout[i];
	orig_dist=info.distances.data[info.bestRout[i]][info.bestRout[i+1]];//需要比较的距离
	i++;
	while(Check_Visit(Head,info.bestRout[i])==false && i<info.bestLength-1)//统计两个优化位置之间的间距
	{
		orig_dist+=info.distances.data[info.bestRout[i]][info.bestRout[i+1]];
		visitCount++;
		i++;
		
	}
	if(visitCount>3){
		left_endNode=info.bestRout[i];
		while(i<info.bestLength-1){//右边路由表的起点为优化间距的右端点
			rightRout[k]=info.bestRout[i];
			right_dist+=info.distances.data[info.bestRout[i]][info.bestRout[i+1]];
			i++;
			k++;
		}
		(Ants+0)->dist=left_dist;
		while((Ants+0)->posi!=start_Position){
			(Ants+0)->rout[(Ants+i)->posi]=leftRout[(Ants+0)->posi];
			(Ants+0)->node[leftRout[(Ants+0)->posi]] = ANT_GOST;        //添加到已经过
			(Ants+0)->posi += 1;
		}
		for(int j=1;j<=k;j++){
			(Ants+0)->tabu[j]=rightRout[j];
			(Ants+0)->node[rightRout[j]] = ANT_TABU;
		}
		rightRout[k]=info.bestRout[i];
		(Ants+0)->rout[(Ants+0)->posi]=leftRout[(Ants+0)->posi];
		(Ants+0)->node[leftRout[(Ants+0)->posi]] = ANT_GOST; 
		(Ants+1)->rout[(Ants+1)->posi]=leftRout[(Ants+0)->posi];
		int best=MAX_DIST_VAL;
		int next=1;
		int i=2;
		int m,n,temp;
		int nextDist;
		int count;
		ofstream fp;
		while(true){
			int position=(Ants+next)->rout[(Ants+next)->posi];
			count = 0;
			nextDist=(Ants+next)->dist;
			for (int j = 0; j < getNode(Head,position,NODE_NUMBER); j++){//获取当前节点的出度值
				temp=0;
				n=0;
				short idata = getSubNode(Head,position,j,SUB_NODE_NEXT);// 获取出度j尾节点的值
				if (Check_Visit(Ants,0,idata)==false && Check_Tabu(Ants,0,idata)==false && Check_Visit(Ants,next,idata)==false){
					if(info.distances.data[position][idata]>0){
						while(n<i){
							n++;
							if(idata==(Ants+next)->rout[(Ants+next)->posi]){
								if(((Ants+next)->dist+info.distances.data[position][idata])<(Ants+i)->dist){
									for(m=0;m<=(Ants+next)->posi;m++){
											(Ants+n)->rout[m]=(Ants+next)->rout[m];
											(Ants+n)->node[(Ants+n)->rout[m]]=ANT_GOST;
										};
										(Ants+n)->dist=(Ants+next)->dist+info.distances.data[position][idata];
										(Ants+n)->rout[m]=idata;
										(Ants+n)->node[idata]=ANT_GOST;
										(Ants+n)->posi=m;
									}
								temp=1;
								break;
							}
						}
							if(temp==0){
								if(count==0){
									(Ants+next)->posi += 1;
									(Ants+next)->rout[(Ants+next)->posi]=idata;
									(Ants+next)->node[idata]=ANT_GOST;
									(Ants+next)->dist += info.distances.data[position][idata];
								}
								else{
									for(m=0;m<(Ants+next)->posi;m++){
										(Ants+i)->rout[m]=(Ants+next)->rout[m];
										(Ants+i)->node[(Ants+i)->rout[m]]=ANT_GOST;
									}
									(Ants+i)->posi=m;
									(Ants+i)->rout[(Ants+i)->posi]=idata;
									(Ants+i)->node[idata]=ANT_GOST;
									(Ants+i)->dist=nextDist+info.distances.data[position][idata];
									i++;
								}
								count++;
							}
							else
								continue;
					}
				}
			}
			n=1;
			best=MAX_DIST_VAL;
			while(n<i){
				if((Ants+n)->dist<best && Check_Visit(Ants,0,(Ants+n)->rout[(Ants+n)->posi])==false){
					next=n;
					best=(Ants+n)->dist;
				}
				n++;
			}
			(Ants+0)->posi += 1;
			(Ants+0)->rout[(Ants+0)->posi]=(Ants+next)->rout[(Ants+next)->posi];
			(Ants+0)->node[(Ants+0)->rout[(Ants+0)->posi]]=ANT_GOST;
			if((Ants+next)->rout[(Ants+next)->posi]==left_endNode)
			{
				break;
			}
		}
		if((Ants+next)->dist<orig_dist){
			int r=1;
			int i;
			for(i=1;i<=(Ants+next)->posi;i++){
				info.bestRout[start_Position+i]=(Ants+next)->rout[i];
			}
			while(rightRout[r]!=-1){
				info.bestRout[start_Position+i]=rightRout[r];
				i++;
				r++;
			}
			info.bestDist=left_dist+(Ants+next)->dist+right_dist;
			info.bestLength=start_Position+i;
		}
	}			
}

void Calculate_Fitness(A_Node Ants,Info_Node& info) {
	double fitness;
	for (int i = 0; i < info.antNumber; i++) {
		if((Ants+i)->rout[(Ants+i)->posi]==info.endNode){
			fitness = (double) (Ants+i)->dist/ (double) info.bestDist;
			/*if (fitness < info.bestFitn) {
				info.bestFitn = fitness;
			} */
			(Ants+i)->fitn = fitness;
		}
	}
}

//计算挥发洗漱
void Pheromone_Evaporates(T_Node Head,Info_Node& info) {
	for (int i = 0; i < info.nodeNumber; i++) {
		for (int j = 0; j < getNode(Head,i,NODE_NUMBER); j++){
			int idata = getSubNode(Head,i,j,SUB_NODE_NEXT);
			if(idata>=0){
				info.pheromone.data[i][idata] = (1 - EVAPORATION_RATE) * info.pheromone.data[i][idata];
			}
		}
	}
}

void Update_Pheromone(T_Node Head,A_Node Ants,Info_Node& info) {
	double pheromone_to_sum;
	int city;
	int next_city;
	for (int i = 0; i < info.antNumber; i++) {
		if((Ants+i)->rout[(Ants+i)->posi]==info.endNode){
			pheromone_to_sum = POSITIVE_CONTS / (Ants+i)->fitn;

			for (int j = 0; j <((Ants+i)->posi); j++) {
				city = (Ants+i)->rout[j];
				next_city = (Ants+i)->rout[j+1];
				if (info.pheromone.data[city][next_city] != INVALID) {
					info.pheromone.data[city][next_city] += pheromone_to_sum;
					//info.pheromone.data[city][next_city] += EVAPORATION_RATE; 
				}
			}
		}
	}
	/*for (int j = 0; j < info.bestLength-1; j++) {
			info.pheromone.data[info.bestRout[j]][info.bestRout[j+1]] += EVAPORATION_RATE;
	}*/
	for (int i = 0; i < info.nodeNumber; i++) {
		for (int j = 0; j < getNode(Head,i,NODE_NUMBER); j++){
			int idata = getSubNode(Head,i,j,SUB_NODE_NEXT);
			if(idata>=0){
				if(info.pheromone.data[i][idata]>MAX_PHEROMONE)
					info.pheromone.data[i][idata]=MAX_PHEROMONE;
				else if(info.pheromone.data[i][idata]<MIN_PHEROMONE)
					info.pheromone.data[i][idata]=MIN_PHEROMONE;
			}
		}
	}
}
int Get_Random(int from, int to) {
	return (from < to) ? (rand() % to) + from : 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//判断节点是否为空 
bool isNodeEmpty(const T_Node Head,const short offset){
	return ((Head+offset)->node == NULL ? 1 : 0); 
}

//返回节点内容 NODE_NUMBER(个数) NODE_TYPE（类型）
short getNode(const T_Node Head,const short i,const short flag){ 
	return (Head+i)->data[flag];
}
//设置节点内容 NODE_NUMBER(个数) NODE_TYPE（类型）
void setNode(T_Node Head,short i,const short data,const short flag){
	(Head+i)->data[flag] = data;
}


//返回当前节点的第j个点 SUB_NODE_INDEX SUB_NODE_NEXT
short getSubNode(const T_Node Head,short i,short j,const short flag){
	return ((Head+i)->node+j)->data[flag];
}
//设置当前节点的第j个点 SUB_NODE_INDEX SUB_NODE_NEXT
void setSubNode(const T_Node Head,const short i,const short j,const short data,const short flag){
	((Head+i)->node+j)->data[flag] = data;
}
//计数器加 链表 偏移值 步进值
void addSubNodeNumber(const T_Node Head,const short i,short steps){
	(Head+i)->data[0] += steps; //记录当前有多少个节点
}
//计数器减 链表 偏移值 步进值
void desSubNodeNumber(const T_Node Head,const short i,short steps){
	(Head+i)->data[0] -= steps; //记录当前有多少个节点
}

//释放头节点
void FreeHead(Topo_Node* Head){
	if(Head != NULL){
		for(int i = 0;i<NODE_MAX_NUM;i++){
			if(((Head+i)->node)!= NULL){
				free((Head+i)->node);
			}
		}
		free(Head);                     //释放头节点
	}
}

//初始化头节点
Topo_Node* Node_Init(const int num){
  Topo_Node* Head = (Topo_Node*)malloc(sizeof(Topo_Node)*num); //创建头节点索引
  for(int i = 0;i<num;i++){
	  (Head+i)->node    = NULL;                   //初始化尾节点索引;
	  setNode(Head,i,0,NODE_NUMBER);              //初始化节点个数
	  setNode(Head,i,NODE_TYPE_COMM,NODE_TYPE);   //初始化类型为通用类型
  }
  return Head;
}
//初始化子节点（尾部节点）
O_Node SubNode_Init(const short num){
  O_Node Node = (Out_Degree*)malloc(sizeof(Out_Degree)*num);    //一次性分配存储空间并初始化
  for(short i = 0; i<num;i++){
	  Node->data[0] = 0;       //节点索引
	  Node->data[1] = 0;       //下一节点
  }
  return Node;
}
//存储待查找点
void Storage_Demand(T_Node Head,               //链表节点
					Info_Node& info,
					map<short,char>& demMap,   //中间节点
					const string str,const string dela,const string delb){//待分割串 分割符1,分隔符2

	string dema = str+delb;                              //扩展字符串 |
	int size = dema.size();                              //得到数据的大小

	int sdela = dela.size();                             //分割串a大小 ，
	int sdelb = delb.size();                             //分割串b大小 |
	int pos = 0;                                         //起始点

	pos = dema.find(dela,0);                             //逗号分割起始位置
	info.startNode = atoi(dema.substr(0,pos).c_str());          //转换成整数
	setNode(Head,info.startNode,NODE_TYPE_START,NODE_TYPE);     //设置节点类型

	int index = pos+sdela-1;                             
	pos = dema.find(dela,index+1);
	info.endNode = atoi(dema.substr(index+1,pos).c_str());
	setNode(Head,info.endNode,NODE_TYPE_END,NODE_TYPE);      //设置节点类型
	info.mustNumber = 0;                                    //初始化必经点计数器
	for(int i = pos+1;i<size;i++){//分割数据
		pos = dema.find(delb,i);
		if(pos < size){
			index = atoi(dema.substr(i,pos-i).c_str());
			info.mustData[info.mustNumber++] = index;      //添加必经点到Info_Node
			
			setNode(Head,index,NODE_TYPE_MID,NODE_TYPE);   //设置节点类型

			i = pos+sdelb-1;
		}
	}
}
//读取有向图进行存储,一行数据的存储
void Storage_Topo(T_Node Head,Info_Node& info,const string str,const string del){
	int pos = 0;                                          //获取位置节点
	short data[NODE_INFO_LINE];                           //获取数据 4个数据
	int index = 0;                                        //分割数据索引
	string topo = str+del;                                //扩展字符串

	int size = topo.size();                               //得到数据的大小
	int stps = del.size();                                //步进大小

	for(int i = 0;i<size;i++){//分割数据
		pos = topo.find(del,i);
		if(pos < size){
			data[index++] = atoi(topo.substr(i,pos-i).c_str());
			i = pos+stps-1;
		}
	}
	//插入数据

	int counter = (data[1] > data[2]) ? data[1] : data[2]; //比较两个顶点的大小
	if(info.nodeNumber <= counter){
		info.nodeNumber = counter + 1;                       //更新顶点总数
	}

	if(info.distances.data[data[1]][data[2]] == INVALID || info.distances.data[data[1]][data[2]] > data[3]){ //如果权重未保存，或路径重复

		info.pheromone.data[data[1]][data[2]] = PHEROMONE_INIT;                                 //加入信息素初始化
		info.distances.data[data[1]][data[2]] = data[3];                                        //初始化权重

		if(Head != NULL){
		int num = getNode(Head,data[1],NODE_NUMBER);                                            //获取节点数量

		if((Head+data[1])->node == NULL){
		  (Head+data[1])->node = SubNode_Init(OUT_DEGREE_NUM);                                  //一次性分配存储空间
		}
		  setSubNode(Head,data[1],num,data[0],SUB_NODE_INDEX);                                  //索引
		  setSubNode(Head,data[1],num,data[2],SUB_NODE_NEXT);                                   //终点
		  addSubNodeNumber(Head,data[1],1);                                                     //记录当前有多少个节点，值增加1
		}
	}
}

//搜索整体初始化
void Search_Init(T_Node Head,       //邻接表初始化
				 Info_Node& info,   //信息表初始化
				 char *topo[5000],  //边的数据
				 int edge_num,      //边的数量
				 char *demand       //必经点数据
				 )
{
	map<short,char> data;                                           //目前没用
	Info_MustData_Init(info);
	Storage_Demand(Head,info,data,demand,",","|");                  //初始化必经点信息 Info 初始化了start end
	for(int i = 0;i<edge_num;i++){
		Storage_Topo(Head,info,topo[i],",");                        //矩阵信息初始化
	}
	Info_Some_Init(info,ANT_NUMBER, MAX_WIGHT_VALUE,MAX_DIST_VAL,ANT_ALPHA,ANT_BETA); //Info_Node信息初始化

}
string number_to_String(double n) {
	stringstream out;
	out << n;
	return out.str();
}
void Print_Route(A_Node Ants,Info_Node& info){
	string tmp ="";
	
	for(int i = 0; i< info.antNumber;i++){
		tmp ="";
		if((Ants+i)->rout[(Ants+i)->posi] == info.endNode){
			cout<<"蚂蚁["<<i<<"]的路径[Distance = "<<(Ants+i)->dist<<"]:"<<endl;
			for (int j = 0; j <=(Ants+i)->posi; j++) {
				cout<<(Ants+i)->rout[j]<<" ";
			}
			cout<<endl;
		}
	}
}

void Print_Route(A_Node Ants,Info_Node& info,ofstream& outfile,T_Node Head){
	string tmp ="";
	
	for(int i = 0; i< info.antNumber;i++){
		tmp ="";
	//	if((Ants+i)->rout[(Ants+i)->posi] == info.endNode){
			outfile<<"蚂蚁["<<i<<"]的路径[Distance = "<<(Ants+i)->dist<<"]:"<<endl;
			for (int j = 0; j <=(Ants+i)->posi; j++) {
				if(getNode(Head,(Ants+i)->rout[j],NODE_TYPE) == NODE_TYPE_MID){
					outfile<<"{"<<(Ants+i)->rout[j]<<"} ";
				}else{
				   outfile<<(Ants+i)->rout[j]<<" ";
				}
			}
			outfile<<endl;
	//	}
	}
	outfile<<"BestWay[Distance = "<<info.bestDist<<"] : "<<endl;
	for(int i = 0;i<info.bestLength;i++){
		if(getNode(Head,info.bestRout[i],NODE_TYPE) == NODE_TYPE_MID){
			outfile<<"{"<<info.bestRout[i]<<"} ";
		}else{
			outfile<<info.bestRout[i]<<" ";
		}
	}
	outfile<<endl;
}

//显示Topo内容
void DisplayTopo(T_Node Head,Info_Node& info){
	for(int i = 0;i < info.nodeNumber;i++){
		cout<<"{节点【"<<i<<"】 :  \n[TYPE] : "<<getNode(Head,i,NODE_TYPE)<<"\n [NUMBER] : "<<getNode(Head,i,NODE_NUMBER)<<endl;
		for(int j = 0;j<getNode(Head,i,NODE_NUMBER);j++)
		{
			cout<<getSubNode(Head,i,j,SUB_NODE_INDEX)<<":"<<getSubNode(Head,i,j,SUB_NODE_NEXT)<<endl;
		}
		cout<<"}"<<endl;
	}
}
//显示Demand内容
void DisplayDemand(Info_Node& info){
	for(int i = 0;i < info.mustNumber;i++){
		cout<<info.mustData[i]<<" ";
	}
	cout<<endl;
}

//显示最好的路径
void DisplayBestWay(Info_Node& info){
	cout<<"BestWay[Distance = "<<info.bestDist<<"] : "<<endl;
	for(int i = 0;i<info.bestLength;i++){
		cout<<info.bestRout[i]<<" ";
	}
	cout<<endl;
}

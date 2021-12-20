#include <iostream>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const int N = 101;
using namespace::std;
int coordinate4[N][2];
int best_path4[N];


const int MAX_ITER = 1250;
const int ANT_NUM = 20;  //蚂蚁数量
const double Alpha = 1.0, Beta = 2.0, Rho = 0.1, q0 = 0.9, Xi = 0.1;
//Alpha:信息素权重,设置为1表示ACS；Beta:启发式信息权重；Rho:全局信息素挥发参数；Xi：局部信息素挥发因子；q0:伪随机因子q0

double Mutual_Distance4[N][N];	//表示两两城市之间的距离
double Cm4;		//信息素局部更新时候使用的的常量，由贪心算法得到

void Read_opt_path4(char a[])
{
	a = strcat(a, ".opt.tour");
	FILE* fp;
	char str[100];  //暂存读取的字符串
	int i0 = 1, j0 = 0; //i0控制从最优路径文件第几行读取，j0为最优路径赋值下标
	int i = 1, j = 0, m = 0;  //i控制从城市坐标文件第几行读取，j控制只读坐标值，不读城市编号,m为城市坐标赋值下标
	fp = fopen(a, "r");
	while (i0 < 6)
	{
		fgets(str, 255, fp);
		i0++;
	}
	while (!feof(fp) && j0 < N)
	{
		fscanf(fp, "%s\n", str);
		best_path4[j0] = atoi(str);
		j0++;
	}
	fclose(fp);
}

void Read_Coordinate4(char a[])
{
	a = strcat(a, ".tsp");
	FILE* fp;
	char str[100];  //暂存读取的字符串
	int i = 1, j = 0, m = 0;  //i控制从城市坐标文件第几行读取，j控制只读坐标值，不读城市编号,m为城市坐标赋值下标
	fp = fopen(a, "r");
	while (i < 7)
	{
		fgets(str, 255, fp);
		i++;
	}
	while (!feof(fp))
	{
		fscanf(fp, "%s\n", str);
		if (j % 3 == 1) {
			coordinate4[m][0] = atoi(str);
		}
		else if (j % 3 == 2) {
			coordinate4[m][1] = atoi(str);
			m++;
		}
		j++;
	}
	fclose(fp);
}

void Print_bestpath4()  //打印已公布的最优解路径及长度
{
	int i;
	cout << "已公布最优路径为：" << endl;
	for (i = 0; i < N - 1; i++)
	{
		cout << best_path4[i] << " ->";
	}
	cout << best_path4[N - 1] << endl;
	cout << "最优解最小路径长度为：629" << endl << endl << endl;
}


int calculateDistance4(int i, int j)//计算两个城市之间的距离(欧氏距离)
{
	return (int)(sqrt((coordinate4[i][0] - coordinate4[j][0])*(coordinate4[i][0] - coordinate4[j][0]) + (coordinate4[i][1] - coordinate4[j][1])*(coordinate4[i][1] - coordinate4[j][1])) + 0.5);
}



void CalculateDistance4()   //用矩阵表示两两城市之间的距离
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (i != j)
			{
				Mutual_Distance4[i][j] = calculateDistance4(i, j);
				Mutual_Distance4[j][i] = Mutual_Distance4[i][j];
			}
		}
	}
}


double calculateSumOfDistance4(int* tour) //获得经过n个城市的路径长度  
{
	double sum = 0;
	for (int i = 0; i< N; i++)
	{
		int row = *(tour + 2 * i);
		int col = *(tour + 2 * i + 1);
		sum += Mutual_Distance4[row][col];
	}
	return sum;
}


int ChooseNextNode4(int currentNode, int visitedNode[])  //根据当前节点和已访问节点返回距离当前点最近的结点
{
	int nextNode = -1;
	double shortDistance = 0.0;
	for (int i = 0; i < N; i++)
	{
		if (1 == visitedNode[i])	//去掉已走过的节点,从剩下节点中选择距离最近的节点
		{
			if (shortDistance == 0.0)
			{
				shortDistance = Mutual_Distance4[currentNode][i];   //当前节点与节点i距离赋值给shortDistance
				nextNode = i;
			}
			if (shortDistance > Mutual_Distance4[currentNode][i])  //若找到更短距离，更新节点nextNode
			{
				shortDistance = Mutual_Distance4[currentNode][i];   //当前节点与节点i距离赋值给shortDistance
				nextNode = i;
			}
		}
	}
	return nextNode;
}//选择下一个节点，配合下面的函数来计算的长度



class ACS4    //蚁群系统
{
private:
	double pheromone[N][N], visible[N][N];//节点之间的信息素量和节点之间的启发式信息量(两点距离的倒数)
public:
	ACS4()
	{}  //构造函数
	void InitParameter4(double value);	//初始化各路径信息素,每边赋值为value	
	double Transition4(int i, int j);	//计算当前节点到下一节点转移的概率
	void UpdateLocalPathRule4(int i, int j);	//局部更新规则	
	void UpdateGlobalPathRule4(int* bestTour, int globalBestLength);	//全局更新规则
	~ACS4()
	{};  //析构函数
};


void ACS4::InitParameter4(double value) 	//初始化各路径上的信息素强度为value
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			pheromone[i][j] = value;   //信息素初始强度
			pheromone[j][i] = value;
			if (i != j)
			{
				visible[i][j] = 1.0 / Mutual_Distance4[i][j];   //启发式信息强度
				visible[j][i] = visible[i][j];
			}
		}
	}
}


double ACS4::Transition4(int i, int j)//计算当前节点到下一节点转移的概率
{
	if (i != j)
	{
		return (pow(pheromone[i][j], Alpha) * pow(visible[i][j], Beta));  //采用的5.1的概率计算公式
	}
	else
	{
		return 0.0;
	}
}


void ACS4::UpdateLocalPathRule4(int i, int j) //局部更新规则
{
	pheromone[i][j] = (1.0 - Xi) * pheromone[i][j] + Xi * (1.0 / (N * Cm4));
	pheromone[j][i] = pheromone[i][j];
}


void ACS4::UpdateGlobalPathRule4(int* bestTour, int globalBestLength)//全局更新规则
{
	for (int i = 0; i < N; i++)
	{
		int row = *(bestTour + 2 * i);
		int col = *(bestTour + 2 * i + 1);
		pheromone[row][col] = (1.0 - Rho) * pheromone[row][col] + Rho * (1.0 / globalBestLength);//只在当前最优路径更新信息素
		pheromone[col][row] = pheromone[row][col];
	}
}


class ACSAnt4                //蚂蚁个体
{
private:
	ACS4* antColony;   //蚁群
protected:
	int startCity, cururentCity;//初始城市编号，当前城市编号
	int allowed[N];//每个蚂蚁维护自己的禁忌表	
	int Tour[N][2];//当前路径，是一个个路径段序列组成，即（currentcity，nextcity），用(Tour[i][0],Tour[i][1])表示
	int currentTourIndex;//当前路径索引，从0开始，存储蚂蚁经过城市的编号
public:
	ACSAnt4(ACS4* acs, int start)  //构造函数
	{
		antColony = acs;
		startCity = start;
	}
	int* Search4();	//开始搜索	
	int Choose4();	//选择下一节点	
	void ant_numoveToNextCity4(int nextCity);	//移动到下一节点
	~ACSAnt4()
	{}; //析构函数
};


int* ACSAnt4::Search4()  //开始搜索
{
	cururentCity = startCity;
	int toCity;
	currentTourIndex = 0;        //当前路径索引，存储蚂蚁经过城市的编号
	for (int i = 0; i < N; i++)
	{
		allowed[i] = 1;          //禁忌表初始化，1表示允许访问
	}
	allowed[cururentCity] = 0;	//cururentCity为当前城市编号，禁忌之
	int endCity;
	int count = 0;
	do
	{
		count++;
		endCity = cururentCity;
		toCity = Choose4();	     //选择下一个节点	
		if (toCity >= 0)
		{
			ant_numoveToNextCity4(toCity);    //移动到下一个节点
			antColony->UpdateLocalPathRule4(endCity, toCity);  //进行局部更新
			cururentCity = toCity;  //更新蚂蚁当前城市编号
		}
	} while (toCity >= 0);
	ant_numoveToNextCity4(startCity);
	antColony->UpdateLocalPathRule4(endCity, startCity);

	return *Tour;  //Tour是一个二维数组，Tour表示首元素地址的地址
				   /*
				   tourPath为指向int数的指针，相当于一维数组tourpath[]；
				   tourpath=*Tour，即将二维数组Tour首元素地址给tourpath；
				   所以tourpath[0]=Tour首元素；
				   tourpath[]={Tour[0][0],Tour[0][1],Tour[1][0],Tour[1][1],...Tour[74][0],Tour[74][1]}
				   tourpath下标：   0			1			2		3				148			149
				   对应路径序列：第1段路径:（*（tourpath），*（tourpath+1））
				   ...
				   第i段路径：（*（tourpath+2*（i-1）），*（tourpath+2*（i-1）+1））
				   */
}


int ACSAnt4::Choose4()  //选择下一节点，返回下一节点标号nextCity
{
	int nextCity = -1;
	double q = rand() / (double)RAND_MAX;    //产生一个0~1之间的随机数q												 
	if (q <= q0)	//如果 q <= q0,按先验知识，否则则按概率转移
	{
		double probability = -1.0;//转移到下一节点的概率
		for (int i = 0; i < N; i++)
		{
			//去掉禁忌表中已走过的节点,从剩下节点中选择最大概率的可行节点
			if (1 == allowed[i])
			{
				double prob = antColony->Transition4(cururentCity, i);  //计算当前节点转移到下一节点的概率
				if (prob  > probability)
				{
					nextCity = i;
					probability = prob;
				}
			}
		}
	}
	else
	{
		//按概率转移			
		double p = rand() / (double)RAND_MAX;	//生成一个随机数,用来判断落在哪个区间段
		double sum = 0.0;  //
		double probability = 0.0;	//概率的区间点，p 落在哪个区间段，则该点是转移的方向										
		for (int i = 0; i < N; i++)	//计算概率公式的分母的值
		{
			if (1 == allowed[i])  //选择未禁忌的城市
			{
				sum += antColony->Transition4(cururentCity, i);
			}
		}
		for (int j = 0; j < N; j++)
		{
			if (1 == allowed[j] && sum > 0)
			{
				probability += antColony->Transition4(cururentCity, j) / sum; //往城市j转移的概率
				if (probability >= p || (p > 0.9999 && probability > 0.9999))
				{
					nextCity = j;
					break;
				}
			}
		}
	}
	return nextCity;
}


void ACSAnt4::ant_numoveToNextCity4(int nextCity) //移动到下一节点
{
	allowed[nextCity] = 0;    //禁忌表，0表示经过，禁忌之
	Tour[currentTourIndex][0] = cururentCity;	//当前路径
	Tour[currentTourIndex][1] = nextCity;
	currentTourIndex++;
	cururentCity = nextCity;
}


double CalAdjacentDistance4(int node)//辅助函数，给一个节点由最近邻距离方法计算长度（计算Cm1）
{
	double sum = 0.0;
	int visitedNode[N];
	for (int j = 0; j < N; j++)
	{
		visitedNode[j] = 1;
	}
	visitedNode[node] = 0;
	int currentNode = node;
	int nextNode;
	do
	{
		nextNode = ChooseNextNode4(currentNode, visitedNode);
		if (nextNode >= 0)
		{
			sum += Mutual_Distance4[currentNode][nextNode];  //Mutual_Distance1为两城市间的距离矩阵
			currentNode = nextNode;
			visitedNode[currentNode] = 0;
		}
	} while (nextNode >= 0);
	sum += Mutual_Distance4[currentNode][node];
	return sum;
}


void Process_eil101()
{
	double start_time, end_time, run_time, sum_time = 0.0,total_length = 0;
	char filename[20] = "eil101";
	char filename1[20] = "eil101";
	Read_Coordinate4(filename);
	Read_opt_path4(filename1);
	CalculateDistance4();	//计算表示两两城市之间的距离
	FILE *total_result;
	if ((total_result = fopen("total_result.txt", "a")) == NULL) {  //fopen(文件名，文件使用方式)，向文件名为total_result.txt的文件中追加（append）内容，若出错fopen返回空指针，则返回报错信息
		printf("Oops,open files error..");
		exit(0);
	}
	fprintf(total_result, "\n\neil101测试结果\ntest_time\trun_time\tbest_result\n");  //追加内容
	for (int test_num = 0; test_num < 30; test_num++)
	{
		start_time = clock();
		srand((unsigned)time(NULL));
		double globalBestLength = 0.0;	 //全局最优长度
		ACS4* acs = new ACS4();	//蚁群系统对象
		ACSAnt4* ants[ANT_NUM];
		for (int k = 0; k < ANT_NUM; k++)	//每只蚂蚁随机选择一个出发城市
		{
			ants[k] = new ACSAnt4(acs, (int)(rand() % N));
		}
		int node = rand() % N;	//随机选择一个节点计算由贪心算法得到的一个长度
		Cm4 = CalAdjacentDistance4(node);
		double initInfo = 1 / (N * Cm4);	//各条路径上初始化的信息素强浓度
		acs->InitParameter4(initInfo);	  //初始化各路径信息素浓度

		int globalTour[N][2];        //全局最优路径，就是路径序列								
		for (int i = 0; i < MAX_ITER; i++)    // MAX_ITER最大循环次数
		{
			int localTour[N][2];	//局部最优路径	
			double localBestLength = 0.0;	//局部最优长度	
			double tourLength;	//当前路径长度
			for (int j = 0; j < ANT_NUM; j++)  //每个蚂蚁并行构建
			{
				int* tourPath = ants[j]->Search4();
				tourLength = calculateSumOfDistance4(tourPath);
				//局部比较，并记录路径和长度
				if (tourLength < localBestLength || abs(localBestLength - 0.0) < 0.000001)
				{
					for (int m = 0; m < N; m++)
					{
						int row = *(tourPath + 2 * m);
						int col = *(tourPath + 2 * m + 1);
						localTour[m][0] = row;
						localTour[m][1] = col;
					}
					localBestLength = tourLength;
				}
			}

			//全局比较，并记录路径和长度
			if (localBestLength < globalBestLength || abs(globalBestLength - 0.0) < 0.000001)
			{
				for (int m = 0; m < N; m++)
				{
					globalTour[m][0] = localTour[m][0];
					globalTour[m][1] = localTour[m][1];
				}
				globalBestLength = localBestLength;
			}
			acs->UpdateGlobalPathRule4(*globalTour, globalBestLength);
		}

		//输出全局最优路径
		end_time = clock();
		run_time = (end_time - start_time) / (double)CLOCKS_PER_SEC;
		sum_time += run_time;
		total_length += globalBestLength;
		fprintf(total_result, "%d\t%5.3fS\t%3.0f\n", test_num + 1, run_time, globalBestLength);
		cout << "全局最优路径长度:" << globalBestLength << endl;
		cout << "全局最优路径:";
		for (int m = 0; m < N; m++)
		{
			cout << globalTour[m][0] + 1 << "	";
		}
		cout << endl << endl;
		delete acs; //释放内存
	}

	fprintf(total_result, "平均值:   \t%5.3fS\t%3.0f\n", sum_time / 30, total_length / 30);
	fprintf(total_result, "已公布最佳路径长度最好结果：629\n");
	fprintf(total_result, "已公布最佳路线：\n");
	for (int d = 0; d < N; d++)
	{
		if ((d + 1) % 10 == 0)
			fprintf(total_result, "%d\n", best_path4[d]);
		else
			fprintf(total_result, "%d\t", best_path4[d]);
	}
	fclose(total_result);
	Print_bestpath4();
}


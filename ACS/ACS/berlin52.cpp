#include <iostream>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const int N = 52;
using namespace::std;
int coordinate1[N][2];
int best_path1[N];


const int MAX_ITER = 2000;
const int ANT_NUM = 10;  //��������
const double Alpha = 1.0, Beta = 5.0, Rho = 0.1, q0 = 0.9, Xi = 0.1;
//Alpha:��Ϣ��Ȩ��,����Ϊ1��ʾACS��Beta:����ʽ��ϢȨ�أ�Rho:ȫ����Ϣ�ػӷ�������Xi���ֲ���Ϣ�ػӷ����ӣ�q0:α�������q0

double Mutual_Distance1[N][N];	//��ʾ��������֮��ľ���
double Cm1;		//��Ϣ�ؾֲ�����ʱ��ʹ�õĵĳ�������̰���㷨�õ�

void Read_opt_path1(char a[])  //��ȡ����������·��
{
	a = strcat(a, ".opt.tour");
	FILE* fp;
	char str[100];  //�ݴ��ȡ���ַ���
	int i0 = 1, j0 = 0; //i0���ƴ�����·���ļ��ڼ��ж�ȡ��j0Ϊ����·����ֵ�±�
	int i = 1, j = 0, m = 0;  //i���ƴӳ��������ļ��ڼ��ж�ȡ��j����ֻ������ֵ���������б��,mΪ�������긳ֵ�±�
	fp = fopen(a, "r");
	while (i0 < 6)
	{
		fgets(str, 255, fp);
		i0++;
	}
	while (!feof(fp) && j0 < N)
	{
		fscanf(fp, "%s\n", str);
		best_path1[j0] = atoi(str);
		j0++;
	}
	fclose(fp);
}

void Read_Coordinate1(char a[])   //��ȡ��������
{
	a = strcat(a, ".tsp");
	FILE* fp;
	char str[100];  //�ݴ��ȡ���ַ���
	int i = 1, j = 0, m = 0;  //i���ƴӳ��������ļ��ڼ��ж�ȡ��j����ֻ������ֵ���������б��,mΪ�������긳ֵ�±�
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
			coordinate1[m][0] = atoi(str);
		}
		else if (j % 3 == 2) {
			coordinate1[m][1] = atoi(str);
			m++;
		}
		j++;
	}
	fclose(fp);
}

void Print_bestpath1()  //��ӡ�ѹ��������Ž�·��������
{
	int i;
	cout << "�ѹ�������·��Ϊ��" << endl;
	for (i = 0; i < N - 1; i++)
	{
		cout << best_path1[i] << " ->";
	}
	cout << best_path1[N - 1] << endl;
	cout << "���Ž���С·������Ϊ��7542" << endl<<endl<<endl;
}

int calculateDistance1(int i, int j)//������������֮��ľ���(ŷ�Ͼ���)
{
	return (int)(sqrt((coordinate1[i][0] - coordinate1[j][0])*(coordinate1[i][0] - coordinate1[j][0]) + (coordinate1[i][1] - coordinate1[j][1])*(coordinate1[i][1] - coordinate1[j][1])) + 0.5);
}



void CalculateDistance1()   //�þ����ʾ��������֮��ľ���
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (i != j)
			{
				Mutual_Distance1[i][j] = calculateDistance1(i, j);
				Mutual_Distance1[j][i] = Mutual_Distance1[i][j];
			}
		}
	}
}


double calculateSumOfDistance1(int* tour) //��þ���n�����е�·������  
{
	double sum = 0;
	for (int i = 0; i< N; i++)
	{
		int row = *(tour + 2 * i);
		int col = *(tour + 2 * i + 1);
		sum += Mutual_Distance1[row][col];
	}
	return sum;
}


int ChooseNextNode1(int currentNode, int visitedNode[])
{
	int nextNode = -1;
	double shortDistance = 0.0;
	for (int i = 0; i < N; i++)
	{
		if (1 == visitedNode[i])	//ȥ�����߹��Ľڵ�,��ʣ�½ڵ���ѡ���������Ľڵ�
		{
			if (shortDistance == 0.0)
			{
				shortDistance = Mutual_Distance1[currentNode][i];
				nextNode = i;
			}
			if (shortDistance > Mutual_Distance1[currentNode][i])
			{
				shortDistance = Mutual_Distance1[currentNode][i];
				nextNode = i;
			}
		}
	}
	return nextNode;
}//ѡ����һ���ڵ㣬�������ĺ���������ĳ���


double CalAdjacentDistance1(int node)//������������һ���ڵ�������ھ��뷽�����㳤�ȣ�����Cm1��
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
		nextNode = ChooseNextNode1(currentNode, visitedNode);
		if (nextNode >= 0)
		{
			sum += Mutual_Distance1[currentNode][nextNode];  //Mutual_Distance1Ϊ�����м�ľ������
			currentNode = nextNode;
			visitedNode[currentNode] = 0;
		}
	} while (nextNode >= 0);
	sum += Mutual_Distance1[currentNode][node];
	return sum;
}



class ACS1    //��Ⱥϵͳ
{
private:
	double pheromone[N][N], visible[N][N];//�ڵ�֮�����Ϣ�����ͽڵ�֮�������ʽ��Ϣ��(�������ĵ���)
public:
	ACS1()
	{}  //���캯��
	void InitParameter1(double value);	//��ʼ����·����Ϣ��	
	double Transition1(int i, int j);	//���㵱ǰ�ڵ㵽��һ�ڵ�ת�Ƶĸ���
	void UpdateLocalPathRule1(int i, int j);	//�ֲ����¹���	
	void UpdateGlobalPathRule1(int* bestTour, int globalBestLength);	//ȫ�ָ��¹���
	~ACS1()
	{};  //��������
};


void ACS1::InitParameter1(double value) 	//��ʼ����·���ϵ���Ϣ��ǿ��Ϊvalue
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			pheromone[i][j] = value;   //��Ϣ�س�ʼǿ��
			pheromone[j][i] = value;
			if (i != j)
			{
				visible[i][j] = 1.0 / Mutual_Distance1[i][j];   //����ʽ��Ϣǿ��
				visible[j][i] = visible[i][j];
			}
		}
	}
}


double ACS1::Transition1(int i, int j)//���㵱ǰ�ڵ㵽��һ�ڵ�ת�Ƶĸ���
{
	if (i != j)
	{
		return (pow(pheromone[i][j], Alpha) * pow(visible[i][j], Beta));  //���õ�5.1�ĸ��ʼ��㹫ʽ
	}
	else
	{
		return 0.0;
	}
}


void ACS1::UpdateLocalPathRule1(int i, int j) //�ֲ����¹���
{
	pheromone[i][j] = (1.0 - Xi) * pheromone[i][j] + Xi * (1.0 / (N * Cm1));
	pheromone[j][i] = pheromone[i][j];
}


void ACS1::UpdateGlobalPathRule1(int* bestTour, int globalBestLength)//ȫ�ָ��¹���
{
	for (int i = 0; i < N; i++)
	{
		int row = *(bestTour + 2 * i);
		int col = *(bestTour + 2 * i + 1);
		pheromone[row][col] = (1.0 - Rho) * pheromone[row][col] + Rho * (1.0 / globalBestLength);//ֻ�ڵ�ǰ����·��������Ϣ��
		pheromone[col][row] = pheromone[row][col];
	}
}


class ACSAnt1                //���ϸ���
{
private:
	ACS1* antColony;   //��Ⱥ
protected:
	int startCity, cururentCity;//��ʼ���б�ţ���ǰ���б��
	int allowed[N];//ÿ������ά���Լ��Ľ��ɱ�	
	int Tour[N][2];//��ǰ·������һ����·����������ɣ�����currentcity��nextcity������(Tour[i][0],Tour[i][1])��ʾ
	int currentTourIndex;//��ǰ·����������0��ʼ���洢���Ͼ������еı��
public:
	ACSAnt1(ACS1* acs, int start)  //���캯��
	{
		antColony = acs;
		startCity = start;
	}
	int* Search1();	//��ʼ����	
	int Choose1();	//ѡ����һ�ڵ�	
	void ant_numoveToNextCity1(int nextCity);	//�ƶ�����һ�ڵ�
	~ACSAnt1()
	{}; //��������
};


int* ACSAnt1::Search1()  //��ʼ����
{
	cururentCity = startCity;
	int toCity;
	currentTourIndex = 0;        //��ǰ·���������洢���Ͼ������еı��
	for (int i = 0; i < N; i++)
	{
		allowed[i] = 1;          //���ɱ���ʼ����1��ʾ��������
	}
	allowed[cururentCity] = 0;	//cururentCityΪ��ǰ���б�ţ�����֮
	int endCity;
	int count = 0;
	do
	{
		count++;
		endCity = cururentCity;
		toCity = Choose1();	     //ѡ����һ���ڵ�	
		if (toCity >= 0)
		{
			ant_numoveToNextCity1(toCity);    //�ƶ�����һ���ڵ�
			antColony->UpdateLocalPathRule1(endCity, toCity);  //���оֲ�����
			cururentCity = toCity;  //�������ϵ�ǰ���б��
		}
	} while (toCity >= 0);
	ant_numoveToNextCity1(startCity);
	antColony->UpdateLocalPathRule1(endCity, startCity);

	return *Tour;  //Tour��һ����ά���飬Tour��ʾ��Ԫ�ص�ַ�ĵ�ַ
				   /*
				   tourPathΪָ��int����ָ�룬�൱��һά����tourpath[]��
				   tourpath=*Tour��������ά����Tour��Ԫ�ص�ַ��tourpath��
				   ����tourpath[0]=Tour��Ԫ�أ�
				   tourpath[]={Tour[0][0],Tour[0][1],Tour[1][0],Tour[1][1],...Tour[74][0],Tour[74][1]}
				   tourpath�±꣺   0			1			2		3				148			149
				   ��Ӧ·�����У���1��·��:��*��tourpath����*��tourpath+1����
				   ...
				   ��i��·������*��tourpath+2*��i-1������*��tourpath+2*��i-1��+1����
				   */
}


int ACSAnt1::Choose1()  //ѡ����һ�ڵ㣬������һ�ڵ���nextCity
{
	int nextCity = -1;
	double q = rand() / (double)RAND_MAX;    //����һ��0~1֮��������q												 
	if (q <= q0)	//��� q <= q0,������֪ʶ�������򰴸���ת��
	{
		double probability = -1.0;//ת�Ƶ���һ�ڵ�ĸ���
		for (int i = 0; i < N; i++)
		{
			//ȥ�����ɱ������߹��Ľڵ�,��ʣ�½ڵ���ѡ�������ʵĿ��нڵ�
			if (1 == allowed[i])
			{
				double prob = antColony->Transition1(cururentCity, i);  //���㵱ǰ�ڵ�ת�Ƶ���һ�ڵ�ĸ���
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
		//������ת��			
		double p = rand() / (double)RAND_MAX;	//����һ�������,�����ж������ĸ������
		double sum = 0.0;  //
		double probability = 0.0;	//���ʵ�����㣬p �����ĸ�����Σ���õ���ת�Ƶķ���										
		for (int i = 0; i < N; i++)	//������ʹ�ʽ�ķ�ĸ��ֵ
		{
			if (1 == allowed[i])  //ѡ��δ���ɵĳ���
			{
				sum += antColony->Transition1(cururentCity, i);
			}
		}
		for (int j = 0; j < N; j++)
		{
			if (1 == allowed[j] && sum > 0)
			{
				probability += antColony->Transition1(cururentCity, j) / sum; //������jת�Ƶĸ���
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


void ACSAnt1::ant_numoveToNextCity1(int nextCity) //�ƶ�����һ�ڵ�
{
	allowed[nextCity] = 0;    //���ɱ���0��ʾ����������֮
	Tour[currentTourIndex][0] = cururentCity;	//��ǰ·��
	Tour[currentTourIndex][1] = nextCity;
	currentTourIndex++;
	cururentCity = nextCity;
}



void Process_berlin52()
{
	double start_time, end_time, run_time, sum_time = 0.0,total_length = 0;
	char filename[20] = "berlin52";
	char filename1[20] = "berlin52";
	Read_Coordinate1(filename);
	Read_opt_path1(filename1);
	CalculateDistance1();	//�����ʾ��������֮��ľ���
	FILE *total_result;
	if ((total_result = fopen("total_result.txt", "a")) == NULL) {  //fopen(�ļ������ļ�ʹ�÷�ʽ)�����ļ���Ϊtotal_result.txt���ļ���׷�ӣ�append�����ݣ�������fopen���ؿ�ָ�룬�򷵻ر�����Ϣ
		printf("Oops,open files error..");
		exit(0);
	}
	fprintf(total_result, "\n\nberlin52���Խ��\ntest_time\trun_time\tbest_result\n");  //׷������
	for (int test_num = 0; test_num < 30; test_num++)
	{
		start_time = clock();
		srand((unsigned)time(NULL));
		double globalBestLength = 0.0;	 //ȫ�����ų���
		ACS1* acs = new ACS1();	//��Ⱥϵͳ����
		ACSAnt1* ants[ANT_NUM];
		for (int k = 0; k < ANT_NUM; k++)	//ÿֻ�������ѡ��һ����������
		{
			ants[k] = new ACSAnt1(acs, (int)(rand() % N));
		}
		int node = rand() % N;	//���ѡ��һ���ڵ������̰���㷨�õ���һ������
		Cm1 = CalAdjacentDistance1(node);
		double initInfo = 1 / (N * Cm1);	//����·���ϳ�ʼ������Ϣ��ǿŨ��
		acs->InitParameter1(initInfo);	  //��ʼ����·����Ϣ��Ũ��

		int globalTour[N][2];        //ȫ������·��������·������								
		for (int i = 0; i < MAX_ITER; i++)    // MAX_ITER���ѭ������
		{
			int localTour[N][2];	//�ֲ�����·��	
			double localBestLength = 0.0;	//�ֲ����ų���	
			double tourLength;	//��ǰ·������
			for (int j = 0; j < ANT_NUM; j++)  //ÿ�����ϲ��й���
			{
				int* tourPath = ants[j]->Search1();
				tourLength = calculateSumOfDistance1(tourPath);
				//�ֲ��Ƚϣ�����¼·���ͳ���
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

			//ȫ�ֱȽϣ�����¼·���ͳ���
			if (localBestLength < globalBestLength || abs(globalBestLength - 0.0) < 0.000001)
			{
				for (int m = 0; m < N; m++)
				{
					globalTour[m][0] = localTour[m][0];
					globalTour[m][1] = localTour[m][1];
				}
				globalBestLength = localBestLength;
			}
			acs->UpdateGlobalPathRule1(*globalTour, globalBestLength);
		}

		//���ȫ������·��
		end_time = clock();
		run_time = (end_time - start_time) / (double)CLOCKS_PER_SEC;
		sum_time += run_time;
		total_length += globalBestLength;
		fprintf(total_result, "%d\t%4.3fS\t%4.0f\n", test_num + 1, run_time, globalBestLength);
		cout << "ȫ������·������:" << globalBestLength << endl;
		cout << "ȫ������·��:";
		for (int m = 0; m < N; m++)
		{
			cout << globalTour[m][0] + 1 << "	";
		}
		cout << endl << endl;
		delete acs; //�ͷ��ڴ�
	}

	fprintf(total_result, "ƽ��ֵ:   \t%4.3fS\t%4.0f\n", sum_time / 30, total_length / 30);
	fprintf(total_result, "�ѹ������·��������ý����7542\n");
	fprintf(total_result, "�ѹ������·�ߣ�\n");
	for (int d = 0; d < N; d++)
	{
		if ((d + 1) % 10 == 0)
			fprintf(total_result, "%d\n", best_path1[d]);
		else
			fprintf(total_result, "%d\t", best_path1[d]);
	}
	fclose(total_result);
	Print_bestpath1();
}

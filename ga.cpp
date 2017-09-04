//
// Created by zhenghui on 2017/6/5.
//
#include <iostream>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <cstring>
#include <iomanip>
#include <queue>
#pragma comment(linker, "/STACK:102400000,102400000");
using namespace std;

const int city_num = 144;          //城市数目
const int individual_num = 2000;   //种群初始值
const int age = 500;              //遗传代数
const double varation_p = 0.1;     //变异因子
typedef struct City_xy          //储存给定的城市数据包含序号和坐标
{
    int order;
    double x;
    double y;
}City_xy;
typedef struct City
{
    int mark;
}City;
typedef struct Population
{
    City city[145];              //城市数组
    double distance;             //个体城市序列距离之和
    double Fitness;              //适应度
    double Fitness_pi;
}Population;

bool  cmp(Population aa,Population bb){
    return aa.Fitness_pi > bb.Fitness_pi;
}
//变量
City_xy city_xy[city_num+1];                //城市数据
Population population[individual_num+1];   //种群中的个体
int randp[city_num+1];                    //随机数序列
int generation;                          //代数
queue <Population> best;                //记录最优

//方法
//void debug();                                 //调试
//void random_array(int n);                     //生成一个长度1-n的随机序列
//double Distance(int ax,int ay,int bx,int by); //求两个城市间的距离
//void  init_Population();                      //初始化种群
//double City_distance(Population &a)           //计算个体城市序列距离之和
//double adaptation(Population &a)              //计算适应度
//void adaptation_rate()                        //适应率---个体适应度占种群的百分比
//accumlation_rate()                            //累积率
//void TSP::Roulette()                          //轮盘赌算法，用于筛选个体
//void Intersect(Population &a, Population &b)  //交叉变异
//void Intersect_all()                          //种群所有个体交叉遗传
//void Varation()                               //变异



void random_array(int n) { //生成一个长度1-n的随机序列
    for (int i = 0; i <= n; ++i) {
        randp[i] = i;
    }
    for (int i = 1; i <= n; ++i) {
        int j = rand() % (n+1);
        if(j != 0)swap(randp[i],randp[j]);
    }
    for (int i = 1; i <= n; ++i) {
        int j = rand() % (n+1);
        if(j != 0)swap(randp[i],randp[j]);
    }
}
//初始化种群
void init_Population() {
    for (int i = 1; i <= individual_num ; ++i) {
        population[i].Fitness    = 0.0;
        population[i].Fitness_pi = 0.0;
        population[i].distance   = 0.0;
        population[i].Fitness_pi = 0.0;
    }
    generation = 0;
    for (int i = 1; i <= individual_num ; ++i) {
        random_array(city_num);
        for (int j = 1; j <=city_num ; ++j) {
            population[i].city[j].mark = randp[j];
        }
    }
}
//距离
double Distance(int ax, int ay, int bx, int by) {
    return (double)sqrt((ax-bx)*(ax-bx)+(ay-by)*(ay-by));
}
double City_distance(Population &a){
    double d;
    a.distance = 0;
    for (int i = 1; i < city_num; ++i) {
        d = Distance(city_xy[(a.city[i].mark)].x,city_xy[(a.city[i].mark)].y, city_xy[(a.city[i+1].mark)].x,city_xy[(a.city[i+1].mark)].y);
        a.distance += d;
    }
    a.distance += Distance(city_xy[(a.city[1].mark)].x,city_xy[(a.city[1].mark)].y, city_xy[(a.city[city_num].mark)].x,city_xy[(a.city[city_num].mark)].y);
    return a.distance;
}
//计算适应度
double adaptation(Population &a){
    double sum = City_distance(a);
    a.Fitness =  ((double)1/sum )*100000;
    return a.Fitness;
}
//适应率---个体适应度占种群的百分比
void adaptation_rate(){
    double d = 0.0;
    for (int i = 1; i <= individual_num; ++i) {
        d += adaptation(population[i]);
    }
    for (int i = 1; i <= individual_num; ++i) {
        population[i].Fitness_pi = (double)population[i].Fitness/(double)d*100;
    }
}

void Copy(Population &a, Population &b)//把a复制给b
{
    b.Fitness_pi = a.Fitness_pi;
    b.Fitness = a.Fitness;
    b.distance = a.distance;
    for (int i = 1; i <= city_num; i++) {
        b.city[i].mark = a.city[i].mark;

    }
}

void sort_pop(){ //根据个体适应率排序
    sort(population+1,population+individual_num+1,cmp);
}

void choose(){
    for (int i = 1; i <= individual_num/2; ++i) {
        Copy(population[i],population[individual_num+1-i]);
    }
}

//交叉
void Intersect(Population &a, Population &b) {
    int k1 = 1 + rand() % (city_num/2);
    int k2 = k1 + rand() % (city_num/2);
    for (int i = k1; i <= k2; ++i) {
        swap(a.city[i].mark,b.city[i].mark);
    }
    int a1[city_num+1],b1[city_num+1];
    memset(a1,0, sizeof(a1));
    memset(b1,0, sizeof(b1));
    for (int i = 1; i <= city_num; ++i) {
        a1[a.city[i].mark]++;
        b1[b.city[i].mark]++;
    }
    for (int i = 1; i <= city_num; ++i) {
        if(i<k1 || i>k2){
            if(a1[a.city[i].mark]>1){
                for (int j = 1; j <= city_num; ++j) {
                    if(j<k1 || j>k2){
                        if(b1[b.city[j].mark]>1){
                            swap(a.city[i].mark,b.city[j].mark);
                            break;
                        }
                    }
                }
            }
        }
    }
}
//种群所有个体交叉遗传
void Intersect_all(){
    for (int i = 1; i <= individual_num/2; i++) {
        Intersect(population[i],population[individual_num/2+i]);
    }
}

//变异
void Varation(){
    for (int i = 1; i <= individual_num; ++i) {
        int t1=rand()%100+1;
        int t2=rand()%(t1+1);
        if((double)(t2/t1) > varation_p){
            for (int j = 0; j < 6; ++j) {
                int k1 = 1 + rand() % (city_num/2);
                int k2 = k1 + rand() % (city_num/2);
                swap(population[i].city[k1].mark,population[i].city[k2].mark);
            }
        }
    }

}


void Reset()
{
    for (int i = 1; i <= individual_num; i++)
    {
        population[i].Fitness = 0.0;
        population[i].Fitness_pi = 0.0;
    }
}


//调试函数
void debug() {
    printf("\n-----------\n");
    for (int i = 1; i <= individual_num; ++i) {
        printf("%f\n",population[i].Fitness_pi);
    }

}

void run_ga() {
    int this_age = age;
    init_Population();
    adaptation_rate();
    sort_pop();
    best.push(population[1]);

    while(this_age--){
        choose();
        Intersect_all();
        Varation();
        adaptation_rate();
        sort_pop();
        best.push(population[1]);
    }
    //debug();
}

int main() {
    freopen("G:\\C++\\project\\algorithm\\city.in","r",stdin);
    freopen("G:\\C++\\project\\algorithm\\ga.out","w",stdout);
    srand((int)time(NULL));
    for (int i = 1; i <= city_num; ++i) {
        scanf("%d %lf %lf",&city_xy[i].order,&city_xy[i].x,&city_xy[i].y);
    }
    run_ga();
    while (!best.empty()){
        //Population temp = best.front();
        printf("%f\n",best.front().distance);
        best.pop();
    }
    return 0;
}

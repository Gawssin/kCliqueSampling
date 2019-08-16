#include <stdlib.h>
#include <cstdio>
#include <stdbool.h>
#include <string.h>
#include <string>
#include <vector>
#include <utility>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <random>
#include <cmath>
#include <time.h>

#define NLINKS 1000000
using namespace std;

class edge
{
public:
	int s;
	int t;
};

class iddeg
{
public:
	int id;
	int degree;
};

class Graph {
public:
	Graph();
	Graph(const Graph& obj);
	~Graph();
	void readedgelist(string edgelist);
	void coreDecomposition();
	void mkGraph();
	int outLargeClique();
	bool isEdge(int, int);
	void mksub(Graph&, int*, int);
	int color(int*);
	void kClique(int, long long*, long long*);
	void kCliqueCount(int, long long*, int*, int*, int*, int*, int*, int**, int**, long long*);
	//bool isEdge(int, int);

	int n;
	int e;
	int maxDeg;
	vector<edge> edges;

	int* deg;
	int* cd;
	int* adj;
	int* coreRank;			//increasing core number order
	int* coreNum;			//coreNum[i] is the core number of node i.
	int* bin;
};

inline int max3(int a, int b, int c) {
	a = (a > b) ? a : b;
	return (a > c) ? a : c;
}

Graph::Graph(void) {}
Graph::~Graph(void)
{
	if (deg != NULL) delete[] deg;
	if (cd != NULL) delete[] cd;
	if (adj != NULL) delete[] adj;
	if (coreRank != NULL) delete[] coreRank;
	if (coreNum != NULL) delete[] coreNum;
	if (bin != NULL) delete[] bin;
}
Graph::Graph(const Graph& obj)
{
	n = obj.n, e = obj.e, maxDeg = obj.maxDeg, edges = obj.edges;
	edges = obj.edges;

	if (deg != NULL) delete[] deg;
	if (obj.deg != NULL) deg = new int[n], memcpy(deg, obj.deg, n * sizeof(int));

	if (cd != NULL) delete[] cd;
	if (obj.cd != NULL) cd = new int[n], memcpy(cd, obj.cd, n * sizeof(int));

	if (adj != NULL) delete[] adj;
	if (obj.adj != NULL) adj = new int[2 * e], memcpy(adj, obj.adj, 2 * e * sizeof(int));

	if (coreRank != NULL) delete[] coreRank;
	if (obj.coreRank != NULL) coreRank = new int[n], memcpy(coreRank, obj.coreRank, n * sizeof(int));

	if (coreNum != NULL) delete[] coreNum;
	if (obj.coreNum != NULL) coreNum = new int[n], memcpy(coreNum, obj.coreNum, n * sizeof(int));

	if (bin != NULL) delete[] bin;
	if (obj.bin != NULL) bin = new int[maxDeg + 2], memcpy(bin, obj.bin, (maxDeg + 2) * sizeof(int));
}

void Graph::readedgelist(string edgelist) {

	int e1 = NLINKS;
	n = 0;
	e = 0;
	edges.resize(e1);
	ifstream file;
	file.open(edgelist);

	while (file >> edges[e].s >> edges[e].t)
	{
		n = max3(n, edges[e].s, edges[e].t);
		e++;
		if (e == e1) {
			e1 += NLINKS;
			edges.resize(e1);
		}
	}
	file.close();
	n++;
	edges.resize(e);
}

void Graph::mkGraph()
{
	deg = new int[n]();
	cd = new int[n + 1];
	adj = new int[2 * e];
	maxDeg = 0;
	for (int i = 0; i < e; i++)
	{
		deg[edges[i].s]++;
		deg[edges[i].t]++;
		maxDeg = max3(maxDeg, deg[edges[i].s], deg[edges[i].t]);
	}
	cd[0] = 0;
	for (int i = 1; i < n + 1; i++) {
		cd[i] = cd[i - 1] + deg[i - 1];
		deg[i - 1] = 0;
	}

	for (int i = 0; i < e; i++) {
		adj[cd[edges[i].s] + deg[edges[i].s]++] = edges[i].t;
		adj[cd[edges[i].t] + deg[edges[i].t]++] = edges[i].s;
	}

	for (int i = 0; i < n; i++) sort(adj + cd[i], adj + cd[i] + deg[i]);
}

bool cmp(const pair<int, int>& a, const pair<int, int>& b)
{
	return a.second > b.second;
}
bool IGCmp(const iddeg& a, const iddeg& b)
{
	return a.degree == b.degree ? (a.id < b.id) : (a.degree > b.degree);
}

bool Graph::isEdge(int a, int b)
{
	if (deg[a] > deg[b]) a = a ^ b, b = a ^ b, a = a ^ b;
	for (int i = cd[a]; i < cd[a] + deg[a]; i++)
		if (adj[i] == b) return true;
	return false;
}
int Graph::outLargeClique()
{
	int CSize = 0;
	for (int i = n - 1; i >= 0; i--)
	{
		int id = coreRank[i];
		if (coreNum[id] >= CSize)
		{
			pair<int, int>* SCore = new pair<int, int>[deg[id]];
			//int *S = new int[deg[id]], cnt = 0;
			int cnt = 0, ind = 0;
			for (int j = cd[id]; j < cd[id] + deg[id]; j++)
				if (coreNum[adj[j]] >= CSize)
				{
					SCore[cnt].first = adj[j];
					SCore[cnt].second = coreNum[adj[j]];
					cnt++;
				}
			//S[cnt++] = adj[j];
			sort(SCore, SCore + cnt, cmp);
			//sort(S, S + cnt, cmp);
			int* C = new int[deg[id]];
			for (int j = 0; j < cnt; j++)
			{
				int flag = 1;
				for (int k = 0; k < ind; k++)
				{
					if (isEdge(SCore[j].first, C[k]) == false)
					{
						flag = 0;
						break;
					}
				}
				if (flag) C[ind++] = SCore[j].first;
			}
			ind++;	//node "id" ?
			if (ind > CSize) CSize = ind;
		}
	}
	return CSize;
}

void Graph::coreDecomposition()
{
	bin = new int[maxDeg + 2]();

	for (int i = 0; i < n; i++)
		bin[deg[i]]++;

	int lastBin = bin[0], nowBin;
	bin[0] = 0;
	for (int i = 1; i <= maxDeg; i++)
	{
		nowBin = lastBin + bin[i - 1];
		lastBin = bin[i];
		bin[i] = nowBin;
	}
	int* vert = new int[n](), * pos = new int[n](), * tmpDeg = new int[n]();
	for (int i = 0; i < n; i++)
	{
		pos[i] = bin[deg[i]];

		vert[bin[deg[i]]++] = i;
		tmpDeg[i] = deg[i];
	}

	bin[0] = 0;
	for (int i = maxDeg; i >= 1; i--)
	{
		bin[i] = bin[i - 1];
	}

	//int core = 0;
	int* cNum = new int[n];
	//int *cNum = (int *)malloc(g->n * sizeof(int));
	for (int i = 0; i < n; i++)
	{
		int id = vert[i], nbr, binFrontId;
		//if (i == bin[core + 1]) ++core;
		cNum[id] = tmpDeg[id];
		for (int i = cd[id]; i < cd[id] + deg[id]; i++)
		{
			nbr = adj[i];

			if (tmpDeg[nbr] > tmpDeg[id])
			{
				binFrontId = vert[bin[tmpDeg[nbr]]];
				if (binFrontId != nbr)
				{

					pos[binFrontId] = pos[nbr];
					pos[nbr] = bin[tmpDeg[nbr]];
					vert[bin[tmpDeg[nbr]]] = nbr;
					vert[pos[binFrontId]] = binFrontId;

				}
				bin[tmpDeg[nbr]]++;
				tmpDeg[nbr]--;

			}

		}

	}

	coreNum = cNum;

	coreRank = vert;

	delete[] tmpDeg;
	delete[] pos;
}

void Graph::mksub(Graph& sg, int* nodes, int NodeNum)
{
	sg.n = NodeNum, sg.e = 0;
	int* newFg = new int[n]();

	for (int i = 0; i < NodeNum; i++) newFg[nodes[i]] = 1;

	for (int i = 0; i < e; i++)
		if (newFg[edges[i].s] == 1 && newFg[edges[i].t] == 1) sg.e++;

	//printf("sg.e = %d\n", sg.e);
	sg.edges.resize(sg.e);



	//sort(nodes, nodes+NodeNum);


	/*
	for (int i = 0; i < sg.n; i++)
	{
		int ind = cd[nodes[i]];
		for (int j = i+1; j < sg.n; j++)
		{
			while (ind < cd[nodes[i]] + deg[nodes[i]])
			{
				if (adj[ind] == nodes[j])
				{
					sg.e++;
					ind++;
					break;
				}
				else if (adj[ind] > nodes[j]) break;
				ind++;
			}
		}
	}
	printf("sg.e = %d\n", sg.e);
	sg.edges.resize(sg.e);
	*/

	sg.e = 0;
	int* lab = new int[n], cnt = 0;
	for (int i = 0; i < sg.n; i++) lab[nodes[i]] = -1;


	for (int i = 0; i < e; i++)
	{
		if (newFg[edges[i].s] == 1 && newFg[edges[i].t] == 1)
		{
			lab[edges[i].s] = (lab[edges[i].s] == -1 ? (cnt++) : lab[edges[i].s]);
			lab[edges[i].t] = (lab[edges[i].t] == -1 ? (cnt++) : lab[edges[i].t]);
			sg.edges[sg.e].s = lab[edges[i].s];
			sg.edges[sg.e].t = lab[edges[i].t];
			sg.e++;
		}
	}


	/*

	for (int i = 0; i < sg.n; i++)
	{
		int ind = cd[nodes[i]];
		for (int j = i + 1; j < sg.n; j++)
		{
			while (ind < cd[nodes[i]] + deg[nodes[i]])
			{
				if (adj[ind] == nodes[j])
				{
					//cout << "nodes[i] = " << nodes[i] << " nodes[j] = " << nodes[j] << endl;
					lab[nodes[i]] = (lab[nodes[i]] == -1 ? (++cnt) : lab[nodes[i]]);
					lab[adj[ind]] = (lab[adj[ind]] == -1 ? (++cnt) : lab[adj[ind]]);
					sg.edges[sg.e].s = lab[nodes[i]];
					sg.edges[sg.e].t = lab[adj[ind]];
					sg.e++;
					ind++;
					break;
				}
				else if (adj[ind] > nodes[j]) break;
				ind++;
			}
		}
	}

	*/
	//printf("sg labeled\n");

	delete[] newFg;
	delete[] lab;
	sg.mkGraph();
}

int Graph::color(int* color)
{
	iddeg* ig = new iddeg[n];
	for (int i = 0; i < n; i++)
	{
		ig[i].id = i;
		ig[i].degree = deg[i];
	}

	sort(ig, ig + n, IGCmp);

	//color = new int[n];
	memset(color, -1, sizeof(int) * n);
	int* C = new int[(ig[0].degree + 1)]();

	color[ig[0].id] = 0;
	int colorNum = 1;

	for (int i = 1; i < n; i++)
	{
		int tmpDeg = ig[i].degree, tmpid = ig[i].id;
		for (int j = 0; j < tmpDeg; j++)
		{
			int now = adj[cd[tmpid] + j];
			if (color[now] != -1)
				C[color[now]] = 1;
		}
		for (int j = 0; j < ig[0].degree + 1; j++)
			if (C[j] == 0)
			{
				color[ig[i].id] = j;
				colorNum = j > colorNum ? j : colorNum;
				break;
			}

		for (int j = 0; j < tmpDeg; j++)
		{
			int now = adj[cd[tmpid] + j];
			if (color[now] != -1)
				C[color[now]] = 0;
		}

	}
	//printf("color number = %d\n", colorNum);
	delete[] ig;
	delete[] C;
	return colorNum + 1;
}

int SelectPoint(double* wpi, int n)
{
	vector<double> wts(wpi, wpi + n);
	//cout << wts[0] << endl;

	random_device rd;  //Will be used to obtain a seed for the random number engine
	mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

	discrete_distribution<size_t> d{ begin(wts), end(wts) };
	return d(gen);

}

void selectNodes(int v, int tolCol, int* cNbrN, int* nC, int k, Graph* g, int* col, int** colSet, int* cols, double **ck)
{
	random_device rd;  //Will be used to obtain a seed for the random number engine
	mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	discrete_distribution<size_t> d;

	int cnt = 0, cptr = tolCol, kptr = k-1, state;

	while (cnt < k - 1)
	{
		d = {ck[cptr-1][kptr-1] * nC[cptr - 1], ck[cptr - 1][kptr] };
		state = d(gen);
		if (!state)
		{
			uniform_int_distribution<> ud(0, nC[cptr - 1] - 1);

			cNbrN[cnt++] = colSet[cptr - 1][ud(gen)];
			cptr--;
			kptr--;
		}
		else cptr--;

	}

}

bool isClique(Graph* g, int* cNbrN, int n)
{
	for (int i = 0; i < n; i++)
		for (int j = i + 1; j < n; j++)
			if (g->isEdge(cNbrN[i], cNbrN[j]) == false)
				return false;

	return true;
}

double calc(int n)
{
	double result = 1;
	for (int i = n; i >= 1; i--)
		result = result * i;
	return result;
}

int main(int argc, char** argv)
{
	random_device rd;  //    be used to obtain a seed for the random number engine
	mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	/*
	double a[5] = { 0.1,0.0,0.0,0.0,0.0 };
	//array<double, 6> ats;// = { 10.0, 10.0, 10.0, 10.0, 10.0, 30.0 };
	vector<double> wts(a, a + 5);
	cout << wts[0] << endl;



	discrete_distribution<size_t> d{ begin(wts), end(wts) };
	wts[4] = 0.1;



	d = { begin(wts), end(wts) };
	//uniform_real_distribution<> dis(1.0, 2.0);

	for (int i = 0; i < 5; i++) printf("wts[%d] = %f\n",i,wts[i]);
	//cout << "wts = " << wts[0] << wts[1] << wts[2] << wts[3] << wts[4] << endl;

	int result[5] = {0,0,0,0,0};
	for (int n = 0; n < 10000; ++n) {
		// Use dis to transform the random unsigned int generated by gen into a
		// double in [1, 2). Each call to dis(gen) generates a new random double
		//cout << "d gen = " << d(gen) << endl;
		result[d(gen)]++;
		//std::cout << d(gen) << endl;
	}

	for (int i = 0; i < 5; i++)
		cout << result[i] << endl;

	*/

	Graph* g = new Graph;
	int k = atoi(argv[1]), r = atoi(argv[3]);
	cout << "Reading edgelist from file " << argv[2] << endl;
	g->readedgelist(argv[2]);
	cout << "Reading edgelist finished!" << endl;
	g->mkGraph();
	cout << "mkGraph finished!" << endl;

	int* col = new int[g->n];
	int tolCol = g->color(col);
	printf("total color = %d\n", tolCol);
	g->coreDecomposition();

	int* corePos = new int[g->n];
	for (int i = 0; i < g->n; i++) corePos[g->coreRank[i]] = i;

	double*** ck = new double** [g->n];

	for (int i = 0; i < g->n; i++)
	{
		ck[i] = new double* [tolCol+1];
		for (int j = 0; j < tolCol + 1; j++)
		{
			ck[i][j] = new double[k + 1];
			ck[i][j][0] = 1;
		}
		for (int j = 0; j < k + 1; j++)
			ck[i][0][j] = 0;

	}

	/*
	//double** ck = new double* [tolCol + 1];

	for (int i = 0; i < tolCol + 1; i++)
	{
		ck[i] = new double[k + 1];
		ck[i][0] = 1;
	}
	for (int i = 1; i < k + 1; i++)
		ck[0][i] = 0;

	*/

	double W = 0, *wpi = new double[g->n];
	int** nbrCol = new int* [g->n];

	//vector<int>** colSet = new vector<int>*[g->n];
	int*** colSet = new int** [g->n], * cPtr = new int[tolCol];
	int* laterNbr = new int[g->n]();

	for (int i = 0; i < g->n; i++)
	{
		nbrCol[i] = new int[tolCol]();

		memset(cPtr, 0, sizeof(int) * tolCol);

		//colSet[i] = new vector<int>[tolCol];
		colSet[i] = new int* [tolCol];
		//for (int j = 0; j < tolCol; j++) colSet[i][j] = NULL;

		for (int j = g->cd[i]; j < g->cd[i + 1]; j++)
		{
			int v = g->adj[j];
			if (corePos[i] < corePos[v])
			{
				laterNbr[i]++;
				nbrCol[i][col[v]]++;

			}
		}

		for (int j = 0; j < tolCol; j++)
			if (nbrCol[i][j]) colSet[i][j] = new int[nbrCol[i][j]];


		for (int j = g->cd[i]; j < g->cd[i + 1]; j++)
		{
			int v = g->adj[j];
			if (corePos[i] < corePos[v])
			{
				colSet[i][col[v]][cPtr[col[v]]++] = v;
			}
		}



		int ptr = 0, cNum = tolCol;

		/*
		while (ptr < tolCol)
		{
			if (nbrCol[ptr])
				nbrCol[cNum++] = nbrCol[ptr];
			ptr++;
		}
		*/

		if (cNum < k - 1)
		{
			wpi[i] = 0;
			continue;
		}

		for (int j = 1; j < cNum + 1; j++)
			for (int p = 1; p < k + 1; p++)
				ck[i][j][p] = ck[i][j - 1][p - 1] * nbrCol[i][j - 1] + ck[i][j - 1][p];

		wpi[i] = ck[i][cNum][k - 1];
		W += wpi[i];
		//printf("i = %d wpi[i] = %f\n", i, wpi[i]);

		//if (wpi[i] < 0)
		//	printf("i = %d wpi[i] = %f\n",i, wpi[i]);


	}

	cout << "W = " << W << endl;

	int degeneracy = 0;
	for (int i = 0; i < g->n; i++) degeneracy = max(degeneracy, laterNbr[i]);

	cout << "degeneracy = " << degeneracy << endl;

	//random_device rd;  //Will be used to obtain a seed for the random number engine
	//mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

	vector<double> wpiv(wpi, wpi + g->n);
	discrete_distribution<size_t> dv{ begin(wpiv), end(wpiv) };

	int* cNbrN = new int[k - 1], * cols = new int[k - 1];
	double cnt = 0;

	double ans = 0, *tmp = new double[1000], tmpu, left;
	for (int run = 0; run < 1; run++)
	{

		for (int i = 0; i < r; i++)
		{
			//int v = SelectPoint(wpi, g->n);



			int v = dv(gen);
			//printf("v = %d\n", v);
			selectNodes(v, tolCol, cNbrN, nbrCol[v], k, g, col, colSet[v], cols, ck[v]);


			//for (int j = 0; j < k - 1; j++) printf("%d ", cNbrN[j]);
			//printf("\n");

			if (isClique(g, cNbrN, k - 1)) cnt += 1;


			//cout << wts[0] << endl;
			//cout << "d(gen) = " << dv(gen) << endl;
		}

		ans += 1.0 * cnt / r * W;
		cnt = 0;
		//tmp[run] = 1.0 * cnt / r * W;
	}

	/*
	double u = ans / 1000, ct = 0;

	for (int i = 0; i < 1000; i++)
	{
		ct += (tmp[i] - u)*(tmp[i]-u);
	}
	*/
	printf("w = %f cnt = %d\n", W, cnt);

	printf("The number of %d-clique is %f\n", k, ans / 1 );

	//cout << "W = " << W << endl;



	cout << "hello world!" << endl;

	return 0;
}

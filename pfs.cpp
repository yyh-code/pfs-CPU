#include <stdio.h>
#include <malloc.h>
#include <windows.h>
#include <math.h>
#include <stdlib.h>

double ofsDist[8];
double xMult = 1.0, yMult = 1.0;
int ni = 0;
struct offsetRec
{
	int ox, oy;
};
struct offsetRec ofs[8] = { 1, 0, -1, 1, 0, -1, 1, 1, -1, 0, 1, -1, 0, 1, -1, -1 };
int nChan, nSink, npq, nVertex;
typedef struct tVertex{
	int ix, iy;
	double zG, zRim, hDist;
	int qi;
	struct tVertex *next;
}pVertex;
pVertex **pq;
pVertex ***adj;
#define onTree 0
double dz = 0.0001;
int M;
int N;
double dx, dy;
double xllcorner;
double yllcorner;
int nodata;
char ncols[15];
char nrows[15];
char xllcorner_label[15];
char yllcorner_label[15];
char cellsize[15];
char NODATA_value[15];
double **z = NULL;
double **zi;
void readzgrid()
{
	int i, j, mul;
	char infile[10];
	FILE *fp;
	printf("输入DEM数据文件名：");
	scanf("%s", infile);
	if ((fp = fopen(infile, "r")) == NULL)
	{
		printf("cannot open file\n");
		return;
	}
	fscanf(fp, "%s %d", &ncols, &M);
	fscanf(fp, "%s %d", &nrows, &N);
	fscanf(fp, "%s %lf", &xllcorner_label, &xllcorner);
	fscanf(fp, "%s %lf", &yllcorner_label, &yllcorner);
	fscanf(fp, "%s %lf", &cellsize, &dx);
	fscanf(fp, "%s %d", &NODATA_value, &nodata);
	z = (double **)malloc(sizeof(double *)* (N+2));//分配指针数组
	z[0] = (double*)malloc((M+2)*(N+2)*sizeof(double));
	for (i = 1; i<=N+1; i++)
	{
		z[i] = z[i - 1] + M+2;
	}
	for (i = 1; i<=N; i++)
	{
		for (j = 1; j<=M; j++)
		{
			fscanf(fp, "%lf", &z[i][j]);
		}
		fscanf(fp, "\n");
	}
	mul = (M+2)*(N+2);
	zi = (double **)malloc(sizeof(double *)* mul);
	pq = (pVertex **)malloc(sizeof(pVertex *)* mul);
	adj = (pVertex ***)malloc(sizeof(pVertex **)* (N+2));
	adj[0] = (pVertex **)malloc((M+2)*(N+2)*sizeof(pVertex *));
	for (i = 1; i<=N+1; i++)
	{
		adj[i] = adj[i - 1] + M+2;//分配每个指针所指向的数组
	}
	fclose(fp);
}

void initHorizOffsets()
{
	int i;
	dy = dx;
	for (i = 0; i <= 7; i++)
	{
		ofsDist[i] = sqrt(pow(ofs[i].ox*dx*xMult, 2) + pow(ofs[i].oy*dy*yMult, 2));
	}
}
void initializeOkPit()
{
	int x, y;
	for (y = 1; y <= N; y++)
	{
		z[y][0] = z[y][1] - 0.5;
		z[y][M + 1] = z[y][M] - 0.5;
	}
	for (x = 0; x <= M + 1; x++)
	{
		z[0][x] = z[1][x] - 0.5;
		z[N + 1][x] = z[N][x] - 0.5;
	}
}
void printpit(int y, int x, FILE *out)
{
	char ch = ' ';
	fputc(ch, out);
	fprintf(out, "%lf", z[y][x]);
}
int pointIsPit(int y, int x, double **z)
{
	int d;
	double z0;
	z0 = z[y][x];

	if (y == 0 || y == N + 1 || x == 0 || x == M + 1 || z0 == -9999)
		return 0;

	for (d = 0; d <= 7; d++)
	{
		if (z[y + ofs[d].oy][x + ofs[d].ox]<z0)
		{
			return 0;
		}
	}
	return 1;
}
bool higherPriority(pVertex *v1, pVertex *v2)
{
	bool result = false;
	if (v1->zRim < v2->zRim)
		result = true;
	else if (v1->zRim == v2->zRim)
	{
		if (v1->zG < v2->zG)
			result = true;
		else if (v1->zG == v2->zG)
		{
			if (v1->hDist < v2->hDist)
				result = true;
		}
	}
	return result;
}
void upHeap(int k)
{
	pVertex *v = NULL;
	v = pq[k];
	while ((k > 1) && higherPriority(v, pq[k / 2]))
	{
		pq[k] = pq[k / 2];
		pq[k]->qi = k;
		k = k / 2;
	}
	pq[k] = v;
	pq[k]->qi = k;
}
void downHeap(int k)
{
	int j = 0;
	pVertex *v = NULL;
	v = pq[k];
	while (k <= npq / 2)
	{
		j = k + k;
		if (j < npq)
		if (higherPriority(pq[j + 1], pq[j]))
			j++;
		if (higherPriority(pq[j], v) == false)
			break;
		pq[k] = pq[j];
		pq[k]->qi = k;
		k = j;
	}
	pq[k] = v;
	pq[k]->qi = k;
}
void PQinsert(pVertex *vOnTree, int x, int y, int d)
{
	pVertex *t;
	nVertex++;
	pq[nVertex] = (pVertex *)malloc(sizeof(pVertex));
	adj[y][x] = pq[nVertex];
	pq[nVertex]->ix = x;
	pq[nVertex]->iy = y;
	pq[nVertex]->zG = z[y][x];
	pq[nVertex]->zRim = vOnTree->zRim;
	if (pq[nVertex]->zG > pq[nVertex]->zRim)
		pq[nVertex]->zRim = pq[nVertex]->zG;
	pq[nVertex]->hDist = vOnTree->hDist + ofsDist[d];
	pq[nVertex]->next = vOnTree;
	npq++;
	t = pq[nVertex];
	pq[nVertex] = pq[npq];
	pq[npq] = t;
	upHeap(npq);
}
pVertex *PQremove()
{
	pVertex *result = NULL;
	pVertex *t;
	result = pq[1];
	pq[1]->qi = onTree;
	t = pq[1];
	pq[1] = pq[npq];
	pq[npq] = t;
	npq--;
	if (npq > 0)
		downHeap(1);
	return result;
}
void pqUpdate(pVertex *vOnTree)
{
	int x = 0, y = 0, d = 0;
	pVertex *v = NULL;
	for (d = 0; d <= 7; d++)
	{
		x = vOnTree->ix + ofs[d].ox;
		y = vOnTree->iy + ofs[d].oy;
		v = adj[y][x];
		if (v == NULL)
		{
			PQinsert(vOnTree, x, y, d);
		}
	}
}
void pqSearch(int y, int x)
{
	pVertex *v = NULL;
	double z0 = 0.0;
	//bool checkOutlet = false;
	while (pointIsPit(y, x, z) == 1)
	{
		//checkOutlet = false;
		z0 = z[y][x];
		pq[1] = (pVertex *)malloc(sizeof(pVertex));
		adj[y][x] = pq[1];
		nVertex = 1;
		pq[1]->ix = x;
		pq[1]->iy = y;
		pq[1]->zG = z0;
		pq[1]->zRim = z0;
		pq[1]->hDist = 0.0;
		pq[1]->qi = onTree;
		pq[1]->next = NULL;
		npq = 0;
		pqUpdate(pq[1]);
		do
		{
			v = PQremove();
			if ((v->zG < z0) || v->iy == 0 || v->iy == N + 1 || v->ix == 0 || v->ix == M + 1)
				break;
			pqUpdate(v);
		} while (npq != 0);
		double slope = 0.0;
	//	nChan++;
		slope = (v->zG - z0) / v->hDist;
		if (slope * ofsDist[0] > -dz)
		{
		//	checkOutlet = true;
			x = v->ix;
			y = v->iy;
			slope = (-dz) / ofsDist[0];
			z[v->iy][v->ix] = z0 + slope * v->hDist;

		}
		do
		{
			v = v->next;
			z[v->iy][v->ix] = z0 + slope * v->hDist;
		} while (v->next != NULL);
		while (nVertex > 0)
		{
			v = pq[nVertex];
			adj[v->iy][v->ix] = NULL;
			free(v);
			nVertex--;
		}
	}
}
void sortDownHeap(double **zi,int k,int high)
{
	int j;
	double * temp = zi[k];
	while (k <= high/2)
	{
		j = k + k;
		if (j<high && *zi[j]>*zi[j + 1])
			j++;
		if (*temp<=*zi[j])
		{
			break;
		}
		else
		{
			zi[k] = zi[j];
			k = j;
		}
	}
	zi[k] = temp;
}
void heapSort(double **zi, int ni)//小根堆
{
	int k;
	double * temp;
	for (k = ni / 2; k >= 1; k--)
		sortDownHeap(zi,k,ni);
	k = ni;
	do{
		temp = zi[1];
		zi[1] = zi[ni];
		zi[ni] = temp;
		ni--;
		sortDownHeap(zi,1,ni);
	} while (ni>1);
	ni = k;
}
void fixinvalidpits(double **zi)
{
	int i = 0, x = 0, y = 0;
	int j;
	heapSort(zi, ni);

	
	for (i = ni; i>=1; i--)
	{
		j = (((int)zi[i]) - ((int)&z[0][0])) / 8;
		x = j % (M+2);
		y = j / (M+2);
		if (pointIsPit(y, x, z) == 1)
		{	
			pqSearch(y, x);
		}
	}
}
void scanGrid()
{
	int x, y, k = 0;
	double *t;
	for (y = 1; y<=N; y++)
	{
		for (x = 1; x<=M; x++)
		{
			if (pointIsPit(y, x, z) == 1)
			{
				ni++;
				zi[ni] = &z[y][x];
			}
		}
	}
	/*	srand(3141621);
	for(int i=1;i<=ni;i++)
	{
	t=zi[i];
	zi[i]=zi[rand() % ni+1];
	zi[rand() % ni+1]=t;
	}*/
}
void getMemory()
{
	int i, x, y;
	int maxVertices = (M+2)*(N+2);
	for (i = 1; i<maxVertices; i++)
	{
		pq[i] = NULL;
	}
	for (y = 0; y<=N+1; y++)
	{
		for (x = 0; x<=M+1; x++)
		{
			adj[y][x] = NULL;
		}
	}
}
void print()
{
	FILE *out;
	int i, j;
	out = fopen("Gridout.txt", "w");
	if (out == NULL)
	{
		printf("无法打开文件\n");
		exit(0);
	}
	fprintf(out, "%s %d\n", "ncols", M);
	fprintf(out, "%s %d\n", "nrows", N);
	fprintf(out, "%s %f\n", "xllcorner", xllcorner);
	fprintf(out, "%s %f\n", "yllcorner", yllcorner);
	fprintf(out, "%s %.10f\n", "cellsize", dx);
	fprintf(out, "%s %d\n", "NODATA_value", nodata);
	for (i = 1; i<=N; i++)
	{
		for (j = 1; j<=M; j++)
		{
			printpit(i, j, out);
		}
		fprintf(out, "\n");
	}
	fclose(out);
}
int main()
{
	int pass = 0;
	readzgrid();
	initHorizOffsets();
	initializeOkPit();
	getMemory();
	do
	{
		pass++;
		nChan = 0;
		nSink = 0;
		ni = 0;
		scanGrid();
		printf("pass:%d 找到%d个洼地\n", pass, ni);
		fixinvalidpits(zi);
	} while (ni != 0);
	free(zi);
	print();
	printf("finished!\n");

	return 0;
}


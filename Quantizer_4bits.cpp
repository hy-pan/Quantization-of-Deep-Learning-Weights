#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#define b1_m  7
#define b1_n  2
#define v1 2
#define N 10

void MatrixOpp(double* A, int m, int n, double* invmat);
void MatrixInver(double* A, int m, int n, double* invmat);
double Surplus(double A[], int m, int n);

void MatrixOpp(double A[], int m, int n, double* invmat)
{
	int i, j, x, y, k;
	double* SP = NULL, * AB = NULL, * B = NULL, X;
	SP = (double*)malloc(m * n * sizeof(double));
	AB = (double*)malloc(m * n * sizeof(double));
	B = (double*)malloc(m * n * sizeof(double));
	X = Surplus(A, m, n);
	X = 1 / X;

	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
		{
			for (k = 0; k < m * n; k++)
				B[k] = A[k];
			{
				for (x = 0; x < n; x++)
					B[i * n + x] = 0;
				for (y = 0; y < m; y++)
					B[m * y + j] = 0;
				B[i * n + j] = 1;
				SP[i * n + j] = Surplus(B, m, n);
				AB[i * n + j] = X * SP[i * n + j];
			}
		}
	MatrixInver(AB, m, n, invmat);
	free(SP);
	free(AB);
	free(B);
}

void MatrixInver(double A[], int m, int n, double* invmat)
{
	int i, j;
	double* B = invmat;

	for (i = 0; i < n; i++)
		for (j = 0; j < m; j++)
			B[i * m + j] = A[j * n + i];
}

double Surplus(double A[], int m, int n) //|A|
{
	int i, j, k, p, r;
	double X, temp = 1, temp1 = 1, s = 0, s1 = 0;
	if (n == 2)
	{
		for (i = 0; i < m; i++)
			for (j = 0; j < n; j++)
				if ((i + j) % 2)
					temp1 *= A[i * n + j];
				else
					temp *= A[i * n + j];
		X = temp - temp1;
	}
	else
	{
		for (k = 0; k < n; k++)
		{
			for (i = 0, j = k; i < m, j < n; i++, j++)
				temp *= A[i * n + j];
			if (m - i)
			{
				for (p = m - i, r = m - 1; p > 0; p--, r--)
					temp *= A[r * n + p - 1];
			}
			s += temp;
			temp = 1;
		}

		for (k = n - 1; k >= 0; k--)
		{
			for (i = 0, j = k; i < m, j >= 0; i++, j--)
				temp1 *= A[i * n + j];
			if (m - i)
			{
				for (p = m - 1, r = i; r < m; p--, r++)
					temp1 *= A[r * n + p];
			}
			s1 += temp1;
			temp1 = 1;
		}
		X = s - s1;
	}
	return X;
}

void swap_float(double* a, double* b)
{
	double t;
	t = *a;
	*a = *b;
	*b = t;
}
void swap_int(int* a, int* b)
{
	int t;
	t = *a;
	*a = *b;
	*b = t;
}
void sort0(double* a)
{
	int j, i;
	double temp;
	int flag;
	for (i = 0; i < 500; i++)
	{
		flag = 0;
		for (j = 499; j > i; --j)
		{
			if (a[j] < a[j - 1])
			{
				flag = 1;
				swap_float(&a[j], &a[j - 1]);
			}
		}
		if (flag == 0)
			break;
	}

}
void equation2(int B[4][500], double w[500], double* v)
{
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 500; j++)
			v[i] += ((double)B[i][j] * w[j]) / 500;
}


int main() {
	double v0[4] = { 0 };
	double v_update[4] = { 0 };
	double weight0[500] = { 0 };//Raw weights
	double weight1[500] = { 0 };//Sorted weights
	double weight2[500] = { 0 };//Quantized weights
	int Biny1[64] = { 0 };
	int Biny2[4][16] = { 0 };
	int Biny3[4][500] = { 0 };
	double eigenvalue[16] = { 0 };
	double temp1[16] = { 0 }; //The result of b times v
	int temp2[16] = { 0 };
	int num[16] = { 0 };
	double quan_err = 0;
	double bou[17] = { 0 };//Boundary

	//load initial v and full precision data and initial B
	errno_t err;
	FILE* fpreader0;
	err = fopen_s(&fpreader0, "E:/研究代码/data1.txt", "r");
	if (err == 1) {
		printf("Error");
		return 0;
	}
	else {
		for (int i = 0; i < 500; i++) {
			fscanf_s(fpreader0, "%lf", &weight0[i]);
			fscanf_s(fpreader0, "%lf", &weight1[i]);
			fscanf_s(fpreader0, "%lf", &weight2[i]);
		}
	}
	fclose(fpreader0);

	FILE* fpreader1;
	err = fopen_s(&fpreader1, "E:/研究代码/Binary_vector4.txt", "r");
	if (err == 1) {
		printf("Error");
		return 0;
	}
	else {
		for (int i = 0; i < 64; i++) {
			fscanf_s(fpreader0, "%d", &Biny1[i]);
		}
	}
	fclose(fpreader1);

	//Sort the raw data
	double* a;
	a = weight1;
	sort0(a);

	bou[0] = weight1[0];
	bou[16] = weight1[499];

	/********************************* Main loop ***************************/
	for (int K = 20; K >= 0; K--) {
		//initialization
		double Biny1_temp[4][16] = { 0 };
		for (int i = 0; i < 16; i++)
		{
			temp1[i] = 0;
			temp2[i] = 0;
		}
		for (int i = 0; i < 4; i++)
			v_update[i] = 0;

		for (int i = 0; i < 16; i++)
			num[i] = 0;
		//Load v and b
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 16; j++) {
				int k = i * 16 + j;
				Biny2[i][j] = Biny1[k];
			}

		FILE* fpreader2;
		err = fopen_s(&fpreader2, "E:/研究代码/data2.txt", "r");
		if (err == 1) {
			printf("Error");
			return 0;
		}
		else {
			for (int i = 0; i < 4; i++) {
				fscanf_s(fpreader2, "%lf", &v0[i]);
			}
		}
		fclose(fpreader2);

		
		for (int i = 0; i < 4; i++)
			printf("v[%d]:%lf\n", i, v0[i]);
       
		//calcute B with initial v
		for (int i = 0; i < 16; i++)
			for (int j = 0; j < 4; j++) {
				temp1[i] += Biny2[j][i] * v0[j];
			}
		for (int i = 0; i < 16; i++)
			temp2[i] = i;
		for (int i = 0; i < 16; i++)
			for (int j = i; j < 16; j++) {
				if (temp1[j] < temp1[i]) {
					swap_int(&temp2[i], &temp2[j]);
					swap_float(&temp1[i], &temp1[j]);
				}
			}
		/*
		for (int i = 0; i < 16; i++)
			printf("temp2[%d]:%d\n", i, temp2[i]);
        */

		/*	for (int i = 0; i < 8; i++)
				for (int j = 0; j < 3; j++) {
					printf("Biny1[%d][%d]=%d\n", j, i, Biny1[j][i]);
				}
		*/
		for (int i = 0; i < 16; i++)
			for (int j = 0; j < 4; j++) {
				Biny1_temp[j][i] = Biny2[j][i];
			}
		for (int i = 0; i < 16; i++)
			for (int j = 0; j < 4; j++) {
				Biny2[j][i] = Biny1_temp[j][temp2[i]];
			}

		//compute bou[i]
		for (int i = 1; i < 16; i++)
			bou[i] = (temp1[i - 1] + temp1[i]) / 2;

		//Extand B
		/*
		for (int j = 0; j < 31; j++)
		{
			for (int i = 0; i < 4; i++)
			{
				for (int n = 0; n < 16; n++)
					Biny3[i][j + 31 * n] = Biny2[i][n];
			}
		}

		for (int j = 496; j < 500; j++)
		{
			for (int i = 0; i < 4; i++)
			{
				Biny3[i][j] = Biny2[i][15];
			}
		}
		*/

		for (int i = 0; i < 500; i++)
			for (int j = 0; j < 16; j++)
				if (weight1[i] > bou[j] && weight1[i] <= bou[j + 1])
				{
					for (int m = 0; m < 4; m++)
						Biny3[m][i] = Biny2[m][j];
					num[j]++;
				}

		for (int i = 0; i < 16; i++)
			printf("num[%d]:%d\n",i, num[i]);

		FILE* fpreader4;
		err = fopen_s(&fpreader4, "E:/研究代码/output_num.txt", "w");
		if (err == 1) {
			printf("Error");
			return 0;
		}
		else {
			for (int i = 0; i < 16; i++) {
				fprintf_s(fpreader4, "%d\n", num[i]);
			}
		}
		fclose(fpreader4);
		/*
		for (int i = 0; i < 4; i++)
		{
			printf("Biny3[%d][0]:%d\n", i, Biny3[i][1]);
			printf("Biny3[%d][31]:%d\n", i, Biny3[i][31]);
			printf("Biny3[%d][62]:%d\n", i, Biny3[i][62]);
			printf("Biny3[%d][93]:%d\n", i, Biny3[i][93]);
			printf("Biny3[%d][124]:%d\n", i, Biny3[i][124]);
		}*/

		/*compute v with b*/
		//compute B*BT
		int temp_B1[4][4] = { 0 };
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				for (int k = 0; k < 500; k++) {
					temp_B1[i][j] += Biny3[i][k] * Biny3[j][k];
				}
		//compute (BBT)-1
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
			{
				printf("%d\t", temp_B1[i][j]);
				if (j == 3)
					printf("\n");
			}

		double temp_B1_double[4][4] = { 0 };
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				temp_B1_double[i][j] = temp_B1[i][j];


		double* arr2 = NULL;
		double* arr = NULL;
		arr = (double*)malloc(sizeof(double) * 4 * 4);
		arr2 = (double*)malloc(sizeof(double) * 4 * 4);
		arr = &temp_B1_double[0][0];

		/*
		for (int i = 0; i < 16; i++)
			printf("arr[%d]:%f\n", i, *(arr + i));
        */

		MatrixOpp(arr, 4, 4, arr2);

		for (int i = 0; i < 16; i++)
		{
			printf("%lf\t", *(arr2 + i));
			if (i % 4 == 3)
				printf("\n");
		}

		//compute update arr2*B
		double arr3[4][500] = { 0 };
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 500; j++)
				for (int k = 0; k < 4; k++)
					arr3[i][j] += (*(arr2 + i * 4 + k)) * Biny3[k][j];
		//compute v
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 500; j++)
				v_update[i] += arr3[i][j] * weight1[j];

		for (int i = 0; i < 4; i++)
			printf("v_update[%d]:%lf\n", i, v_update[i]);

		//write updated v into file
		FILE* fpreader3;
		err = fopen_s(&fpreader3, "E:/研究代码/data2.txt", "w");
		if (err == 1) {
			printf("Error");
			return 0;
		}
		else {
			for (int i = 0; i < 4; i++) {
				fprintf_s(fpreader3, "%lf\n", v_update[i]);
			}
		}
		fclose(fpreader3);	
	}
	/********************************* Main loop ***************************/

	//Compute eigenvalue
	for (int i = 0; i < 16; i++)
	{
		eigenvalue[i] = temp1[i];
		printf("eigenvalue[%d]:%lf\n", i, eigenvalue[i]);
	}

	//output
	FILE* fpreader5;
	err = fopen_s(&fpreader5, "E:/研究代码/output_feature.txt", "w");
	if (err == 1) {
		printf("Error");
		return 0;
	}
	else {
		for (int i = 0; i < 16; i++) {
			fprintf_s(fpreader5, "%lf\n", eigenvalue[i]);
		}
	}
	fclose(fpreader5);

	FILE* fpreader6;
	err = fopen_s(&fpreader6, "E:/研究代码/Edges.txt", "w");
	if (err == 1) {
		printf("Error");
		return 0;
	}
	else {
		for (int i = 0; i < 17; i++) {
			fprintf_s(fpreader6, "%lf\n", bou[i]);
		}
	}
	fclose(fpreader6);

	//verify result
	double temp;
	for (int i = 0; i < 500; i++)
		for (int j = 0; j < 16; j++)
			if (weight1[i] > bou[j] && weight1[i] <= bou[j + 1])
			{
				temp = weight1[i] - bou[j];
				quan_err = temp * temp;
			}

		printf("quantization_err:%lf\n", quan_err);

		//Write quantized weights into file
		/*
		FILE* fpreader4;
		err = fopen_s(&fpreader4, "E:/研究代码/Quantized_data.txt", "w");
		if (err == 1) {
			printf("Error");
			return 0;
		}
		else {
			for (int i = 0; i < 500; i++) {
				fprintf_s(fpreader4, "%lf\n", weight2[i]);
			}
		}
		fclose(fpreader4);
	    */
	return 0;
}

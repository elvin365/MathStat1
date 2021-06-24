#define _CRT_SECURE_NO_WARNINGS
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define N 100
# define M_PI		3.14159265358979323846	/* pi */
#define partition 10
const double x0 = 0.65;

void qsortx(float *a, int low, int  high) {

	int i = 0, j = 0;
	float tmp = 0, medianofthree = 0;

	i = low;
	j = high;
	float k = 0, l = 0, m = 0;
	k = *(a + 0);
	l = *(a + high);
	m = *(a + (high + low) / 2);
	/*Median-of-three*/
	if (k > l && k < m || k < l && k > m)
		medianofthree = k;
	else
		if (l > k && l < m || l < k && l > m)
			medianofthree = l;
		else
			medianofthree = m;
	/*Median-of-three*/
	//printf("\n%d\n", medianofthree);
	do {
		while (*(a + i) < medianofthree)
		{
			i++;
		}
		while (*(a + j) > medianofthree)
		{
			j--;
		}
		if (i <= j)
		{
			if (*(a + i) > *(a + j))
			{
				tmp = *(a + i);
				*(a + i) = *(a + j);
				*(a + j) = tmp;
			}
			i++;
			if (j > 0)
			{
				j--;
			}
		}
	} while (i <= j);

	if (i < high)
	{
		qsortx(a, i, high);
	}
	if (j > low)
	{
		qsortx(a, low, j);
	}

}
void show(float *x)
{
	int i = 0;
	for (i = 0; i < N; i++)
		printf("%f  ", x[i]);
}

void Input(float* x)
{
	int i = 0;
	x[0] = 0.65;
	//i++;
	int t = 6;
	int a = 2;
	int K = 10 * t + a;
	for (i=1; i < N; i++)
	{
		//x[i] = K * x[i - 1]-int (K * x[i - 1]); 	
		x[i] = (int((K * x[i - 1] - int(K * x[i - 1])) * 1000)) / 1000.0; // три занка после запятой

	}
}

/*void Input1(float* Uniq, float *x,float* y)
{
	int counter = 0;
	float temp = x[0];
	int local = 0;
	for (int i = 0; i < N; i++)
	{
		//temp = x[i];
		counter=0;
		while (x[i]==temp)
		{
			counter++;

		}
			
	}
}*/

int LookingForUniq(float* x, float* Uniq)
{
	int i = 0;
	Uniq[0] = x[0];
	int counter = 0;
	for (i = 1; i < N; i++)
	{
		if (x[i] != x[i - 1])
		{
			counter++;
			Uniq[counter] = x[i];
		}
	}
		return counter;

}

void FillinFrequancy(float* x,float* Uniq,float* y)
{
	int k = 0;
	int freq = 0;
	for (int i = 0; i <N;i++ )//идем по х
	{
		if (Uniq[k] == x[i])
		{
			freq++;
			//y[k] = freq / 100.0;
			y[k] = ((double)freq) / N;
		}
		else
		{
			//y[k] = freq / 100;
			k++;
			freq = 1;
			//y[k] = freq / 100.0;
			y[k] = ((double)freq) / N;
		}
	}

}
float VyborSred(float* x)
{
	float sum = 0;
	for (int temp = 0; temp < N; temp++)
	{
		sum = sum + x[temp];
	}
	return ((1.0/N)*sum);
}
float VyborDispersia(float* x,float ChosenX)
{
	float sum = 0;
	for (int temp = 0; temp < N; temp++)
	{
		sum = sum + (x[temp]-ChosenX)*(x[temp] - ChosenX);
	}
	return ((1.0/N)*sum);

}
float FixedDispersia(float* x, float ChosenD)
{
	//ChosenD = ChosenD * ChosenD;
	float final = 0;
	final = (N / (N - 1))*ChosenD;
	//return sqrt(final);
	return final;
}
void whiteinfile(float* Uniq,float* F,int n)
{
	FILE* file;
	file = fopen("C:\\Users\\Elvin\\source\\repos\\MathStat1\\result.txt", "w+t");
	for (int i = 0; i < n+1; i++)
	{

		fprintf(file, "%f ", Uniq[i]);
		fprintf(file, "\n");
	}
	fprintf(file, "\n");
	for (int i = 0; i < n+1; i++)
	{

		fprintf(file, "%f ", F[i]);
		fprintf(file, "\n");
	}
}
void saveinfile(float* x)
{
	FILE* file;
	file = fopen("C:\\Users\\Elvin\\source\\repos\\MathStat1\\1.txt", "w+t");
	for (int i = 0; i < N; i++)
	{

		fprintf(file, "%.3f ", x[i]);
		fprintf(file, "\n");
	}
	fclose(file);

}
void Input2(float* x)
{
	int i = 0;
	x[0] = 0.65;
	//i++;
	
	for (i = 1; i < N; i++)
	{
		x[i] = (int(((11 * x[i - 1]+M_PI) - int(11 * x[i - 1]+ M_PI)) * 1000)) / 1000.0; // три занка после запятой

	}
}
double MakeSumTill(int m,float* Bernulli,int n)
{
	float locsum = 0.0;
	if (m <=n)
	{
		
		for (int i = 0; i <m; i++)
		{
			locsum = locsum + Bernulli[i];
		}
	}
	return locsum;
}
float freqfunc(float x1, float* Uniq, float* y1,int n1)
{
	float temploc = 0;
	int i1 = 0;
	if (x1 <= Uniq[0])
	{
		return 0;
	}
	if (x1 > Uniq[n1])
	{
		for (int i = 0; i < n1; i++)
		{
			temploc = temploc + y1[i];
		}
		return 1;
	}
	if (x1>Uniq[0] && x1<=Uniq[n1])
	{
		while (Uniq[i1]<=x1  && i1<=n1)
		{
			

			temploc = temploc + y1[i1];
			i1++;

		}
		return temploc;
	}

	
}
void saveinfilegist(float* partx4_2,float* party4_2)
{
	FILE* file;
		//printf("[%f  ;  %f] , %f \n", partx4_2[i], partx4_2[i + 1], party4_2[i]);// гистограмма
	file = fopen("C:\\Users\\Elvin\\.spyder\\Gist.txt", "w+t");
	for (int i = 0; i < partition; i++)
	{

		//fprintf(file, "%.3f ", [i]);
		//fprintf(file, "\n");
		fprintf(file,"[%f  ;  %f] , %f \n", partx4_2[i], partx4_2[i + 1], party4_2[i]);// гистограмма
	}
	fclose(file);
}

int main()
{
	float x[N];
	//int y[N];
	//Input(x);
	Input2(x);
	//show(x);

	/*for (int i = 0; i < N - 1; i++)
		y[i] = int(1000*x[i]);
	//qsortx(y, 0, N - 1);*/
	printf("In file viborka\n");
	show(x);
	saveinfile(x);


	qsortx(x, 0, N - 1); //вариационный ряд
	printf("\nvariacionnii rrd\n");
	show(x);


	float y[N];//массив частот
	float Uniq[N];//уникальные
	int n = LookingForUniq(x,Uniq);
	puts("\n");
	//show(Uniq);
	for (int i = 0; i<n + 1; i++)
		printf("%f ", Uniq[i]);

	//Input1(Uniq, x,y);
	FillinFrequancy(x,Uniq,y);
	puts("\n");

	for(int i=0;i<n+1;i++)
	printf("%f ", y[i]);//массив частот

	/*float temp1 = 0;
	for (int i = 0; i < n+1; i++)
		temp1 = temp1 + y[i];
	printf("\n %f \n", temp1);

	for (int i = 0; i<n+1; i++)
		printf("%f ", y[i]);//массив частот


	float F[N];
	F[0] = y[0];//первый скачок
	for (int i = 1; i < n+1; i++)
	{
		F[i] = y[i]+F[i-1];
		
	}
	puts("\n");
	for (int i = 0; i < n+1; i++)
	{
		printf("%f   %f \n", Uniq[i],F[i] );//массив частот
	}

	whiteinfile(Uniq, F,n);*/
	float countofdots = 30.0;
	printf("cumulative function or func of raspredelenie\n");
	for (int i = 0; i < countofdots; i++)
	{

		printf("%f , %f \n", i*(Uniq[n]) / countofdots, freqfunc(i*(Uniq[n]) / countofdots, Uniq, y, n));// первый аргумент - х

	}
	


	float ChosenX=VyborSred(x);
	printf("\n%f is Vyborochnoe Srednie, theoretic 0.5",ChosenX);
	
	float ChosenD = VyborDispersia(x,ChosenX);
	printf("\n%f is Vyborochnoe Dispersia, theortic %f ", ChosenD,1/12.0);

	float FixedD = FixedDispersia(x, ChosenD);
	printf("\n%f is Fixed Dispersia", FixedD);

	float mediana = 0;
	mediana = (x[49] + x[50]) / 2;
	printf("\n%f is mediana, theoretic 0.5\n", mediana);
	//медиана для 
	//ожидаемое



	//--------------------------------------------------
	/*
	2задание
	*/
	printf("2 task");
	float Bernulli[N];
	
	Input(x);
	//Input2(x); //вторая выборка !!!!!!!!!!!!!!!!и для 3
	for (int i = 0; i < N; i++)
	{
		if (x[i] > x0)
			Bernulli[i] = 0;
		else
			Bernulli[i] = 1;
	}

	show(Bernulli);
	
	double Mu[10];
	//printf("\n%f\n", MakeSumTill(20, Bernulli, N));
	for (int j = 0; j < 10; j++)
	{
		Mu[j] = (MakeSumTill((j+1) * 10, Bernulli, N) )/ ((j+1) * 10.0);

	}
	puts("\n");
	for (int i = 0; i < 10; i++)
	{
		printf("%d   %f \n", (i+1)*10, Mu[i]);//от него построить график

	}
	//int amount=MakeSumTill;
	//----------------------------------------------------------------
	//3 задание
	printf("\n3 task\n");
	float Expon1[N];
	float u = 0.418;
	float y1[N]; //массив с частотами


	for (int i = 0; i < N; i++)
	{
			Expon1[i] = (-1/u)*log(x[i]);
			//Expon1[i] = int((Expon1[i] - int(Expon1[i]) )* 1000) / 1000.0; // три занка после запятой
	}

	for (int i = 0; i<N; i++)
		printf("%f ", Expon1[i]);
	saveinfile(Expon1);
	qsortx(Expon1, 0, N - 1); //вариационный ряд

	int n1 = LookingForUniq(Expon1,Uniq);

	FillinFrequancy(Expon1, Uniq, y1);
	//float x1=0;
	//float height =freqfunc(x1,Uniq, y1,n1);
	puts("\n");


	//for (int i = 0; i<n1+1; i++)
    //printf("%f ", Expon1[i]);

	puts("\n");


	//for(int i=0;i<n1+1;i++)
	//printf("%f ", Uniq[i]);
	puts("\n");

	// [0-Uniq[n1]]/20
	float changablepar = 30;
	printf("cumulative function or func of raspredelenie\n");
	for (int i = 0; i < changablepar; i++)
	{
		
		printf("%f , %f \n", i*(Uniq[n1]) / changablepar,freqfunc(i*(Uniq[n1]) / changablepar, Uniq, y1, n1) );// 

	}
	//теоритическая функция 1-е^(-ux)
	
	//гисторграмма 
	
	
	int part = 10;
	float partx[11];//полученный от разбиения отрезка на 10
	float party[10];
	for (int i = 0; i < 11; i++)
	{
		partx[i] =0 +((Uniq[n1] - 0)*i) / 10.0;
	}

	for (int i = 0; i <10; i++)
	{
		party[i] = freqfunc(partx[i+1], Uniq, y1, n1)- freqfunc(partx[i], Uniq, y1, n1);
	}
	puts("\n");
	printf("Gistogramma\n");
	for (int i = 0; i < 10; i++)
	printf("[%f  ;  %f] , %f \n", partx[i], partx[i+1] , party[i]);// гистограмма 
	saveinfilegist(partx, party);


	 ChosenX = VyborSred(Expon1);
	printf("\n%f is Vyborochnoe Srednie, theoretic %f", ChosenX,1.0/u);

	 ChosenD = VyborDispersia(Expon1, ChosenX);
	printf("\n%f is Vyborochnoe Dispersia, theortic %f ", ChosenD, 1.0 / (u*u) );

	 FixedD = FixedDispersia(Expon1, ChosenD);
	printf("\n%f is Fixed Dispersia", FixedD);

	 mediana = 0;
	mediana = (Expon1[49] + Expon1[50]) / 2;
	printf("\n%f is mediana, theoretic %f\n", mediana,log(2)/u);
	


	//для k
	float Expon2[N];
	float k = 3.18;
	float y2[N]; //массив с частотами


	for (int i = 0; i < N; i++)
	{
		Expon2[i] = (-1 / k)*log(x[i]);
		//Expon2[i] = int((Expon2[i] - int(Expon2[i])) * 1000) / 1000.0; // три занка после запятой

	}
	//выборка
	for (int i = 0; i < N; i++)
		printf("%f ",Expon2[i]);
	saveinfile(Expon2);
	qsortx(Expon2, 0, N - 1); //вариационный ряд

	 n1 = LookingForUniq(Expon2, Uniq);

	FillinFrequancy(Expon2, Uniq, y2);
	puts("\n");


	//for (int i = 0; i<n1 + 1; i++)
	//	printf("%f ", Expon2[i]);

	puts("\n");


	//for (int i = 0; i<n1 + 1; i++)
	//	printf("%f ", Uniq[i]);
	puts("\n");

	// [0-Uniq[n1]]/100
	float changable = 30;
	for (int i = 0; i < changable; i++)
	{
		printf("%f , %f \n", i*Uniq[n1] / changable, freqfunc(i*Uniq[n1] / changable, Uniq, y2, n1));// 

	}
	//теоритическая функция 1-е^(-ux)

	//гисторграмма 


	 part = 10;
	 
	for (int i = 0; i < 11; i++)
	{
		partx[i] = 0 + ((Uniq[n1] - 0)*i) / 10.0;
	}

	for (int i = 0; i < 10; i++)
	{
		party[i] = freqfunc(partx[i + 1], Uniq, y2, n1) - freqfunc(partx[i], Uniq, y2, n1);
	}
	puts("\n");
	for (int i = 0; i < 10; i++)
		printf("[%f  ;  %f] , %f \n", partx[i], partx[i + 1], party[i]);// гистограмма 
	saveinfilegist(partx, party);

	//-------------------------------------------------------------------------------




	ChosenX = VyborSred(Expon2);
	printf("\n%f is Vyborochnoe Srednie, theoretic %f", ChosenX, 1.0 / k);

	ChosenD = VyborDispersia(Expon2, ChosenX);
	printf("\n%f is Vyborochnoe Dispersia, theortic %f ", ChosenD, 1.0 / (k*k));

	FixedD = FixedDispersia(Expon2, ChosenD);
	printf("\n%f is Fixed Dispersia", FixedD);

	mediana = 0;
	mediana = (Expon2[49] + Expon2[50]) / 2;
	printf("\n%f is mediana, theoretic %f\n", mediana, log(2) / k);






	//------------------------4
	printf("\n4 task\n");
	const float E = 5.0;
	const float d = 4.18;

	float x4_1[N];//кси
	float x4_2[N];
	Input(x4_1);
	Input2(x4_2);
	//этта
	float Eta1[N];
	float Eta2[N];
	for (int i=0;i<N;i++)
	{
		Eta1[i] = E+cos(2 * M_PI*x4_1[i])*sqrt(-2 * log(x4_2[i]))*sqrt(d);
		Eta2[i]= E+sin(2 * M_PI*x4_1[i])*sqrt(-2 * log(x4_2[i]))*sqrt(d);
		//Eta1[i]=int((Eta1[i] - int(Eta1[i])) * 1000) / 1000.0;
		//Eta2[i] = int((Eta2[i] - int(Eta2[i])) * 1000) / 1000.0;
	}
	printf("Viborki");
	for (int i = 0; i < N; i++)
	{
		printf("%f ", Eta1[i]);
	}
	saveinfile(Eta1);
	float y4[N];
	qsortx(Eta1, 0, N - 1); //вариационный ряд

	 n1 = LookingForUniq(Eta1, Uniq);

	FillinFrequancy(Eta1, Uniq, y4);
	//float x1=0;
	//float height =freqfunc(x1,Uniq, y1,n1);
	puts("\n");


	//for (int i = 0; i<n1 + 1; i++)
	//	printf("%f ", Eta1[i]);

	puts("\n");


	//for (int i = 0; i<n1 + 1; i++)
	//	printf("%f ", Uniq[i]);
	puts("\n");

	// [0-Uniq[n1]]/20
	printf("cumulative function or func of raspredelenie\n");
	changable = 20;
	for (int i = 0; i < changable; i++)
	{

		printf("%f , %f \n", Uniq[0]+i*(Uniq[n1]) / changable, freqfunc(Uniq[0] + i*(Uniq[n1]) / changable, Uniq, y4, n1)); 

	}

	//теоритическая функция 1-е^(-ux)

	//гисторграмма 


	float partx4[partition+1];//полученный от разбиения отрезка на 10
	float party4[partition];
	for (int i = 0; i < partition+1; i++)
	{
		partx4[i] = Uniq[0] + ((Uniq[n1] - Uniq[0])*i) / partition;
	}

	for (int i = 0; i <partition; i++)
	{
		party4[i] = freqfunc(partx4[i + 1], Uniq, y4, n1) - freqfunc(partx4[i], Uniq, y4, n1);
	}
	puts("\n");
	printf("Gistogramma\n");
	for (int i = 0; i < partition; i++)
		printf("[%f  ;  %f] , %f \n", partx4[i], partx4[i + 1], party4[i]);// гистограмма
	saveinfilegist(partx4, party4);






	ChosenX = VyborSred(Eta1);
	printf("\n%f is Vyborochnoe Srednie, theoretic %f", ChosenX,E);

	ChosenD = VyborDispersia(Eta1, ChosenX);
	printf("\n%f is Vyborochnoe Dispersia, theortic %f ", ChosenD, d);

	FixedD = FixedDispersia(Eta1, ChosenD);
	printf("\n%f is Fixed Dispersia", FixedD);

	mediana = 0;
	mediana = (Eta1[49] + Eta1[50]) / 2;
	printf("\n%f is mediana, theoretic %f\n", mediana, E);












	//------------для эта 2
	float y4_2[N];
	for (int i = 0; i<n1 + 1; i++)
		printf("%f ", Eta2[i]);
	saveinfile(Eta2);
	qsortx(Eta2, 0, N - 1); //вариационный ряд

	n1 = LookingForUniq(Eta2, Uniq);

	FillinFrequancy(Eta2, Uniq, y4_2);
	//float x1=0;
	//float height =freqfunc(x1,Uniq, y1,n1);
	puts("\n");


	//for (int i = 0; i<n1 + 1; i++)
	//	printf("%f ", Eta2[i]);

	puts("\n");


	//for (int i = 0; i<n1 + 1; i++)
	//	printf("%f ", Uniq[i]);
	puts("\n");

	// [0-Uniq[n1]]/20
	printf("cumulative function or func of raspredelenie\n");
	changable = 20;
	for (int i = 0; i < changable; i++)
	{

		printf("%f , %f \n", Uniq[0] + i * (Uniq[n1]) / changable, freqfunc(Uniq[0] + i * (Uniq[n1]) / changable, Uniq, y4_2, n1));// 

	}

	//теоритическая функция 1-е^(-ux)

	//гисторграмма 


	float partx4_2[partition + 1];//полученный от разбиения отрезка на 10
	float party4_2[partition];
	for (int i = 0; i < partition + 1; i++)
	{
		partx4_2[i] = Uniq[0] + ((Uniq[n1] - Uniq[0])*i) / partition;
	}

	for (int i = 0; i <partition; i++)
	{
		party4_2[i] = freqfunc(partx4_2[i + 1], Uniq, y4_2, n1) - freqfunc(partx4_2[i], Uniq, y4_2, n1);
	}
	puts("\n");
	printf("Gistogramma\n");
	for (int i = 0; i < partition; i++)
		printf("[%f  ;  %f] , %f \n", partx4_2[i], partx4_2[i + 1], party4_2[i]);// гистограмма
	saveinfilegist(partx4_2, party4_2);







	ChosenX = VyborSred(Eta2);
	printf("\n%f is Vyborochnoe Srednie, theoretic %f", ChosenX, E);

	ChosenD = VyborDispersia(Eta2, ChosenX);
	printf("\n%f is Vyborochnoe Dispersia, theortic %f ", ChosenD, d);

	FixedD = FixedDispersia(Eta2, ChosenD);
	printf("\n%f is Fixed Dispersia", FixedD);

	mediana = 0;
	mediana = (Eta2[49] + Eta2[50]) / 2;
	printf("\n%f is mediana, theoretic %f\n", mediana, E);



	return 0;
}


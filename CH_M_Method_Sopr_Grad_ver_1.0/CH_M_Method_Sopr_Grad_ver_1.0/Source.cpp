// 4M_Puasson.cpp : Defines the entry point for the console application.
//


#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <conio.h>
#include <iomanip>

using namespace std;

double setfunc(double x, double y)
{
	return 1 - x * x - y * y;
}

void proverka(int n, int m)
{

	if ((n == 1) || (n == 2) || (m == 1) || (m == 2))
	{
		cout << "Cетка слишком маленькая, введите ещё раз" << endl;
		cin >> n >> m;
		proverka(n, m);
	}
}


int main()
{
	char c;
	while (1) {

		setlocale(0, "rus");
		int Nmax = 0;
		int S = -1;
		double eps = 0.00000000000001;
		double eps_max = 0;
		double eps_cur = 0;
		double a2, k2, h2;
		int n = 0, m = 0;
		double x = 0, y = 0;
		double **v, **vist;
		double **f;
		//
		double **A_matrix;
		double **B_vector;
		double *B_another_size;
		double *v_another_size;
		int key = 0;
		int temp = 0;
		double tmp = 0;
		int size = 0;
		double *R; // невязка
		double *H; // коэффициент
		double *Help;
		double alfa = 0;
		double beta = 0;
		double c1 = 0;
		double c2 = 0;
		//
		double a = -1, b = 1, c = -1, d = 1;
		double v_old;
		double v_new;
		bool flag = false;

		//double du = -4;

		int text = 0;
		cout << endl;
		cout << "Введите 1 если хотите решить систему 2*2: (3 1; 1 8) = (7 10) с начальным приближением (0 1) \n\nВведите 2, если хоите решить Вариант 4 на произвольной сетке" << endl;
		cin >> text;
		if (text == 1)
		{
			cout << "Введите максимальное число шагов" << endl;
			cin >> Nmax;
			size = 2;

			A_matrix = new double *[2];
			for (int i = 0; i < 2; i++)
				A_matrix[i] = new double[2];

			B_vector = new double *[2];
			for (int i = 0; i < 2; i++)
				B_vector[i] = new double[2];

			B_another_size = new double[2];
			v_another_size = new double[2];
			H = new double[2];
			R = new double[2];
			Help = new double[2];
			
			A_matrix[0][0] = 3; // Зададим матрицу
			A_matrix[0][1] = 1;
			A_matrix[1][0] = 1;
			A_matrix[1][1] = 8;
			B_another_size[0] = 7;
			B_another_size[1] = 10;
			v_another_size[0] = 0;
			v_another_size[1] = 1;
			double aa = 0;
			double aabb = 0;
			double aacc = 0; 
				while (!flag)
				{

					for (int i = 0; i < size; i++)
					{
						R[i] = Help[i] = 0;
					}
					c1 = 0;
					c2 = 0;
					if (S > -1) {
						//посчитаем невязку
						for (int i = 0; i < size; i++)
						{
							for (int j = 0; j < size; j++)
							{

								R[i] += A_matrix[i][j] * v_another_size[j];

							}
							R[i] -= B_another_size[i];
						}
					}
					if (S == 0)
					{
						for (int i = 0; i < size; i++)
						{
							H[i] = -R[i];
						}

						for (int i = 0; i < size; i++)
						{
							for (int j = 0; j < size; j++)
							{
								Help[i] += A_matrix[i][j] * H[j];
							}

						}

						for (int i = 0; i < size; i++)
						{
							c1 += (R[i] * H[i]);
							c2 += (Help[i] * H[i]);
						}
						alfa = -(c1 / c2);
						for (int i = 0; i < size; i++)
						{
							v_another_size[i] = v_another_size[i] + alfa * H[i];
						}
					}
					c1 = 0;
					c2 = 0;
					if (S > 0)
					{
						{
							for (int i = 0; i < size; i++)
							{
								for (int j = 0; j < size; j++)
								{
									Help[i] += A_matrix[i][j] * H[j];
								}
							}
						}

						for (int i = 0; i < size; i++)
						{
							c1 += (Help[i] * R[i]);
							c2 += (Help[i] * H[i]);
						}
						beta = (c1 / c2);

						for (int i = 0; i < size; i++)
						{
							H[i] = H[i] * beta - R[i];
						}

						for (int i = 0; i < size; i++)
						{
							Help[i] = 0;
						}


						for (int i = 0; i < size; i++)
						{
							for (int j = 0; j < size; j++)
							{
								Help[i] += A_matrix[i][j] * H[j];
							}
						}

						c1 = 0;
						c2 = 0;

						for (int i = 0; i < size; i++)
						{
							c1 += (R[i] * H[i]);
							c2 += (Help[i] * H[i]);

						}


						alfa = -(c1 / c2);
						for (int i = 0; i < size; i++)
						{
							v_another_size[i] = v_another_size[i] + alfa * H[i];
						}

					}



					eps_cur = fabs(v_another_size[0] - 2);
					if(eps_cur<fabs(v_another_size[1] - 1))
					eps_cur = fabs(v_another_size[1] - 1);

					

					S = S + 1;

					if ((eps_cur < eps) || (S > Nmax)) { flag = true; }
				}

			cout << endl;

			
			for (int i = 0; i < size; i++)
			{
				cout << endl;
				for (int j = 0; j < size; j++)
				{
					cout << A_matrix[i][j] << "  ";
				}
				cout << " = " << B_another_size[i] << endl;
			}

			for (int i = 0; i < size; i++)
			{
				cout << endl <<"V"<<i<<" = "<< v_another_size[i];
			}

			cout << endl;
			cout << endl;

			cout << "Обеспечивающее точность= " << eps_cur << endl;


			S = 0;

			cout << endl;
		}

		else {
			cout << endl;
			cout << "Введите размерность сетки nxm" << endl;
			cin >> n >> m;
			proverka(n, m);
			cout << "Введите максимальное число шагов" << endl;
			cin >> Nmax;

			
			v = new double *[n + 1];
			for (int i = 0; i < n + 1; i++)
				v[i] = new double[m + 1];

			vist = new double *[n + 1];
			for (int i = 0; i < n + 1; i++)
				vist[i] = new double[m + 1];
			//
			size = (n - 1) * (m - 1);
			A_matrix = new double *[size];
			for (int i = 0; i < size; i++)
				A_matrix[i] = new double[size];

			B_vector = new double *[n + 1];
			for (int i = 0; i < n + 1; i++)
				B_vector[i] = new double[m + 1];

			B_another_size = new double[size];
			v_another_size = new double[size];
			H = new double[size];
			R = new double[size];
			Help = new double[size];
			//

			f = new double *[n + 1];
			for (int i = 0; i < n + 1; i++)
				f[i] = new double[m + 1];


			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
				{
					f[i][j] = 4;
					v[i][j] = 0;
					vist[i][j] = 0;
					B_vector[i][j] = 0;
				}
			}

			for (int i = 0; i < size; i++)
			{
				v_another_size[i] = 0;
				B_another_size[i] = 0;
				H[i] = 0;
				R[i] = 0;
				Help[i] = 0;
				for (int j = 0; j < size; j++)
				{
					A_matrix[i][j] = 0;
				}
			}

			double h = ((b - a) / n);
			double k = ((d - c) / m);

			for (int i = 0; i < n + 1; i++)
			{
				x = a + i * h;
				v[i][0] = setfunc(x, c);
				v[i][m] = setfunc(x, d);
			}

			for (int j = 0; j < m + 1; j++)
			{
				y = c + j * k;
				v[0][j] = setfunc(a, y);
				v[n][j] = setfunc(b, y);
			}

			

			h2 = (1 / (h*h));
			k2 = (1 / (k*k));
			a2 = -2 * (h2 + k2);
			key = n - 1;

			for (int i = 0; i < size; i++)
			{
				if (i%key == 0 && i > 0)
					temp++;
				for (int j = 0; j < size; j++)
				{

					if (i == j)
					{
						A_matrix[i][j] = a2;
						if ((j + 1) < (temp + 1)*key)
							A_matrix[i][j + 1] = h2;
						if ((j + key) < size)
							A_matrix[i][j + key] = k2;
						if ((j - 1) > -1 + temp * key)
							A_matrix[i][j - 1] = h2;
						if ((j - key) > -1)
							A_matrix[i][j - key] = k2;
					}
				}
			}

			for (int j = 1; j < m; j++)
			{
				for (int i = 1; i < n; i++)
				{
					B_vector[i][j] += -f[i][j];
					if (i == 1)
					{
						B_vector[i][j] -= v[0][j] * k2;
					}
					if (i == n - 1)
					{
						B_vector[i][j] -= v[n][j] * k2;
						tmp = B_vector[i][j];
					}
					if (j == 1)
					{
						B_vector[i][j] -= v[i][0] * h2;
					}
					if (j == m - 1)
					{
						B_vector[i][j] -= v[i][m] * h2;
					}

				}
			}



			for (int j = 1; j < m; j++)
			{
				for (int i = 1; i < n; i++)
				{
					B_another_size[(i - 1)*(n - 1) + j - 1] = B_vector[i][j];
					v_another_size[(i - 1)*(n - 1) + j - 1] = v[i][j];
				}
			}

			
			while (!flag)
			{
				for (int i = 0; i < size; i++)
				{
					R[i] = Help[i] = 0;
				}
				c1 = 0;
				c2 = 0;
				if (S > -1) {
					//посчитаем невязку
					for (int i = 0; i < size; i++)
					{
						for (int j = 0; j < size; j++)
						{

							R[i] += A_matrix[i][j] * v_another_size[j];

						}
						R[i] -= B_another_size[i];
					}
				}
				if (S == 0)
				{
					for (int i = 0; i < size; i++)
					{
						H[i] = -R[i];
					}

					for (int i = 0; i < size; i++)
					{
						for (int j = 0; j < size; j++)
						{
							Help[i] += A_matrix[i][j] * H[j];
						}

					}

					for (int i = 0; i < size; i++)
					{
						c1 += (R[i] * H[i]);
						c2 += (Help[i] * H[i]);
					}
					alfa = -(c1 / c2);
					for (int i = 0; i < size; i++)
					{
						v_another_size[i] = v_another_size[i] + alfa * H[i];
					}
				}
				c1 = 0;
				c2 = 0;
				if (S > 0)
				{
					{
						for (int i = 0; i < size; i++)
						{
							for (int j = 0; j < size; j++)
							{
								Help[i] += A_matrix[i][j] * H[j];
							}
						}
					}

					for (int i = 0; i < size; i++)
					{
						c1 += (Help[i] * R[i]);
						c2 += (Help[i] * H[i]);
					}
					if (c2 != 0)
					beta = (c1 / c2);

					for (int i = 0; i < size; i++)
					{
						H[i] = H[i] * beta - R[i];
					}

					for (int i = 0; i < size; i++)
					{
						Help[i] = 0;
					}


					for (int i = 0; i < size; i++)
					{
						for (int j = 0; j < size; j++)
						{
							Help[i] += A_matrix[i][j] * H[j];

						}
					}

					c1 = 0;
					c2 = 0;

					for (int i = 0; i < size; i++)
					{
						c1 += (R[i] * H[i]);
						c2 += (Help[i] * H[i]);

					}

					if(c2!=0)
					alfa = -(c1 / c2);
					for (int i = 0; i < size; i++)
					{
						v_another_size[i] = v_another_size[i] + alfa * H[i];
					}

				}

				for (int i = 0; i < size; i++)
				{
					eps_cur = alfa * H[i];
				}
				
				S = S + 1;

				if (((eps_cur)>0 && (eps_cur<eps)) || (S >= Nmax)) { flag = true; }
			}

			cout << endl;

			cout << "При решении разностной схемы с помощью метода Сопряженных градиентов с параметрами:" << endl;
			cout << "u = 1 - x^2 - y^2" << endl;
			cout << "a = " << a << ". b = " << b << endl;
			cout << "c = " << c << ". d = " << d << endl;
			cout << "Nmax = " << Nmax << ". Eps = " << eps << endl;
			cout << "h = " << h << ". k = " << k << endl;
			cout << "За число итераций S = " << S << endl;
			cout << "Получено решение " << endl;
			//Можно посмотреть разностную схему ввиде СЛАУ
		/*	for (int i = 0; i < size; i++)  
			{
				cout << endl;
				for (int j = 0; j < size; j++)
				{
					cout << A_matrix[i][j] << "  ";
				}
				cout << " = " << B_another_size[i];
			}*/
			cout << endl;
			cout << endl;

			for (int i = size-1; i >-1; i--)
			{
				if((i+1)%key==0)
					cout << endl;
				cout << "          " << v_another_size[i] ;
			}

			cout << endl;
			cout << "i=0.." << n << ". j=0.." << m << endl;
			cout << "Обеспечивающее точность= " << eps_cur << endl;

			cout << "Сравнение истинного и численного решений" << endl;


			for (int i = 0; i < n + 1; i++)
			{
				x = a + i * h;
				for (int j = 0; j < m + 1; j++)
				{
					y = c + j * k;
					vist[i][j] = setfunc(x, y);
				}

			}

			for (int j = m-1; j > 0; j--)
			{
				for (int i = n-1; i > 0; i--)
				{
					cout << "          " << vist[j][i];
				}
				cout << endl;
			}

			S = 0;

			cout << endl;
			cout << "Нажмите ESC - выхода" << endl << "Остальные клавиши - продолжить" << endl;
			c = _getch();

			if (c == 27) exit(0);


		}
	}

	system("pause");
	return 0;
}



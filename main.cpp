#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <random>
#include <chrono>
#include <cmath>
#include <fstream>

using namespace std;
using namespace chrono;

typedef vector<vector<double>> matrix;
typedef vector<double> vec;


vec MulMV(const matrix &m, const vec &x)
{
    vec res = vec(m.size(), 0.);
    for (size_t i = 0; i < m.size(); i++)
    {
        for (size_t j = 0; j < m.size(); j++)
        {
            res[i] += m[i][j] * x[j];
        }
    }
    return res;
}

vec operator-(const vec &x, const vec &y)
{
    vec ans(x.size(), 0);
    for (size_t i = 0; i < x.size(); i++)
    {
        ans[i] = x[i] - y[i];
    }
    return ans;
}

vec operator+(const vec &x, const vec &y)
{
    vec ans(x.size(), 0);
    for (size_t i = 0; i < x.size(); i++)
    {
        ans[i] = x[i] + y[i];
    }
    return ans;
}

vec operator*(const vec &x, double &y)
{
    vec ans(x.size(), 0);
    for (size_t i = 0; i < x.size(); i++)
    {
        ans[i] = x[i] * y;
    }
    return ans;
}

struct Equasions
{
    long N = 0; //Размерность матрицы
    matrix Matrix; //Матрица
    vec x; //Правильный ответ
    vec res; //Полученный ответ
    vec f; //Соответствущая x правая часть
    double h1, h2;//ограничение на спектор матрицы
    vector<uint64_t> it_par = {1};//вектор итерационных параметров
    vec error;//вектор среднеквадратической нормы погрешности

    Equasions(string filename)
    {

        string line;
        ifstream in(filename);
        if (in.is_open())
        {
            while (getline(in, line))
            {
                if (line.empty())
                    break;
                Matrix.emplace_back();
                for (size_t i = 0; i < line.length(); i++)
                {
                    string x = "";
                    while (line[i] != ',' && i != line.length())
                    {
                        x += line[i];
                        ++i;
                    }
                    Matrix[N].push_back(stod(x, nullptr));
                }
                ++N;
            }
        }
        in.close();  
    }

    void CalculateRandX(const unsigned int seed)
    {
        x = vec(N, 0);

        default_random_engine rand(seed);
        uniform_real_distribution<double> dist(0.0, 1.0);

        for (size_t i = 0; i < N; i++)
        {
            x[i] = dist(rand);
        }
    }

    void CalculateF()
    {
        f = MulMV(Matrix, x);
    }

    void CalculateDiag(size_t i)
    {
        double sum = 0;
        for (size_t j = 0; j < i && i > 0; j++) {
            sum += Matrix[j][i] * Matrix[j][i];
        }

        Matrix[i][i] =  sqrt(Matrix[i][i] - sum);
    }

    //метод Холецкого
    void CalculateNonDiag(size_t k, size_t i) 
    {
        double sum = 0;
        for (size_t j = 0; j < k && k > 0; j++) {
            sum += Matrix[j][k] * Matrix[j][i];
        }
        Matrix[k][i] = (Matrix[k][i] - sum) / Matrix[k][k];
    }

    void ProcessRow(size_t k) 
    {
        for (size_t i = 0; i < k; i++)
            Matrix[k][i] = 0.;
        
        CalculateDiag(k);

        for (size_t i = k + 1; i < N; i++)
            CalculateNonDiag(k, i);
    }

    void CalculateL() 
    {
        for (size_t i = 0; i < N; i++)
        {
            ProcessRow(i);
        }
    }

    void SecondGaussBackward()
    {
        for (size_t i = N; i > 0; i--)
        {
            double d = 0;
            for (size_t j = i + 1; j <= N; j++)
            {
                d += Matrix[i - 1][j - 1] * res[j - 1];
            }
            res[i - 1] = (res[i - 1] - d) / Matrix[i - 1][i - 1];
        }
    }

    void FirstGaussBackward()
    {
        res = vec(N, 0.);
        for (size_t i = 0; i < N; i++)
        {
            double d = 0;
            for (size_t j = 0; j < i; j++)
            {
                d += Matrix[j][i] * res[j];
            }
            res[i] = (f[i] - d) / Matrix[i][i];
        }
    }

    double MaxNorm()
    {
        double ans = 0;
        for (int i = 0; i < N; ++i) 
        {
            ans = max(ans, fabs(x[i] - res[i]));
        }
        return ans;
    }

    // Метод Чебышева
    double RMSE()
    {
        double ans = 0;
        for (int i = 0; i < N; ++i) 
        {
            ans += (res[i] - x[i]) * (res[i] - x[i]);
        }
        return sqrt(ans / N);
    }

    double Relative_Error()
    {
        double ans = 0;
        double nx = 0;
        for (int i = 0; i < N; ++i) 
        {
            ans += (res[i] - x[i]) * (res[i] - x[i]);
            nx += x[i] * x[i];
        }
        return sqrt(ans) / sqrt(nx);
    }

    void Interval_Characteristics()
    {
        double sum = 0;
        double min = numeric_limits<double>::max();
        double max = 0;
        for (size_t i = 0; i < N; i++)
        {
            double center = fabs(Matrix[i][i]);
            sum = 0;
            for (size_t j = 0; j < N; j++)
            {
                sum += fabs(Matrix[i][j]);
            }
            sum -= center;
            if (center - sum < min)
            {
                min = center - sum ;
            }
            if (center + sum > max)
            {
                max = center + sum ;
            }
        }
        h1 = min;
        h2 = max;
    }

    void Iteration_Parameters()
    {
        vector<uint64_t> new_par;
        uint64_t m = it_par.size();
        for (size_t i = 0; i < it_par.size(); i++)
        {
            new_par.push_back(it_par[i]);
            uint64_t new_p = 4 * m - it_par[i];
            new_par.push_back(new_p);
        }
        it_par = new_par;
    }
    
    void Chebyshev_method()
    {
        double tau0 = 2 / (h1 + h2);
        double r0 = (h2 - h1) / (h2 + h1);
        res = vec(N, 0);
        uint64_t m = it_par.size();
        error = vec(m, 0);
        for (size_t i = 0; i < m; i++)
        {
            double tauk = tau0 / (1 - r0 * cos(M_PI * it_par[i]/ (2 * m)));
            res = (f - MulMV(Matrix, res)) * tauk + res;
            error[i] = RMSE();
        }
    }
};



bool CheckMatrix(const Equasions m)
{
    matrix mt = m.Matrix;
    bool res = true;
    for (size_t i = 0; i < m.N; i++)
    {
        for (size_t j = 0; j < m.N; j++)
        {
            if (mt[j][i] != mt[i][j])
            {
                res = false;
                break;
            };
        }
    }
    return res;
}

int main()
{
    Equasions m = Equasions("SLAU_var_2.csv"); //Считывает матрицу из файла
    for (size_t i = 0; i < m.N; i++)//Аx+x=f <=> (A+I)x=f
    {
        m.Matrix[i][i] += 1;
    }
    
    if (!CheckMatrix(m)) //Проверка на симметричность
    {
        cout << "Матрица не является симметричной" << endl;
        return 0;
    }

    Equasions copy_m = m;// Сохраняет копию матрицы
    double norm = 0; //Переменная для подсчета средней среднеквадратической нормы погрешности
    for (size_t i = 0; i < 100; i++) //Производим 100 итераций вычисления решения системы уравнений
    {
        m.CalculateRandX(i);
        m.CalculateF();m.CalculateL();

        m.FirstGaussBackward();
        m.SecondGaussBackward();
        norm += m.RMSE(); //Считаем среднеквадратическую погрешность
        m = copy_m; //возвращаем матрице ее первоначальные значения
    }
    
    double Kholetsky_error = norm / 100;
    cout << "Среднеквадратическую погрешность методом Холецкого = " << Kholetsky_error << endl;
    
    m.Interval_Characteristics();
    cout << "Оценки спектра: h1 = " << m.h1 << " h2 = " << m.h2 << endl;

    copy_m = m;// Сохраняет копию матрицы
    double Chebyshev_error = numeric_limits<double>::max();
    int n = -1;//степень двойки количества итераций
    while(Chebyshev_error > Kholetsky_error)
    {
        ++n;
        Chebyshev_error = 0;
        for (size_t i = 0; i < 10; i++) 
        {
            m.CalculateRandX(i);
            m.CalculateF();
            m.Chebyshev_method();
            Chebyshev_error += m.RMSE();
        }
        Chebyshev_error /= 10;
        cout << "Среднеквадратической отклонение при числе итераций " << (1 << n) << " = "<< Chebyshev_error << endl;
        m.Iteration_Parameters();
    }
    
    cout << "Наименьший показатель степени двойки = " << n << endl;
    cout << "Относительная погрешность решения, полученного методом Чебышева = " << m.Relative_Error() << endl;
    ofstream file;
    file.open("Errors.txt");
    file << '[';
    for (auto elem : m.error) {
        file << std::fixed << std::setprecision(20) << elem << ", ";
    }
    file << ']';
    file.close();
    return 0;
}
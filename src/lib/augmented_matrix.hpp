#include <vector>
#include <iostream>
#include <iomanip>
#include "prompt.hpp"
#include <random>
#include <fstream>
#include <ctime>
#include <string>

using namespace std;

const double EPS = 1E-9;

enum am_coeff_type
{
    variable_coefficient,
    free_member,
};

enum am_print_mode
{
    current,
    original,
};

// singleton
class augmented_matrix
{
private:
    int m;
    int n;
    int rank;
    int augmented_rank;
    vector<vector<double>> matrix;
    vector<vector<double>> current_matrix;
    vector<double> roots;
    ofstream fout;

public:
    augmented_matrix()
    {
    }

    augmented_matrix(int m, int n)
    {
        this->m = m;

        this->n = n + 1;

        this->roots.resize(m);

        this->matrix.resize(m);
        this->current_matrix.resize(m);

        for (int i = 0; i < m; i++)
        {
            this->current_matrix.at(i).resize(this->n);
            this->matrix.at(i).resize(this->n);
        }
    }

    augmented_matrix(const char *filepath)
    {
        try
        {
            ifstream fin(filepath);

            if (fin.is_open())
            {
                char symbol;

                string word;

                am_coeff_type coeff_type = am_coeff_type::variable_coefficient;

                bool at_spaces = false;
                bool free_member_read = false;
                am_coeff_type mode = am_coeff_type::variable_coefficient;
                int m = 0, n = 0, columns = 0;

                while (fin.get(symbol))
                {
                    if (symbol == ' ')
                    {
                        if (!at_spaces && mode == am_coeff_type::variable_coefficient)
                        {
                            n += 1;
                        }

                        at_spaces = true;

                        word.clear();
                    }
                    else if (symbol == '\n')
                    {
                        mode = am_coeff_type::variable_coefficient;

                        if (columns == 0)
                        {
                            columns = n;
                        }

                        if (columns != 0 && n != columns)
                        {
                            throw runtime_error("Уравнения из файла имеют разную размерность");
                        }

                        m += 1;

                        n = 0;

                        free_member_read = false;
                    }
                    else if (symbol == '|')
                    {
                        mode = am_coeff_type::free_member;
                    }
                    else
                    {
                        at_spaces = false;

                        word += symbol;
                    }

                    if (word.length() > 0)
                    {
                        try
                        {
                            stof(word);
                        }
                        catch (exception &e)
                        {
                            throw runtime_error("Встречено не число");
                        }

                        if (mode == am_coeff_type::free_member)
                        {
                            free_member_read = true;
                        }

                        if (free_member_read)
                        {
                            throw runtime_error("Встречено больше одного свободного члена в уравнении");
                        }

                        mout << "\t(" << word << ")"
                             << "\tn = " << n << "\tm = " << m << endl;
                        ;
                    }
                }

                if (fin.eof())
                {
                    m += 1;

                    mout << "\t(" << word << ")"
                         << "\tn = " << n << "\tm = " << m << endl;
                    ;
                }

                if (m != n)
                {
                    throw runtime_error("СЛАУ с количеством неизвестных отличных от количества уравнений не поддерживаются");
                }

                mout << "Предоставленный файл матрицы соответствует шаблону" << endl
                     << endl;

                // validation loop
                // while (fin >> word)
                // {
                //     if (word == "|")
                //     {
                //         coeff_type = am_coeff_type::variable_coefficient;
                //     }
                //     else
                //     {
                //         try
                //         {
                //             if(coeff_type = am_coeff_type::)
                //             m += 1;
                //             stof(word);
                //         }
                //         catch (exception &e)
                //         {
                //             throw runtime_error("При чтении коэффициентов СЛАУ из файла встречено не число");
                //         }
                //     }

                //     mout << line << endl;
                // }

                // word.clear();

                // fin.close();

                // fin.open(filepath);

                // // writing loop
                // while (getline(fin, line))
                // {
                //     mout << line << endl;
                // }
            }
            else
            {
                throw runtime_error("Не удалось открыть файл для чтения");
            }
        }
        catch (exception &e)
        {
            cout << "[Ошибка]:\t" << e.what() << endl;
        }
    }

    void random_filling()
    {
        random_device rd;
        mt19937 rng{rd()};

        uniform_int_distribution<int> uid(1, 64);

        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                this->current_matrix.at(i).at(j) = this->matrix.at(i).at(j) = uid(rng);
            }
        }
    }

    void scan()
    {
        for (int i = 0; i < this->m; i++)
        {
            for (int j = 0; j < this->n; j++)
            {
                char *prompt_text = new char[256];

                if (j == m)
                {
                    sprintf(prompt_text, "Введите значение уравнения для уравнения №%d: ", i + 1);

                    int value = prompt_double((const char *)prompt_text);

                    this->matrix.at(i).at(j) = value;
                }
                else
                {
                    sprintf(prompt_text, "Введите значение A[%d, %d]: ", i + 1, j + 1);

                    int value = prompt_double((const char *)prompt_text);

                    this->matrix.at(i).at(j) = value;
                }

                delete[] prompt_text;
            }
        }

        this->current_matrix = this->matrix;
    }

    void calculate_ranks()
    {
        this->rank = this->get_rank(false);
        this->augmented_rank = this->get_rank(true);
    }

    void print(am_print_mode mode)
    {
        vector<vector<double>> *source = nullptr;

        if (mode == am_print_mode::original)
        {
            source = &this->matrix;
        }

        if (mode == am_print_mode::current)
        {
            source = &this->current_matrix;
        }

        for (int i = 0; i < this->m; i++)
        {
            for (int j = 0; j < this->n; j++)
            {
                if (j == m)
                {
                    mout << setw(6) << "| " << left << setw(12) << (*source).at(i).at(j);
                }
                else
                {
                    mout << setw(12) << right << (*source).at(i).at(j);
                }
            }

            mout << endl;
        }

        mout << endl;
    }

    void print_roots()
    {
        mout << "Корни СЛАУ:" << endl
             << endl;

        for (int i = 0; i < this->m; i++)
        {
            mout << "X[" << i + 1 << "] = " << roots.at(i);
            mout << endl;
        }

        mout << endl;
    }

    void set_value_at_matrix(int m, int n, double value)
    {
        this->matrix.at(m).at(n) = value;
    }

    augmented_matrix(vector<vector<double>> slae)
    {
        const int m = slae.size(), n = slae.at(0).size();

        this->matrix.resize(m);

        for (int i = 0; i < m; i++)
        {
            this->matrix.at(i).resize(n + 1);
        }

        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n + 1; j++)
            {
                matrix.at(i).at(j) = slae.at(i).at(j);
            }
        }
    }

    vector<vector<double>> clone_matrix()
    {
        return this->matrix;
    }

    vector<vector<double>> get_triangle_form()
    {
        vector<vector<double>> t = this->clone_matrix();

        int _n = t.size();
        int _m = t[0].size();

        int rank = 0;

        vector<bool> row_selected(_n, false);

        for (int i = 0; i < _m; ++i)
        {
            int j;
            for (j = 0; j < _n; ++j)
            {
                if (!row_selected[j] && abs(t[j][i]) > EPS)
                    break;
            }

            if (j != _n)
            {

                row_selected[j] = true;

                for (int p = i + 1; p < _m; ++p)
                {
                    t[j][p] /= t[j][i];
                }

                for (int k = 0; k < _n; ++k)
                {
                    if (k != j && abs(t[k][i]) > EPS)
                    {
                        for (int p = i + 1; p < _m; ++p)
                            t[k][p] -= t[j][p] * t[k][i];
                    }
                }
            }
        }

        return t;
    }

    int get_rank(bool for_augmented)
    {
        vector<vector<double>> A = this->clone_matrix();

        if (for_augmented)
        {
            mout << "Получение ранга для дополненной матрицы A|B: ";

            for (int i = 0; i < A.size(); i++)
            {
                A[i] = vector<double>(A[i].begin(), A[i].end() - 1);
            }
        }
        else
        {
            mout << "Получение ранга для матрицы A: ";
        }

        int _n = A.size();
        int _m = A[0].size();

        int rank = 0;

        vector<bool> row_selected(_n, false);

        for (int i = 0; i < _m; ++i)
        {
            int j;
            for (j = 0; j < _n; ++j)
            {
                if (!row_selected[j] && abs(A[j][i]) > EPS)
                    break;
            }

            if (j != _n)
            {
                ++rank;

                row_selected[j] = true;

                for (int p = i + 1; p < _m; ++p)
                {
                    A[j][p] /= A[j][i];
                }

                for (int k = 0; k < _n; ++k)
                {
                    if (k != j && abs(A[k][i]) > EPS)
                    {
                        for (int p = i + 1; p < _m; ++p)
                            A[k][p] -= A[j][p] * A[k][i];
                    }
                }
            }
        }

        mout << rank << endl
             << endl;

        return rank;
    }

    void solve()
    {
        try
        {
            mout << "Ранг матрицы A = " << this->rank << ", ранг расширенной матрицы A|B = " << this->augmented_rank << endl
                 << endl;

            if (this->rank != this->augmented_rank)
            {
                throw "Ранги не равны, решения нет";
            }

            if (this->rank == this->augmented_rank)
            {
                mout << "Ранги равны, решение одно:" << endl
                     << endl;

                double max = 0;
                int k = 0, index = 0;

                vector<vector<double>> *a = &this->current_matrix;
                vector<double> y = {};

                y.resize(this->m);

                vector<double> *x = &this->roots;

                int n = this->m;

                for (int i = 0; i < n; i++)
                {
                    y[i] = (*a)[i][n];
                }

                while (k < n)
                {
                    mout << "[Действ.]\tПроход №" << k + 1 << endl
                         << endl;

                    max = abs((*a)[k][k]);

                    mout << "[Действ.]\tЭлемент для элиминации "
                         << "A[" << k + 1 << "][" << k + 1 << "] = " << max << endl
                         << endl;

                    index = k;

                    for (int i = k + 1; i < n; i++)
                    {
                        if (abs((*a)[i][k]) > max)
                        {
                            max = abs((*a)[i][k]);
                            index = i;
                        }
                    }

                    if (max < EPS)
                    {
                        mout << "Решение получить невозможно" << endl
                             << endl;
                        return;
                    }

                    for (int j = 0; j < n; j++)
                    {
                        double temp = (*a)[k][j];
                        (*a)[k][j] = (*a)[index][j];
                        (*a)[index][j] = temp;
                    }

                    double temp = y[k];
                    y[k] = y[index];
                    y[index] = temp;

                    if (k != index)
                    {
                        mout << "[Действ.]\tПерестановка строк " << k + 1 << " и " << index + 1 << endl
                             << endl;
                    }

                    for (int i = k; i < n; i++)
                    {
                        double temp = (*a)[i][k];

                        if (abs(temp) < EPS)
                        {
                            continue; // для нулевого коэффициента пропустить
                        }

                        for (int j = k; j < n; j++)
                        {
                            (*a)[i][j] = (*a)[i][j] / temp;
                        }

                        this->current_matrix[i][this->m] = y[i] = y[i] / temp;

                        if (temp != 1)
                        {
                            mout << "[Действ.]\tДеление строки №" << k + 1 << " на " << temp << endl
                                 << endl;

                            this->print(am_print_mode::current);
                        }

                        if (i == k)
                        {
                            continue; // уравнение не вычитать само из себя
                        }

                        mout << "[Действ.]\tВычитание строки №" << i + 1 << " из строки № " << k + 1 << endl
                             << endl;

                        for (int j = 0; j < n; j++)
                        {
                            mout << "[Действ.]\tA[" << i + 1 << "][" << j + 1 << "] (" << (*a)[i][j] << ") - A[" << k + 1 << "][" << j + 1 << "] (" << (*a)[k][j] << ") = " << (*a)[i][j] - (*a)[k][j] << endl
                                 << endl;

                            (*a)[i][j] = (*a)[i][j] - (*a)[k][j];

                            this->print(am_print_mode::current);
                        }

                        this->current_matrix[i][this->m] = y[i] = y[i] - y[k];

                        mout << "[Действ.]\tA[" << i + 1 << "][" << (this->m + 1) << "] (" << (*a)[i][this->m] << ") - A[" << k + 1 << "][" << (this->m + 1) << "] (" << (*a)[k][this->m] << ") = " << (*a)[i][this->m] - (*a)[k][this->m] << endl
                             << endl;
                    }

                    this->print(am_print_mode::current);

                    k++;
                }

                mout << "[Действ.]\tОбратный ход" << endl
                     << endl;

                for (k = n - 1; k >= 0; k--)
                {
                    mout << "[Действ.]\tX[" << k + 1 << "] = Y[" << k << "] = " << y[k + 1] << endl
                         << endl;

                    (*x)[k] = y[k];

                    for (int i = 0; i < k; i++)
                        y[i] = y[i] - (*a)[i][k] * (*x)[k];
                }

                return;
            }

            if (this->rank < this->augmented_rank)
            {
                throw "Множество решений";
                // a lot of solutions
            }
        }
        catch (const char *error)
        {
            mout << error << endl
                 << endl;
        }
    }

    ~augmented_matrix()
    {
        mout << "[Сист.]: Освобождение памяти" << endl;
    }
};
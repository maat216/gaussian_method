#include <iostream>
#include <string>
#include "lib/augmented_matrix.hpp"
#include <random>

using namespace std;

static random_device rd;
static mt19937 rng{rd()};

enum operating_mode
{
    manual_mode = 1,
    random_mode,
    file_mode
};

int main()
{
    operating_mode mode = (operating_mode)prompt_int("Выберите режим работы приложения: \n\n(1) Ручной режим\n(2) Случайное задание СЛАУ\n(3) Чтение из файла input.txt: ", 1, 3);

    if (mode == operating_mode::manual_mode)
    {
        mout << "Выбран ручной режим:" << endl
             << endl;

        int m1 = prompt_int("Введите прядок СЛАУ (размерность квадратной матрицы А) (минимум 2): ", 2);

        augmented_matrix *am1 = new augmented_matrix(m1, m1);

        am1->scan();

        am1->calculate_ranks();

        am1->print(am_print_mode::current);

        am1->solve();

        am1->print_roots();

        delete am1;
    }
    else if (mode == operating_mode::random_mode)
    {
        mout << "Выбран режим случайного задания СЛАУ:" << endl
             << endl;

        uniform_int_distribution<int> uid(2, 8);

        int m2 = uid(rng);

        mout << "Случайная размерность M = " << m2 << endl
             << endl;

        augmented_matrix *am2 = new augmented_matrix(m2, m2);

        am2->random_filling();

        am2->calculate_ranks();

        am2->print(am_print_mode::current);

        am2->solve();

        am2->print_roots();
    }
    else if (mode == operating_mode::file_mode)
    {
        mout << "Выбран режим чтения из файла input.txt:" << endl
             << endl;

        augmented_matrix *am3 = new augmented_matrix("input.txt");
    }

    cin.get();

    return 0;
}
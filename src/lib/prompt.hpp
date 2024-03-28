#include <iostream>
#include <stdio.h>
#include <string.h>
#include <limits>
#include <format>
#include "mout.hpp"

using namespace std;

double prompt_double(const char *text)
{
reprompt:
    try
    {
        double value;

        mout << text;

        cin >> value;

        mout << endl;

        if (cin.fail())
        {
            cin.clear();

            cin.ignore(numeric_limits<streamsize>::max(), '\n');

            throw "Error while reading input";
        }

        return value;
    }
    catch (const char *error)
    {
        mout << error << endl
             << endl;

        goto reprompt;
    }
}

int prompt_int(const char *text, int min)
{
reprompt:
    try
    {
        int value;

        mout << text;

        cin >> value;

        mout << endl;

        if (cin.fail())
        {
            cin.clear();

            cin.ignore(numeric_limits<streamsize>::max(), '\n');

            throw "Error while reading input";
        }

        if (value < min)
        {
            cin.clear();

            cin.ignore(numeric_limits<streamsize>::max(), '\n');

            char *error_template = new char[256];

            strcpy(error_template, "Введёное значение менее ожидаемого (%d). Повторите ввод");

            sprintf(error_template, error_template, min);

            throw (const char *)error_template;
        }

        return value;
    }
    catch (const char *error)
    {
        mout << error << endl
             << endl;

        goto reprompt;
    }
}

int prompt_int(const char *text, int min, int max)
{
reprompt:
    try
    {
        int value;

        mout << text;

        cin >> value;

        mout << endl;

        if (cin.fail())
        {
            cin.clear();

            cin.ignore(numeric_limits<streamsize>::max(), '\n');

            throw "Error while reading input";
        }

        if (value < min)
        {
            cin.clear();

            cin.ignore(numeric_limits<streamsize>::max(), '\n');

            char *error_template = new char[256];

            strcpy(error_template, "Введёное значение менее ожидаемого (%d). Повторите ввод");

            sprintf(error_template, error_template, min);

            throw (const char *)error_template;
        }

        if (value > max)
        {
            cin.clear();

            cin.ignore(numeric_limits<streamsize>::max(), '\n');

            char *error_template = new char[256];

            strcpy(error_template, "Введёное значение больше ожидаемого (%d). Повторите ввод");

            sprintf(error_template, error_template, min);

            throw (const char *)error_template;
        }

        return value;
    }
    catch (const char *error)
    {
        mout << error << endl
             << endl;

        goto reprompt;
    }
}
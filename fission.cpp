#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <iterator>
#include <utility>
#include <typeinfo>
#include <fstream>
#include <cstring>
#include "evo.h"

using namespace std;

int main()
{
    Seed(mut::s);
    population environment;
    if (mut::invasion)
    {
        environment = invasion();
    }
    else
    {
        cell one;
        cell two;
        if (mut::start_broken == true)
        {

            ///DO NOT PUT VALUES IN HERE!!!! This is to test if equilibrium is stable
            one.set_cell('A', 1-mut::a_pg, mut::a_pg, mut::a_switch, 1-mut::b_pg, mut::b_pg, mut::b_switch, 0.5);
            two.set_cell('B', 1-mut::a_pg, mut::a_pg, mut::a_switch, 1-mut::b_pg, mut::b_pg, mut::b_switch, 0.5);
        }
        else
        {
            one.set_cell('A', 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5);
            two.set_cell('B', 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5);
        }


        organism new_organism;
        vector<cell> new_vector;
        new_vector.resize(mut::start_size);
        if (mut::start_broken == true)
        {
            for (unsigned int i = 0; i < new_vector.size(); ++i)
            {
                ///ensure unequal frequency for equilibrium check
                {
                    if (i % 5 == 0)
                        new_vector.at(i) = one;
                    else
                        new_vector.at(i) = two;
                }
            }
        }
        else 
        {
            for (unsigned int i = 0; i < new_vector.size(); ++i)
            {
                ///ensure equal frequency of types for initial conditions
                {
                    if (i % 2 == 0)
                        new_vector.at(i) = one;
                    else
                        new_vector.at(i) = two;
                }
            }
        }
        new_organism.set_organism(new_vector);

        vector<organism> total;
        total.resize(mut::start_pop);
        for (unsigned int i = 0; i < total.size(); ++i)
        {
            total.at(i) = new_organism;
        }
        environment.set_population(total);
    }


    environment.process(mut::time_steps);


    return 0;
}

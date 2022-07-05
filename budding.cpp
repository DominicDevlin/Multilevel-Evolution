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
    cell one;
    cell two;
    one.set_code_vector();
    two.set_code_vector();


    if (mut::start_broken == true)
    {
        ///DO NOT PUT VALUES IN HERE!!!! This is to test if equilibrium is stable
        one.set_cell('A', 0.99, 0.01, 0.4, 0.4, 0.6, 0.0005, 0.5);
        two.set_cell('B', 0.99, 0.01, 0.4, 0.4, 0.6, 0.0005, 0.5);
    }
    else
    {
        one.set_cell('A', 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5);
        two.set_cell('B', 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5);
    }

    ///set cells to initial values so that states 'A' and 'B' are equal

    organism new_organism;
    vector<cell> new_vector;
    new_vector.resize(mut::start_size);
    for (unsigned int i = 0; i < new_vector.size(); ++i)
    ///ensure equal frequency of cell states at initiation
    {
        if (i % 2 == 0)
            new_vector.at(i) = one;
        else
            new_vector.at(i) = two;
    }
    new_organism.set_organism(new_vector);

    vector<organism> total;
    total.resize(mut::start_pop);
    for (unsigned int i = 0; i < total.size(); ++i)
    {
        total.at(i) = new_organism;
    }

    population environment;

    environment.set_population(total);
    environment.budding_process(mut::time_steps);

    return 0;
}

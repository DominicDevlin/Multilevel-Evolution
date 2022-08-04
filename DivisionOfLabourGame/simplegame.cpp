#include <iostream>
#include <cmath>
#include <string>
#include <array>
#include <vector>
#include <algorithm>
#include <random>
#include <ctime>
#include <iterator>
#include <utility>
#include <typeinfo>
#include <chrono>
#include <fstream>
#include <cstring>
#include <math.h>
#include <map>


namespace params
{

    const int group_size{20};
    double init_freq{0.16};

    double mut_d{0.05};
    double mut_k_1{1};
    double mut_k_2{0};

    double wt_d{0.047619};
    double wt_k_1{1};
    double wt_k_2{0};


    double alpha{10};



}


auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
std::mt19937 mersenne( static_cast<std::mt19937::result_type>(seed) );

std::uniform_int_distribution<> distro(0, 10);



int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}



class game
{

    std::vector<double> group_frequencies{};

    std::vector<double> mutant_fitness;
    std::vector<double> wt_fitness;


public:



    double generation()
    {
        double mut_freq = params::init_freq;
        int counter{};
        double mk1 = params::mut_k_1;
        double mk2  = params::mut_k_2;
        double md  = params::mut_d; 
        double wk1  = params::wt_k_1;
        double wk2  = params::wt_k_2;
        double wd  = params::wt_d;
        double group = params::group_size;


        for (int i = 0; i < params::group_size + 1; ++i)
        {
            double it = i;

            if (static_cast<int>(group_frequencies.size()) > params::group_size + 1)
                std::cout << "ERROR   " << group_frequencies.size() << std::endl;

            //mutant fitness within group i   - note that trait sociality = 1, does this small difference from my model have an effect?
            double wt_public_good = (wd * pow(wk1, params::alpha) + (1 - wd) * pow(wk2, params::alpha)) * (it / group);

            double mut_public_good = (md * pow(mk1, params::alpha) + (1 - md) * pow(mk2, params::alpha)) * ((group - it) / group);

            double m_fitness = md * (1-mk1) * ((md * pow(mk1, params::alpha) + (1-md) * pow(mk2, params::alpha) - 
            (pow(mk1, params::alpha) * md) / group) * ((group - it) / group)  + wt_public_good) 
            + (1-md) * (1-mk2) * ((md * pow(mk1, params::alpha) + (1-md) * pow(mk2, params::alpha) - 
            (pow(mk2, params::alpha) * md) / group) * ((group - it) / group)  + wt_public_good);
            
            mutant_fitness.push_back(m_fitness);

            //wild type fitness within group i

            double w_fitness = wd * (1-wk1) * ((wd * pow(wk1, params::alpha) + (1-wd) * pow(wk2, params::alpha) - 
            (pow(wk1, params::alpha) * wd) / group) * (it / group)  + mut_public_good) 
            + (1-wd) * (1-wk2) * ((wd * pow(wk1, params::alpha) + (1-wd) * pow(wk2, params::alpha) - 
            (pow(wk2, params::alpha) * wd) / group) * (it / group)  + mut_public_good);

            wt_fitness.push_back(w_fitness);

            std::cout << "group number:  " << i << "   mt fitness:  " << m_fitness << "  wt fitness:  " << w_fitness << std::endl;
            

        }


        while (mut_freq > 0.001 && mut_freq < 0.999)
        {
            double total_mutant_f{};
            double total_wt_f{};

            double sum_m{};
            double sum_wt{};

            //Get total population fitness for mutant and wild type
            for (int i = 0; i < params::group_size + 1; ++i)
            {
                double coefficient = factorial(group) / (static_cast<double>((factorial(i) * factorial(params::group_size - i))));
                double freq = coefficient * pow(mut_freq, group - i) * pow(1 - mut_freq, i);
                
                group_frequencies.push_back(freq);

                if (static_cast<int>(group_frequencies.size()) > params::group_size + 1)
                    std::cout << "ERROR   " << group_frequencies.size() << std::endl;

                total_mutant_f = total_mutant_f + static_cast<double>(params::group_size - (i)) * group_frequencies.at(i) * mutant_fitness.at(i);

                sum_m = sum_m + static_cast<double>(params::group_size - i) * group_frequencies.at(i);

                total_wt_f = total_wt_f + static_cast<double>(i) * group_frequencies.at(i) * wt_fitness.at(i);

                sum_wt = sum_wt + static_cast<double>(i) * group_frequencies.at(i);
            }

            total_mutant_f = total_mutant_f / sum_m;
            total_wt_f = total_wt_f / sum_wt;



            //calculate next generation fitness

            mut_freq = (total_mutant_f * mut_freq) / (total_mutant_f * mut_freq + total_wt_f * (1-mut_freq));

            group_frequencies.clear();
            
            ++ counter;

            if (counter % 1 == 0)
            {
                if (mut_freq > 1)
                    mut_freq = 1;
                if (mut_freq < 0)
                    mut_freq = 0;         
                std::cout << "mutant frequency is  " << mut_freq << std::endl;
                std::string file_name = "data.dat";
                std::ofstream outfile;
                outfile.open(file_name, std::ios::app);
                outfile << counter << "\t" << mut_freq << std::endl;
                outfile.close();
            }

            if (counter > 5000)
                break;


        }
        if (mut_freq > 1)
            mut_freq = 1;
        if (mut_freq < 0)
            mut_freq = 0;


        if (mut_freq > 0.5)
        {
            std::cout << "mutant won   " << mut_freq << std::endl;
            return mut_freq;
        }
        else
        {
            std::cout << "wild type won     " << mut_freq << std::endl;
            return mut_freq;
        }    

    }



};


int main()
{
    // double i = 0.01;
    // while (i < 1)
    // {
    //     params::init_freq = i;
    //     game game1;
    //     double end_freq = game1.generation();

    //     std::string file_name = "init_end.dat";
    //     std::ofstream outfile;
    //     outfile.open(file_name, std::ios::app);
    //     outfile << i << "\t" << end_freq << std::endl;
    //     outfile.close();
    //     i = i + 0.01;

    // }

    game game1;
    game1.generation();



    return 0;
}
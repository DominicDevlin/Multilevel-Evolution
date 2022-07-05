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
#include <sys/timeb.h>
#include "evo.h"

using namespace std;

///RNG
auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
mt19937 mersenne( static_cast<mt19937::result_type>(seed) );
std::normal_distribution<double> curve(0, 0.02);

namespace mut
{
    /// prob of mutation
    const double mut_rate{0.01};

    ///the number of cells in an organism required for fission
    const int threshold{562};

    //starting parameters
    const int start_size{mut::threshold / 2};
    const int start_pop{150};

    ///total cells in the system
    long total_cells{threshold * 80};

    int out_freq{1000};

    ///run steps
    const long time_steps{100000000};

    ///returns on cooperation
    const float alpha{1.0};

    ///returns on personal reproduction
    const float beta{1.0};

        
    /// trait essentiality
    const double trait_e{1.0};

    ///the time step divisor to approximate a continuous model
    const double delta_factor{0.4};

    ///maximum organisms in population - moran process, turn on or off, used for case where there is no selection so there is collective dynamics.
    const bool moran_on{0};

    /// 0 = all traits mutate, 1 = d mutates only, 2 = all traits mutate
    const int selection_type{0};

    //initial conditions if needed. 
    const bool start_broken{false};
    const double a_pg{0};
    const double a_switch{0.5};
    const double b_pg{0.8};
    const double b_switch{0.0005};
    
    ///what data files to send out
    const bool output_rel{true};
    const bool output_var{false};
    const bool output_values{true};
    const bool output_type_diff{false};
    
    // change file names to alpha and beta (returns)
    const bool output_returns{false};
    string al =  to_string(alpha);
    string be = to_string(beta);
    string a = 'a' + al.substr(0,4);
    string b = 'b' + be.substr(0,4);
    string returns =  a + b; 


    ///death rate is calculated based on current cells from the previous time step and an average of the "r" values
    double death_rate{};
    double new_death_rate{};
    long current_cells{};
    long new_current_cells{};

    const int s{-1};

    const bool invasion{false};
    /// 0 = 50/50, 1 = suppressor invades, 2 = fittest invades, 3 = choose frequency
    const int inv_type{3};
    const int inv_freq{5};

    // for budding simulation
    long propagules{};    
}



static int idum = -1;

int Seed(int z)
{
  if (z < 0) 
  {
    int rseed=Randomize();
    return rseed;
  } 
  else 
  {
    idum = static_cast<int>(z);
    for (int i=0; i <100; i++)
        RANDOM();
    return z;
  }
}

int Randomize(void) {
  
  // Set the seed according to the local time

  int val=abs((int)((seed)%655337));
  Seed(val);
  fprintf(stderr,"Random seed is %d\n",val);
  return val;
}


double RANDOM(void)
/* Knuth's substrative method, see Numerical Recipes */
{
  static int inext,inextp;
  static long ma[56];
  static int iff=0;
  long mj,mk;
  int i,ii,k;

  if (idum < 0 || iff == 0) {
    iff=1;
    mj=MSEED-(idum < 0 ? -idum : idum);
    mj %= MBIG;
    ma[55]=mj;
    mk=1;
    i=1;
    do {
      ii=(21*i) % 55;
      ma[ii]=mk;
      mk=mj-mk;
      if (mk < MZ) mk += MBIG;
      mj=ma[ii];
    } while ( ++i <= 54 );
    k=1;
    do {
      i=1;
      do {
        ma[i] -= ma[1+(i+30) % 55];
        if (ma[i] < MZ) ma[i] += MBIG;
      } while ( ++i <= 55 );
    } while ( ++k <= 4 );
    inext=0;
    inextp=31;
    idum=1;
  }
  if (++inext == 56) inext=1;
  if (++inextp == 56) inextp=1;
  mj=ma[inext]-ma[inextp];
  if (mj < MZ) mj += MBIG;
  ma[inext]=mj;
  return mj*FAC;
}

/*! Returns a random integer value between 1 and 'max'
  \param The maximum value (long)
  \return A random integer (long)
**/
long RandomNumber(long max)
{
   return((long)(RANDOM()*max)+1);
}

///mutate the chosen genotype, accounting for boundaries
double normal_mutation(double gene)
{
    gene = gene + curve(mersenne);
    if (gene > 1000)
    {
        gene = 1000;
    }
    if (gene < 0)
    {
        gene = 0;
    }
    return gene;
}

///mutate the chosen genotype log normal, accounting for boundaries
double mutation(double gene)
{

    gene = gene * exp(-curve(mersenne));
    if (gene > 1.0)
    {
        gene = 1.0;
    }
    if (gene < 0.0)
    {
        gene = 0.0;
    }
    return gene;
}

///calculate poisson distributed death rate based on the previous generation cell average, adapted to the size of the organism
int death_calculator(unsigned int organism_size)
{
    double death_calc = mut::death_rate * (static_cast<double>(organism_size) / static_cast<double>(mut::current_cells));  ///calculate death rate based on the previous generation cels, adapted to the size of the organism
    poisson_distribution<int> p_distribution(death_calc); ///poisson distribution for organism deaths
    int deaths = p_distribution(mersenne);
    return deaths;
}

///Deprecated 
int death_calculator_deterministic(unsigned int organism_size)
{
    int deaths = ceil( mut::death_rate * (static_cast<double>(organism_size) / static_cast<double>(mut::current_cells)));  ///calculate death rate based on the previous generation cels, adapted to the size of the organism
    return deaths;
}



void cell::set_code()  ///set a random <double> code at the front of the cell, and remove the one at the back. Used to determine close relatives.
{
    double code = RANDOM();
    codes.insert(codes.begin(), code);
    codes.pop_back();
}


void cell::mutate()
{
    ///if mutation is hit, mutate whole genotype.
    double prob = RANDOM();
    if (prob < mut::mut_rate)
    {
        switch (mut::selection_type)
        {
            case 0:
            {
                m_dummy = normal_mutation(m_dummy);

                ma_replication = mutation(ma_replication);
                ma_catalysis = 1.0 - ma_replication;

                ma_conversion = mutation(ma_conversion);

                mb_replication = mutation(mb_replication);
                mb_catalysis = 1.0 - mb_replication;

                mb_conversion = mutation(mb_conversion);
                break;
            }
            case 1:
            {
                m_dummy = normal_mutation(m_dummy);
                ma_conversion = mutation(ma_conversion);
                break;                
            }
            case 2:
            {
                m_dummy = normal_mutation(m_dummy);
                break;
            
            }
        }
    }
}


void organism::set_organism(vector<cell>& total)
{
    m_current = total;
}

vector<cell>& organism::get_organism()
{
    return m_current;
}

int organism::get_organism_population()
{
    return static_cast<int>(m_current.size());
}



int organism::growth()
{
    /// add to get N
    mut::new_current_cells += m_current.size();

    ///this is public good
    double functionality{};
    ///personal reproduction
    double replication{};

    ///iterate through all cells to determine cellular replication 'rate' = (r/<r> * (1-<r>) * N)
    ///create vector to be used for determining which cells will replicate based on personal fitness
    vector<double> fitness;
    fitness.resize(m_current.size());

    for (unsigned int i = 0; i < m_current.size(); ++i)
    {
        switch (m_current.at(i).get_type())
        {
            case 'A':
                fitness.at(i) = pow(m_current.at(i).get_areplication(), mut::beta);
                functionality += pow(m_current.at(i).get_acatalysis(), mut::alpha);
                replication += pow(m_current.at(i).get_areplication(), mut::beta);
                break;
            case 'B':
                fitness.at(i) = pow(m_current.at(i).get_breplication(), mut::beta);
                functionality += pow(m_current.at(i).get_bcatalysis(), mut::alpha);
                replication += pow(m_current.at(i).get_breplication(), mut::beta);
                break;
            default:
                cout << "Switch Error" << endl;
                break;
        }
    }

    functionality = functionality / m_current.size();
    replication = replication / m_current.size();

    ///discrete distribution based on vector
    double divisions{};
    discrete_distribution<int> total_fitness(fitness.begin(), fitness.end());

    divisions = m_current.size() * (1 - mut::trait_e + mut::trait_e * functionality) * replication * mut::delta_factor;
    mut::new_death_rate += divisions;
    poisson_distribution<int> p_distribution(divisions);
    ///TRYING DIVISIONS HERE
    int replicants = p_distribution(mersenne);
    for (int i = 0; i < replicants; ++i)
    {
        int number = total_fitness(mersenne);
        m_current.push_back(m_current.at(number));
        ///convert the cell state of the push back cell at its conversion probability.
        double conv = RANDOM();
        switch (m_current.back().get_type())
        {
            case 'A':
                if (conv < m_current.back().get_aconversion())
                {
                    m_current.back().set_type('B');
                }
                break;
            case 'B':
                if (conv < m_current.back().get_bconversion())
                {
                    m_current.back().set_type('A');
                }
                break;
            default:
                cout << "Switch Error" << endl;
                break;
        }
        m_current.back().mutate();
    }
    ///proportionate killing of cells in organism
    int deaths = death_calculator(m_current.size());
    for (int i = 0; i < deaths; ++i)
    {
        ///if the organism size is 1, kill the organism.
        if (m_current.size() < 3)
        {
            return 2;
        }
        int number = RandomNumber(m_current.size() - 1)-1;
        m_current.erase(m_current.begin() + number);
    }

    if (static_cast<int>(m_current.size()) >= mut::threshold)
    {
        return 1;
    }
    return 0;

}

vector<cell> organism::split()
{
    vector<cell> new_organism{};
    vector<cell> new_copy{};
    ///randomly allocate each cell to old organism or new organism
    for (cell &i : m_current)
    {
        double number = RandomNumber(2)-1;
        if (number == 1)
        {
            new_copy.push_back(i);
        }
        else
        {
            new_organism.push_back(i);
        }
    }
    m_current = new_copy;
    return new_organism;
}


int organism::budding_growth()
{
    mut::new_current_cells += m_current.size();

    ///guess how many cells will be over the threshold to calculate number of propagules
    int guess = death_calculator(m_current.size());
    int prop_calc = guess + mut::threshold;

    double functionality{};
    double replication{};
    ///iterate through all cells and add their (1-r), then average

    ///this will become the (1-<r>) used to calculate organism growth rate
    std::vector<double> fitness;
    fitness.resize(m_current.size());

    for (unsigned int i = 0; i < m_current.size(); ++i)
    {
        switch (m_current.at(i).get_type())
        {
            case 'A':
                fitness.at(i) = m_current.at(i).get_areplication();
                functionality += m_current.at(i).get_acatalysis();
                replication += m_current.at(i).get_areplication();
                break;
            case 'B':
                fitness.at(i) = m_current.at(i).get_breplication();
                functionality += m_current.at(i).get_bcatalysis();
                replication += m_current.at(i).get_breplication();
                break;
            default:
                std::cout << "Switch Error" << std::endl;
                break;
        }
    }
    functionality = functionality / m_current.size();
    replication = replication / m_current.size();

    std::discrete_distribution<int> total_fitness(fitness.begin(), fitness.end());
    double divisions{};
    divisions = m_current.size() * functionality * replication * mut::delta_factor;
    std::poisson_distribution<int> p_distribution(divisions);
    int replicants = p_distribution(mersenne);
    mut::new_death_rate += replicants;

    for (int i = 0; i < replicants; ++i)
    {
        int number = total_fitness(mersenne);
        switch (static_cast<int>(m_current.size()) > prop_calc)
        {
            ///if the size of the organism is greater than the threshold that is predicted after deaths, cells replicate into propagules instead.
            case true:
            {
                ///create a new organism and push it onto the propagule vector, mutate and convert back cell
                std::vector<cell> new_organism;
                new_organism.push_back(m_current.at(number));
                new_organism.push_back(m_current.at(number));
                double conv = RANDOM();
                switch (new_organism.back().get_type())
                {
                    case 'A':
                        if (conv < new_organism.back().get_aconversion())
                        {
                            new_organism.back().set_type('B');
                        }
                        break;
                    case 'B':
                        if (conv < new_organism.back().get_bconversion())
                        {
                            new_organism.back().set_type('A');
                        }
                        break;
                    default:
                        std::cout << "Switch Error" << std::endl;
                        break;
                }
                new_organism.back().mutate();
                ///erase the cell from the host organism as it moves into a propagule
                m_current.erase(m_current.begin() + number);
                fitness.erase(fitness.begin() + number);
                std::vector<double> &checks = new_organism.at(0).get_codes();
                int sizer = m_current.size();   ///iterate through organism to find related individuals within 60-80 generations. THis number may be too large.
                for (int x = (sizer - 1); x > -1; --x)
                {

                    std::vector<double> &codes = m_current.at(x).get_codes();
                    for (int z = 60; z < 80; ++z)
                    {
                        for (int a = 60; a < 80; ++a)
                        {
                            if (codes.at(a) == checks.at(z))   ///look for identical codes.
                            {
                                new_organism.push_back(m_current.at(x));  
                                m_current.erase(m_current.begin() + x);
                                a = 81;
                                z = 81;
                            }
                            else
                                continue;
                        }
                    }
                    if (new_organism.size() > 50)
                    {
                        x = -1;
                    }

                }
                new_organism.at(0).set_code();
                new_organism.at(1).set_code();
                fitness.resize(m_current.size());
                for (unsigned int q = 0; q < m_current.size(); ++q)   ///redo fitness vector
                {
                    switch (m_current.at(q).get_type())
                    {
                        case 'A':
                            fitness.at(q) = m_current.at(q).get_areplication();
                            break;
                        case 'B':
                            fitness.at(q) = m_current.at(q).get_breplication();
                            break;
                    }
                }
                /// redo fitness distribution
                std::discrete_distribution<int> redo(fitness.begin(), fitness.end());
                total_fitness = redo;
                organism prop;
                ///set the new organism with x related cells
                prop.set_organism(new_organism);
                propagules.push_back(prop);
                break;
            }
            case false:    ///if the organism is not dividing into propagules
            {
                m_current.at(number).set_code();
                ///copy the cell and push it to the back of the vector, then mutate both.
                m_current.push_back(m_current.at(number));
                m_current.back().set_code();
                //m_current.at(number).mutate();

                ///convert the cell state of the push back cell at its conversion probability.
                double conv = RANDOM();
                switch (m_current.back().get_type())
                {
                    case 'A':
                        if (conv < m_current.back().get_aconversion())
                        {
                            m_current.back().set_type('B');
                        }
                        break;
                    case 'B':
                        if (conv < m_current.back().get_bconversion())
                        {
                            m_current.back().set_type('A');
                        }
                        break;
                    default:
                        std::cout << "Switch Error" << std::endl;
                        break;
                }
                m_current.back().mutate();
                break;
            }
        }
    }

    ///proportionate killing of cells in organism
    int deaths = death_calculator(m_current.size());
    for (int i = 0; i < deaths; ++i)
    {
        ///if the organism size is 5(?) or less, kill the organism.
        if (m_current.size() < 3)
        {
            return -1;
        }
        int number = RandomNumber(m_current.size() - 1)-1;
        m_current.erase(m_current.begin() + number);
    }

    if (propagules.size() > 0)
    {
        return propagules.size();
    }
    else
        return 0;

}

///spit the propagules from the host organism.
std::vector<cell> organism::budding_split()
{
    int x = propagules.size();
    for (int i = 0; i < x; ++ i)
    {
        std::vector<cell> org = propagules.at(i).get_organism();
        propagules.erase(propagules.begin() + i);
        return org;
    }
    std::cout << "SPITTING ERROR" << std::endl;
    std::vector <cell> error;
    return error;
}



double get_within_variance(vector<organism>& m_population, double average, double total_cell_count, int pick)
{
    double variance{};

    for (unsigned int i = 0; i < m_population.size(); ++i)
    {
        vector<cell> test = m_population.at(i).get_organism();
        double mean{};
        double xi_xi{};
        for (cell &c : test)
        {
            switch (pick)
            {
            case 1:
                mean += c.get_acatalysis();
                break;
            case 2:
                mean += c.get_aconversion();
                break;
            case 3:
                mean += c.get_bcatalysis();
                break;
            case 4:
                mean += c.get_bconversion();
                break;
            case 5:
                mean += c.get_dummy();
                break;
            }
        }
        mean = mean / test.size();
        for (cell &c : test)
        {
            switch (pick)
            {
            case 1:
                xi_xi += pow(c.get_acatalysis() - mean, 2);
                break;
            case 2:
                xi_xi += pow(c.get_aconversion() - mean, 2);
                break;
            case 3:
                xi_xi += pow(c.get_bcatalysis() - mean, 2);
                break;
            case 4:
                xi_xi += pow(c.get_bconversion() - mean, 2);
                break;
            case 5:
                xi_xi += pow(c.get_dummy() - mean, 2);
                break;
            }
        }
        variance += xi_xi;

    }
    variance = variance / total_cell_count;
    return variance;
}

double get_between_variance(vector<organism>& m_population, double mean, double total_cell_count, int pick)
{
    double variance{};

    for (unsigned int i = 0; i < m_population.size(); ++i)
    {
        double x_x{};
        double between_mean{};
        vector<cell> test = m_population.at(i).get_organism();
        for (cell &c : test)
        {
            switch (pick)
            {
            case 1:
                between_mean += c.get_acatalysis();
                break;
            case 2:
                between_mean += c.get_aconversion();
                break;
            case 3:
                between_mean += c.get_bcatalysis();
                break;
            case 4:
                between_mean += c.get_bconversion();
                break;
            case 5:
                between_mean += c.get_dummy();
                break;
            }
        }
        between_mean = between_mean / test.size();
        x_x += pow(between_mean - mean, 2);

        x_x = x_x * test.size();
        variance += x_x;
    }
    variance = variance / total_cell_count;
    return variance;
}


double get_total_variance(vector<organism>& m_population, double mean, double total_cell_count, int pick)
{
    double variance{};
    double xi_xi{};

    for (unsigned int i = 0; i < m_population.size(); ++i)
    {
        vector<cell> test = m_population.at(i).get_organism();
        for (cell &c : test)
        {
            switch (pick)
            {
            case 1:
                xi_xi += pow(c.get_acatalysis() - mean, 2);
                break;
            case 2:
                xi_xi += pow(c.get_aconversion() - mean, 2);
                break;
            case 3:
                xi_xi += pow(c.get_bcatalysis() - mean, 2);
                break;
            case 4:
                xi_xi += pow(c.get_bconversion() - mean, 2);
                break;
            case 5:
                xi_xi += pow(c.get_dummy() - mean, 2);
            }
        }
    }
    variance = xi_xi / total_cell_count;
    return variance;
}


///calculate averages to send to output. Averages are of all cell values & cell numbers
void get_averages(vector<organism>& m_population, long generations)
{
    double a_catal_average{};
    double b_catal_average{};

    double a_total_average{};
    double b_total_average{};
    double a_conv_average{};
    double b_conv_average{};

    double dummy_average{};

    double c_total_average{};
    long total_cell_count{};

    int total_organisms{};
    long largest_org = {};

    double transmitter_k1{};
    double utiliser_k1{};
    double transmitter_k2{};
    double utiliser_k2{};

    // for invasion analysis
    double fit_freq{};
    double suppress_freq{};

    for (unsigned int i = 0; i < m_population.size(); ++i)
    {
        vector<cell> test = m_population.at(i).get_organism();
        ++ total_organisms;

        if (static_cast<int>(test.size()) > largest_org)
            largest_org = test.size();

        for (cell &c : test)
        {
            switch (c.get_type())
            {
                case 'A':
                {
                    ++ a_total_average;
                    transmitter_k1 += c.get_acatalysis();
                    transmitter_k2 += c.get_bcatalysis();
                    break;
                }
                case 'B':
                {
                    ++ b_total_average;
                    utiliser_k1 += c.get_acatalysis();
                    utiliser_k2 += c.get_bcatalysis();
                    break;
                }
                default:
                {
                    cout << "type error" << endl;
                }

            }
            ++ total_cell_count;
            a_catal_average += c.get_acatalysis();
            a_conv_average += c.get_aconversion();
            b_catal_average += c.get_bcatalysis();
            b_conv_average += c.get_bconversion();
            dummy_average += c.get_dummy();

            fit_freq += 1 - c.get_tag();
            suppress_freq += c.get_tag();
        }
    }
    a_catal_average = a_catal_average / total_cell_count;
    b_catal_average = b_catal_average / total_cell_count;
    a_conv_average = a_conv_average / total_cell_count;
    b_conv_average = b_conv_average / total_cell_count;
    dummy_average = dummy_average / total_cell_count;

    transmitter_k1 = transmitter_k1 / a_total_average;
    transmitter_k2 = transmitter_k2 / a_total_average;
    utiliser_k1 = utiliser_k1 / b_total_average;
    utiliser_k2 = utiliser_k2 / b_total_average;

    a_total_average = a_total_average / total_organisms;
    b_total_average = b_total_average / total_organisms;

    double catal_average = ( b_catal_average * b_total_average + a_catal_average * a_total_average ) / ( a_total_average + b_total_average);


    double k1_within = get_within_variance(m_population, a_catal_average, total_cell_count, 1);
    double k1_between = get_between_variance(m_population, a_catal_average, total_cell_count, 1);
    double k1_total = get_total_variance(m_population, a_catal_average, total_cell_count, 1);
    double k1_relatedness = k1_between / (k1_within + k1_between);

    double d1_within = get_within_variance(m_population, a_conv_average, total_cell_count, 2);
    double d1_between = get_between_variance(m_population, a_conv_average, total_cell_count, 2);
    double d1_relatedness = d1_between / (d1_within + d1_between);

    double k2_within = get_within_variance(m_population, b_catal_average, total_cell_count, 3);
    double k2_between = get_between_variance(m_population, b_catal_average, total_cell_count, 3);
    double k2_relatedness = k2_between / (k2_within + k2_between);

    double d2_within = get_within_variance(m_population, b_conv_average, total_cell_count, 4);
    double d2_between = get_between_variance(m_population, b_conv_average, total_cell_count, 4);
    double d2_relatedness = d2_between / (d2_within + d2_between);


    double within_dummy = get_within_variance(m_population, dummy_average, total_cell_count, 5);
    double between_dummy = get_between_variance(m_population, dummy_average, total_cell_count, 5);
    double full_dummy = get_total_variance(m_population, dummy_average, total_cell_count, 5);
    double relatedness = between_dummy / (between_dummy + within_dummy);

    double global_average_r = catal_average * (1 - catal_average);



    if (mut::output_values == true)
    {
        string file_name{};
        if (mut::output_returns == true)
        {
            file_name = to_string(mut::threshold) + mut::returns + ".dat";
        }
        else
            file_name = to_string(mut::threshold) + ".dat";
        ofstream outfile;
        outfile.open(file_name, ios::app);
        outfile << generations << '\t' << a_total_average / (a_total_average + b_total_average) << '\t' << a_total_average << '\t' << a_catal_average << '\t' << a_conv_average << '\t'
        << b_total_average << '\t' << b_catal_average << '\t' << b_conv_average << endl;
        outfile.close();
    }


    if (mut::output_rel == true)
    {
        string var_name{};
        if (mut::output_returns == true)
            var_name = to_string(mut::threshold) + mut::returns + "rel.dat";
        else
            var_name = to_string(mut::threshold) + "rel.dat";
        ofstream outfile;
        outfile.open(var_name, ios::app);
        outfile << global_average_r << '\t' << relatedness << '\t' << k1_relatedness << '\t' << d1_relatedness 
        << '\t' << k2_relatedness << '\t' << d2_relatedness << endl;
        outfile.close();
    }
    if (mut::output_var == true)
    {
        string var_name = to_string(mut::threshold) + "var.dat";
        ofstream outfile;
        outfile.open(var_name, ios::app);
        outfile << full_dummy << '\t' << within_dummy << '\t' << between_dummy << '\t' << d1_within << '\t' << d1_between << endl;
        outfile.close();
    }

    if (mut::output_type_diff == true)
    {
        string var_name = to_string(mut::threshold) + "type_diff.dat";
        ofstream outfile;
        outfile.open(var_name, ios::app);
        outfile << a_total_average / (a_total_average + b_total_average) << '\t' << transmitter_k1 << '\t' 
        << transmitter_k2 << '\t' << utiliser_k1 << '\t' << utiliser_k2 << endl;
        outfile.close();
    }

    if (mut::invasion)
    {
        fit_freq = fit_freq / total_cell_count;
        suppress_freq = suppress_freq / total_cell_count;

        string f_name = to_string(mut::threshold) + "invasion.dat";
        ofstream outfile;
        outfile.open(f_name, ios::app);
        outfile << fit_freq << '\t' << suppress_freq << std::endl;
    }

    cout << largest_org <<  "     population size: " << total_cell_count << "   number of orgs: " << m_population.size() <<
    "   approx number of cells before death: " << mut::current_cells + mut::new_death_rate << endl;


}


void population::set_population(vector<organism>& total)
{
    m_population = total;
    if (mut::invasion)
    {
        mut::out_freq = 10;
    }
}

vector<organism> population::get_population()
{
    return m_population;
}

void population::process(long time_steps)
{
    const long total_steps = time_steps;
    ///loop for total time steps
    while (time_steps > 0)
    {
        ///reset the death calculators which changes each time step
        mut::new_current_cells = 0;
        mut::new_death_rate = 0;
        int x = m_population.size();
        for (int i = 0; i < x; ++i)
        {
            ///if check split returns 1, split the organism and push the new one to the back of m_popuation
            ///if organism size is 1 or less, erase the organism from the vector
            int check_split = m_population.at(i).growth();
            switch (check_split)
            {
                case 2:
                {
                    m_population.erase(m_population.begin() + i);
                    -- x;
                    -- i;
                    break;
                }
                case 1:
                {
                    vector<cell> new_vector = m_population.at(i).split();
                    organism new_organism;
                    new_organism.set_organism(new_vector);
                    ///push new organisms to the back
                    m_population.push_back(new_organism);
                    break;
                }
            }
        }
        /// moran process for case with no selection.
        if (mut::new_current_cells > (mut::total_cells * 0.8) && mut::moran_on == true)
        {
            int num = RandomNumber(m_population.size()-1)-1;
            m_population.erase(m_population.begin() + num);
        }

        ///calculate next time step death rate ((<r> * N^2) / K)
        mut::current_cells = mut::new_current_cells;
        mut::death_rate =  (mut::new_death_rate * mut::current_cells) / mut::total_cells;

        // if (time_steps % 200 == 0)
        // {
        //     cout << total_steps - time_steps << " time steps have passed" << endl;
        // }
        
        if (time_steps % mut::out_freq == 0)
        {
            get_averages(m_population, total_steps - time_steps);
        }
        --time_steps;
    }
}


void population::budding_process(long time_steps)
{
    const long total_steps = time_steps;
    ///loop for total time steps
    while (time_steps > 0)
    {
        ///reset the death calculators which changes each time step
        mut::new_current_cells = 0;
        mut::new_death_rate = 0;
        mut::propagules = 0;
        int x = m_population.size();
        for (int i = 0; i < x; ++i)
        {
            int check_split = m_population.at(i).budding_growth();
            ///if -1 is returned, erase organism, if more than 1 is returned, spit the propagules which have reached greater than 10 in size
            switch (check_split)
            {
                case -1:
                    m_population.at(i) = m_population.at(x - 1);
                    m_population.erase(m_population.begin() + (x - 1));
                    -- x;
                    -- i;
                    break;
                case 0:
                    break;
                default:
                    ///split propagules, iterate for how many in the organism have reached 10+
                    for (int y = 0; y < check_split; ++y)
                    {
                        std::vector<cell> new_vector = m_population.at(i).budding_split();
                        organism new_organism;
                        new_organism.set_organism(new_vector);
                        ///push new organisms to the back
                        m_population.push_back(new_organism);
                        // std::cout << "Propagule has " << new_vector.size() << " cells" << std::endl;
                    }
                    break;
            }
        }
        ///calculate next time step death rate ((<r> * N^2) / K)
        mut::current_cells = mut::new_current_cells;
        mut::death_rate = (mut::new_death_rate * mut::current_cells) / mut::total_cells;

        if (time_steps % 1000 == 0)
        {
            get_averages(m_population, total_steps - time_steps);
        }
        --time_steps;
    }
}




population invasion()
{
    cell fit_one;
    cell fit_two;
    cell suppress_one;
    cell suppress_two;

    fit_one.set_cell('A', 0.99, 0.01, 0.5, 0.01, 0.99, 0.0005, 0.5, false);
    fit_two.set_cell('B', 0.99, 0.01, 0.5, 0.01, 0.99, 0.0005, 0.5, false);

    suppress_one.set_cell('A', 0.99, 0.01, 0.5, 0.47, 0.53, 0.0005, 0.5, true);
    suppress_two.set_cell('B', 0.99, 0.01, 0.5, 0.47, 0.53, 0.0005, 0.5, true);

    organism fit;
    organism suppress;

    std::vector<cell> fit_vector;
    fit_vector.resize(mut::start_size);
    for (unsigned int i = 0; i < mut::start_size; ++i)
    {
        ///ensure equal frequencies for max fitness
        {
            if (i % 2 == 0)
                fit_vector.at(i) = fit_one;
            else
                fit_vector.at(i) = fit_two;
        }
    }    
    std::vector<cell> suppress_vector;
    suppress_vector.resize(mut::start_size);
    for (unsigned int i = 0; i < mut::start_size; ++i)
    {
        ///ensure equal frequencies for max fitness
        {
            if (i % 10 == 0)
                suppress_vector.at(i) = suppress_one;
            else
                suppress_vector.at(i) = suppress_two;
        }
    }
    fit.set_organism(fit_vector);
    suppress.set_organism(suppress_vector);

    /// want total number of cells to be just below carrying capacity.
    std::vector<organism> org_list;


    int n_orgs = static_cast<int>(ceil((mut::total_cells * 0.85) / mut::start_size));

    std::vector<organism> total;
    total.resize(n_orgs);

    for (int i=0; i<n_orgs; ++i)
    {
        switch(mut::inv_type)
        {   
            ///50/50 population
            case 0:
            {
                if (i % 2 == 0)
                    total.at(i) = suppress;
                else
                    total.at(i) = fit;
                break;
            }
            /// suppressor invades
            case 1:
            {
                if (i == 0)
                    total.at(i) = suppress;
                else
                    total.at(i) = fit;
                break;
            }
            /// fit invades
            case 2:
            {
                if (i == 0)
                    total.at(i) = fit;
                else
                    total.at(i) = suppress;
                break;
            }
            /// specific frequency (must be integer)
            case 3:
            {
                if (i % mut::inv_freq == 0)
                    total.at(i) = suppress;
                else
                    total.at(i) = fit;
                break;
            }
        }
    }
    population env;
    env.set_population(total);
    return env;
}
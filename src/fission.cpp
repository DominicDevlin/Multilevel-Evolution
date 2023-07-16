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
#include "omp.h"
#include <sys/timeb.h>
#include <sstream>
#include <sys/stat.h>



using namespace std;


namespace mut
{


    /// mut_rate is the chance of mutation upon reproduction
    double mut_rate{0.00316};
    double m{0};
    double sd{0.02};

    ///the number of cells in an organism required for fission
    int threshold{1778};

    //starting parameters
    int start_size{mut::threshold / 2};
    constexpr int start_pop{150};

    // max cells 
    int capacity{threshold * 300};
    // max cells (carrying capacity if starting parallel)
    int array_tc{300};

    ///run steps
    constexpr long time_steps{500000000};

    //number of threads
    const int n_threads{6};

    // save simulation at set time points.
    const bool savestates{true};
    const long savet{10000};
    //load from a save state
    const bool loadstate{false};
    const string SaveName{"1000ss40000.dat"};

    ///returns on cooperation
    float alpha{0.9};

    ///returns on personal reproduction
    constexpr float beta{1.0};

    //change the file name to include alpha/beta parameters
    constexpr bool output_returns{false};

    // sweep through alpha instead of mutation rate, used for alpha-V phase.
    constexpr bool alpha_sweep{true};

    
    string al =  to_string(alpha);
    string be = to_string(beta);
    string a = 'a' + al.substr(0,4);
    string b = 'b' + be.substr(0,4);
    string returns =  a + b; 

    /// trait essentiality
    constexpr double trait_e{1.0};

    ///the time step divisor to approximate a continuous model
    constexpr double delta_factor{0.4};

    ///maximum organisms in population - moran process, turn on or off, used for case where there is no selection so there is collective dynamics.
    constexpr bool moran_on{0};

    /// 0 = all traits mutate, 1 = d mutates only, 2 = only dummy mutates, 3 = public good mutate
    constexpr int selection_type{0};


    bool choose_start{true};
    double a_pg{0.24};
    double a_switch{0.27};
    double b_pg{0.47};
    double b_switch{0.0005};   
    
    constexpr bool start_broken_frequencies{false};
    
    bool sweep = true;
    string output_name = to_string(threshold);

    ///what data files to send out
    constexpr bool output_rel{true};
    constexpr bool output_var{true};
    constexpr bool output_values{true};
    constexpr bool output_type_diff{false};


    ///death rate is calculated based on current cells from the previous time step and an average of the "r" values
    double death_rate{};
    double new_death_rate{};
    long current_cells{};
    long new_current_cells{};



    const bool invasion{false};
    /// 0 = 50/50, 1 = suppressor invades, 2 = fittest invades, 3 = choose frequency
    const int inv_type{0};
    // ratio of a to b
    const int inv_freq{1};

    int out_freq{1000};


}

///RNG
auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
mt19937 mersenne( static_cast<mt19937::result_type>(seed) );


class xoshiro
{
private:
    uint64_t state[4];
public:

    void seed(uint64_t n)
    {
        // Need to seed correctly. 
        state[0] = splitmix(n);
        state[1] = splitmix(state[0]);
        state[2] = splitmix(state[1]);
        state[3] = splitmix(state[2]);
    }


    uint64_t splitmix(uint64_t x) 
    {
        uint64_t z = (x += 0x9e3779b97f4a7c15);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
        z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
        return z ^ (z >> 31);
    }


    // thread safe xooshiro256**
    uint64_t rotl(uint64_t x, int k)
    {
        return (x << k) | (x >> (64 - k));
    }


    uint64_t operator()() 
    {
        const uint64_t result = rotl(state[1] * 5, 7) * 9;

        const uint64_t t = state[1] << 17;

        state[2] ^= state[0];
        state[3] ^= state[1];
        state[1] ^= state[2];
        state[0] ^= state[3];

        state[2] ^= t;

        state[3] = rotl(state[3], 45);

        return result;
    }

};

    double double_num(xoshiro& rng) 
    {
        uint64_t x = rng();
        const union { uint64_t i; double d; } u = {.i = UINT64_C(0x3FF) << 52 | x >> 12 };
        return u.d - 1.0;
    }


    int RandomInt(xoshiro& rng, int min, int max)
    {
        return (std::floor(double_num(rng)*(max-min+1)+min));
    }

    // int RandomInt(xoshiro& rng, int min, int max)
    // {
    //     uniform_int_distribution<int> distribution(min, max);
    //     return distribution(rng);
    // }

    double curve(xoshiro &rng)
    {
        double U1 = double_num(rng);
        double U2 = double_num(rng);
        double Z1 = sqrt(-2*log(U1)) * cos(2 * M_PI * U2);
        return mut::m + mut::sd * Z1;
    }


    double normal_dist(xoshiro& rng, double mean, double sdev)
    {
        double U1 = double_num(rng);
        double U2 = double_num(rng);
        double Z1 = sqrt(-2*log(U1)) * cos(2 * M_PI * U2);
        return mean + sdev * Z1;
    }


    int poisson(double mean, xoshiro& rng) 
    {
        double L = exp(-mean);
        double p = 1.0;
        int k = 0;
        do
        {
            p*= double_num(rng);
            k++;
        } while (p > L);

        return k - 1;
    }




///mutate the chosen genotype, only used for dummy variable (1000 boundary)
void normal_mutation(double &gene, xoshiro& rng)
{
    gene = gene + curve(rng);
    if (gene > 1000)
    {
        gene = 1000;
    }
    if (gene < 0)
    {
        gene = 0;
    }
}

///mutate the chosen genotype log normal, accounting for boundaries
void mutation(double &gene, xoshiro& rng)
{

    gene = gene * exp(-curve(rng));
    if (gene > 1.0)
    {
        gene = 1.0;
    }
}


///mutate the chosen genotype linearly, accounting for wide boundaries
void lin_mutation(double &gene, xoshiro& rng)
{
    gene = gene + curve(rng);
    if (gene > 1.1)
    {
        gene = 1.1;
    }
    else if (gene < -0.1)
    {
        gene = -0.1;
    }
}


void spread_mutation(double &gene, xoshiro& rng)
{
    gene = gene * exp(normal_dist(rng, 0, 0.05));
    if (gene > 1.0)
    {
        gene = 1.0;
    }    
}


///calculate poisson distributed death rate based on the previous generation cell average, normalised to the size of the organism
int death_calculator(unsigned int organism_size, xoshiro& rng)
{
    double death_calc = mut::death_rate * (static_cast<double>(organism_size) / static_cast<double>(mut::current_cells));  ///calculate death rate based on the previous generation cels, adapted to the size of the organism
    int deaths = poisson(death_calc, rng);
    return deaths;
}

///depracated
int death_calculator_deterministic(unsigned int organism_size)
{
    int deaths = floor( mut::death_rate * (static_cast<double>(organism_size) / static_cast<double>(mut::current_cells)));  ///calculate death rate based on the previous generation cels, adapted to the size of the organism
    return deaths;
}


class cell
{
    /// cell type can either be 'A' or 'B', each has 3 values (replication = '1-k', catalysis = 'k', conversion = 'd')
    char m_type;
    double ma_replication;
    double ma_catalysis;
    double ma_conversion;
    double mb_replication;
    double mb_catalysis;
    double mb_conversion;

    bool m_tag;

    double m_dummy;


    double alpha_ka;
    double alpha_kb;
    double beta_ra;
    double beta_rb;

    double xa_replication;
    double xa_conversion;
    double xb_replication;
    double xb_conversion;


public:
    void set_cell(char type, double replication1, double catalysis1, double conversion1, double replication2, double catalysis2, double conversion2, double dummy, bool tag=false)
    {
        // parameters
        m_type = type;
        ma_replication = replication1;
        ma_conversion = conversion1;
        ma_catalysis = catalysis1;
        mb_replication = replication2;
        mb_conversion = conversion2;
        mb_catalysis = catalysis2;
        m_dummy = dummy;
        m_tag = tag;

        // variables used only for linear mutation to allow a trait range of < 0 or >1
        xa_replication = replication1;
        xa_conversion = conversion1;
        xb_replication = replication2;
        xb_conversion = conversion2;

        // if alpha or beta !=1, the powers are done here.
        alpha_ka = pow(ma_catalysis, mut::alpha);
        beta_ra = pow(ma_replication, mut::beta);
        alpha_kb = pow(mb_catalysis, mut::alpha);
        beta_rb = pow(mb_replication, mut::beta);

    }

    double get_tag()
    {
        return m_tag;
    }

    double get_dummy()
    {
        return m_dummy;
    }

    void set_type(char type)
    {
        m_type = type;
    }

    double get_aconversion()
    {
        return ma_conversion;
    }

    double get_areplication()
    {
        return ma_replication;
    }

    double get_acatalysis()
    {
        return ma_catalysis;
    }


    double get_bconversion()
    {
        return mb_conversion;
    }

    double get_breplication()
    {
        return mb_replication;
    }

    double get_bcatalysis()
    {
        return mb_catalysis;
    }

    char get_type()
    {
        return m_type;
    }

    double get_alpha_ka()
    {
        return alpha_ka;
    }

    double get_alpha_kb()
    {
        return alpha_kb;
    }

    double get_beta_ra()
    {
        return beta_ra;
    }

    double get_beta_rb()
    {
        return beta_rb;
    }

    void mutate(xoshiro& s)
    {
        ///if mutation is hit, mutate whole genotype.
        double prob = double_num(s);
        if (prob < mut::mut_rate)
        {
            switch (mut::selection_type)
            {
                case 0: // regular simulation
                {
                    normal_mutation(m_dummy, s);

                    mutation(ma_replication, s);
                    ma_catalysis = 1.0 - ma_replication;

                    alpha_ka = pow(ma_catalysis, mut::alpha);
                    beta_ra = pow(ma_replication, mut::beta);

                    mutation(ma_conversion, s);

                    mutation(mb_replication, s);
                    mb_catalysis = 1.0 - mb_replication;

                    alpha_kb = pow(mb_catalysis, mut::alpha);
                    beta_rb = pow(mb_replication, mut::beta);

                    mutation(mb_conversion, s);
                    break;
                }
                case 1: // differentiation test
                {
                    normal_mutation(m_dummy, s);
                    mutation(ma_conversion, s);
                    break;                
                }
                case 2: // all traits equal, except dummy to measure relatedness. Need moran process with this for turnover
                {
                    normal_mutation(m_dummy, s);
                    break;
                
                } 
                case 3: // mutating k instead of 1-k
                {
                    
                    normal_mutation(m_dummy, s);

                    mutation(ma_catalysis, s);
                    ma_replication = 1.0 - ma_catalysis;
                    mutation(ma_conversion, s);

                    alpha_ka = pow(ma_catalysis, mut::alpha);
                    beta_ra = pow(ma_replication, mut::beta);

                    mutation(mb_catalysis, s);
                    mb_replication = 1.0 - mb_catalysis;
                    mutation(mb_conversion, s);

                    alpha_kb = pow(mb_catalysis, mut::alpha);
                    beta_rb = pow(mb_replication, mut::beta);

                    break;
                }
                case 4: // mutating 1 - d
                {
                    normal_mutation(m_dummy, s);

                    mutation(ma_replication, s);
                    ma_catalysis = 1.0 - ma_replication;

                    double ma_1d = 1.0 - ma_conversion;
                    mutation(ma_1d, s);
                    ma_conversion = 1.0 - ma_1d;

                    mutation(mb_replication, s);
                    mb_catalysis = 1.0 - mb_replication;

                    double mb_1d = 1.0 - mb_conversion;
                    mutation(mb_1d, s);
                    mb_conversion = 1.0 - mb_1d;
                    
                    alpha_ka = pow(ma_catalysis, mut::alpha);
                    beta_ra = pow(ma_replication, mut::beta);
                    alpha_kb = pow(mb_catalysis, mut::alpha);
                    beta_rb = pow(mb_replication, mut::beta);
                    break;
                }
                case 5: // linear mutation
                {
                    normal_mutation(m_dummy, s);

                    lin_mutation(xa_replication, s);

                    ma_replication = xa_replication;
                    if (ma_replication > 1.)
                        ma_replication = 1.;
                    else if (ma_replication < 0.)
                        ma_replication = 0.;

                    ma_catalysis = 1.0 - ma_replication;

                    lin_mutation(xa_conversion, s);
                    ma_conversion = xa_conversion;
                    if (ma_conversion > 1.)
                        ma_conversion = 1.;
                    else if (ma_conversion < 0.)
                        ma_conversion = 0.;

                    lin_mutation(xb_replication, s);

                    mb_replication = xb_replication;
                    if (mb_replication > 1.)
                        mb_replication = 1.;
                    else if (mb_replication < 0.)
                        mb_replication = 0.;

                    mb_catalysis = 1.0 - mb_replication;

                    lin_mutation(xb_conversion, s);
                    mb_conversion = xb_conversion;
                    if (mb_conversion > 1.)
                        mb_conversion = 1.;
                    else if (mb_conversion < 0.)
                        mb_conversion = 0.;

                    alpha_ka = pow(ma_catalysis, mut::alpha);
                    beta_ra = pow(ma_replication, mut::beta);
                    alpha_kb = pow(mb_catalysis, mut::alpha);
                    beta_rb = pow(mb_replication, mut::beta);
                    break;
                }
            }
        }
    }
    void spread(xoshiro& s) // variation of traits at beginning of simulation.
    {
        normal_mutation(m_dummy, s);
        spread_mutation(ma_replication, s);
        ma_catalysis = 1.0 - ma_replication;

        spread_mutation(ma_conversion, s);

        spread_mutation(mb_replication, s);
        mb_catalysis = 1.0 - mb_replication;

        spread_mutation(mb_conversion, s);
    }
};



///cells in organism are contained within 'm_current'
class organism
{
    vector<cell> m_current;
    /// check how many mutations are happening per life cycle
    long mut_num{};
    xoshiro rng;
    int n_births{};


public:

    organism()
    {
        m_current.reserve(mut::threshold);
    }


    void seed_org(uint64_t s)
    {
        rng.seed(s);
    }


    void set_organism(vector<cell>& total)
    {
        m_current = total;
    }

    vector<cell>& get_organism()
    {
        return m_current;
    }

    int get_organism_population()
    {
        return static_cast<int>(m_current.size());
    }

    int get_rN()
    {
        return n_births;
    }

    int growth()
    {

        /// add to get N
        int org_size = m_current.size();
        ///this is public good
        double functionality{};
        ///personal reproduction
        double replication{};

        ///iterate through all cells to determine cellular replication 'rate' = (r/<r> * (1-<r>) * N)
        ///create vector to be used for determining which cells will replicate based on personal fitness
        vector<double> fitness{};
        fitness.resize(org_size);
        for (unsigned int i = 0; i < org_size; ++i)
        {
            switch (m_current.at(i).get_type())
            {
                case 'A':
                    fitness.at(i) = m_current.at(i).get_beta_ra();
                    functionality += m_current.at(i).get_alpha_ka();
                    replication += m_current.at(i).get_beta_ra();
                    break;
                case 'B':
                    fitness.at(i) = m_current.at(i).get_beta_rb();
                    functionality += m_current.at(i).get_alpha_kb();
                    replication += m_current.at(i).get_beta_rb();
                    break;
                default:
                    cout << "Switch Error" << endl;
                    break;
            }
        }

        functionality = functionality / (double)org_size;
        replication = replication / (double)org_size;

        double sum_numbers = accumulate(fitness.begin(), fitness.end(), 0.);
        vector<double> probabilities(org_size);
        for (int i = 0; i < org_size; ++i)
        {
            probabilities[i] = fitness[i] / sum_numbers;
        }

        vector<double> cdf(org_size);
        partial_sum(probabilities.begin(), probabilities.end(), cdf.begin());
        double divisions = org_size * (1 - mut::trait_e + mut::trait_e * functionality) * replication * mut::delta_factor;

        int replicants = poisson(divisions, rng);
        n_births = replicants;

        for (int i = 0; i < replicants; ++i)
        {
            // int number = total_fitness(rng);
            double dnum = double_num(rng);
            auto it_num = upper_bound(cdf.begin(), cdf.end(), dnum);
            int number = distance(cdf.begin(), it_num);
            if (number >= org_size)
            {
             	cout << "ERROR WITH DISCRETE DISTRIBUTION" << endl;
            }

            m_current.push_back(m_current.at(number));
            ///convert the cell state of the push back cell at its conversion probability.
            double conv = double_num(rng);
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
            m_current.back().mutate(rng);
        }

        ///proportionate killing of cells in organism
        int deaths = death_calculator(m_current.size(), rng); 
        // cout << deaths << endl;

        // int deaths = death_calculator_deterministic(m_current.size());
        for (int i = 0; i < deaths; ++i)
        {
            ///if the organism size is 1, kill the organism.
            if (m_current.size() < 2)
            {
                return 2;
            }
            int number = RandomInt(rng, 0, m_current.size() - 1);
            auto it = m_current.begin() + number;
            *it = move(m_current.back());
            m_current.pop_back();
        }
        if (static_cast<int>(m_current.size()) >= mut::threshold)
        {
            return 1;
        }
        return 0;

    }


    ///split the organism into two. return the new organism
    vector<cell> split()
    {

        vector<cell> new_organism{};
        vector<cell> new_copy{};
        ///randomly allocate each cell to old organism or new organism
        for (cell &i : m_current)
        {
            double number = RandomInt(rng, 0, 1);
            if (number == 1)
            {
                new_copy.push_back(i);
            }
            else
            {
                new_organism.push_back(i);
            }
        }
        if (new_copy.empty())
        {
            new_copy.push_back(m_current[0]);
        }
        else if (new_organism.empty())
        {
            new_organism.push_back(m_current[0]);
        }
        
        m_current = new_copy;
        return new_organism;
    }
};


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

    double pg{};
    double repl{};

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
                    
                    pg += c.get_alpha_ka();
                    repl += c.get_beta_ra();
                    break;
                }
                case 'B':
                {
                    ++ b_total_average;
                    utiliser_k1 += c.get_acatalysis();
                    utiliser_k2 += c.get_bcatalysis();

                    pg += c.get_alpha_kb();
                    repl += c.get_beta_rb();
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

    pg = pg / total_cell_count;
    repl = repl / total_cell_count;


    double global_r = (1 - mut::trait_e + mut::trait_e * pg) * repl;


    double prop1 = a_total_average / (a_total_average + b_total_average);
    double prop2 = 1 - prop1;


    a_total_average = a_total_average / total_organisms;
    b_total_average = b_total_average / total_organisms;



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



    double sup_gr{};
    double fit_gr{};
    double sup_rel{};
    double fit_rel{};

    if (mut::invasion)
    {
        vector<organism> cp_fit;
        vector<organism> cp_sup;

        double sa_catal_average{};
        double sb_catal_average{};
        double sa_total_average{};
        double sb_total_average{};

        double fa_catal_average{};
        double fb_catal_average{};
        double fa_total_average{};
        double fb_total_average{};

        double fa_rep{};
        double fb_rep{};
        double sa_rep{};
        double sb_rep{};


        double fit_avg{};
        long n_fit{};
        double sup_avg{};
        long n_sup{};
        for (unsigned int i=0; i < m_population.size(); ++i)
        {
            bool ctag=false;
            vector<cell> org = m_population.at(i).get_organism();
            organism new_org;
            new_org.set_organism(org);
            if (org.at(0).get_tag())
            {
                cp_sup.push_back(new_org);
                ctag=true;
            }
            else
                cp_fit.push_back(new_org);   

            if (ctag)
                for (cell &c : org)
                {
                    fit_avg += c.get_dummy();
                    ++n_fit;

                    switch (c.get_type())
                    {
                        case 'A':
                        {
                            ++ sa_total_average;
                            sa_catal_average += pow(c.get_acatalysis(), mut::alpha);
                            sa_rep += 1. - c.get_acatalysis();
                            break;
                        }
                        case 'B':
                        {
                            ++ sb_total_average;
                            sb_catal_average += pow(c.get_bcatalysis(), mut::alpha);
                            sb_rep += 1. - c.get_bcatalysis();
                            break;
                        }
                        default:
                        {
                            cout << "type error" << endl;
                        }

                    }
                } 
            else
                for (cell &c : org)
                {
                    sup_avg += c.get_dummy();            
                    ++n_sup;
                    switch (c.get_type())
                    {
                        case 'A':
                        {
                            ++ fa_total_average;
                            fa_catal_average += pow(c.get_acatalysis(), mut::alpha);
                            fa_rep += 1. - c.get_acatalysis();
                            break;
                        }
                        case 'B':
                        {
                            ++ fb_total_average;
                            fb_catal_average += pow(c.get_bcatalysis(), mut::alpha);
                            fb_rep += 1. - c.get_bcatalysis();
                            break;
                        }
                        default:
                        {
                            cout << "type error" << endl;
                        }

                    }
                } 
        }

        sa_rep = sa_rep / sa_total_average;
        sb_rep = sb_rep / sb_total_average;
        fa_rep = fa_rep / fa_total_average;
        fb_rep = fb_rep / fb_total_average;

        // cout << sa_rep << "  " << sb_rep << "  " <<  fa_rep << "  " << fb_rep << endl;
        
        sa_catal_average = sa_catal_average / sa_total_average;
        sb_catal_average = sb_catal_average / sb_total_average;

        fa_catal_average = fa_catal_average / fa_total_average;
        fb_catal_average = fb_catal_average / fb_total_average;

        sa_total_average = sa_total_average / static_cast<double>(cp_sup.size());
        sb_total_average = sb_total_average / static_cast<double>(cp_sup.size());

        fa_total_average = fa_total_average / static_cast<double>(cp_fit.size());
        fb_total_average = fb_total_average / static_cast<double>(cp_fit.size());


        double scatal_average = ( sb_catal_average * sb_total_average + sa_catal_average * sa_total_average ) / ( sa_total_average + sb_total_average);
        double fcatal_average = ( fb_catal_average * fb_total_average + fa_catal_average * fa_total_average ) / ( fa_total_average + fb_total_average);

        double srep = ( sa_rep * sa_total_average + sb_rep * sb_total_average) / (sa_total_average + sb_total_average);
        double frep = ( fa_rep * fa_total_average + fb_rep * fb_total_average) / (fa_total_average + fb_total_average);


        sup_gr = srep * scatal_average;
        fit_gr = frep * fcatal_average;



        fit_avg = fit_avg / n_fit;
        sup_avg = sup_avg / n_sup;

        double w_dfit = get_within_variance(cp_fit, fit_avg, n_fit, 5);
        double b_dfit = get_between_variance(cp_fit, fit_avg, n_fit, 5);
        fit_rel = b_dfit / (w_dfit + b_dfit);

        double w_dsup = get_within_variance(cp_sup, sup_avg, n_sup, 5);
        double b_dsup = get_between_variance(cp_sup, sup_avg, n_sup, 5);
        sup_rel = b_dsup / (w_dsup + b_dsup);   

        // frequencies from first loop
        fit_freq = fit_freq / total_cell_count;
        suppress_freq = suppress_freq / total_cell_count;
        
        double gr_average = (sup_gr * n_sup + fit_gr * n_fit) / static_cast<double>(total_cell_count);

        // double org_gen_time = log(2) / ()


        // double avg_m = (fit_freq * fit_gr + suppress_freq * sup_gr);
        
        double sup_fitness = exp((sup_gr-gr_average));// / suppress_freq;
        double fit_fitness = exp((fit_gr-gr_average));// / fit_freq;

        string f_name = to_string(mut::threshold) + to_string(mut::alpha) + "invasion.dat";
        ofstream outfile;
        outfile.open(f_name, ios::app);
        outfile << fit_freq << '\t' << suppress_freq << '\t' << fit_rel << '\t' << sup_rel << '\t' << fit_gr << '\t' << sup_gr << endl;
        
    } 






    // total population calculations. 
    double within_dummy = get_within_variance(m_population, dummy_average, total_cell_count, 5);
    double between_dummy = get_between_variance(m_population, dummy_average, total_cell_count, 5);
    double full_dummy = get_total_variance(m_population, dummy_average, total_cell_count, 5);
    double relatedness = between_dummy / (between_dummy + within_dummy);


    if (mut::output_values == true)
    {
        string file_name{};
        if (mut::output_returns == true)
        {
            file_name = to_string(mut::threshold) + mut::returns + ".dat";
        }
        else
            file_name = mut::output_name + ".dat";
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
            var_name = mut::output_name + "rel.dat";
        ofstream outfile;
        outfile.open(var_name, ios::app);
        outfile << global_r << '\t' << relatedness << '\t' << k1_relatedness << '\t' << d1_relatedness 
        << '\t' << k2_relatedness << '\t' << d2_relatedness << endl;
        outfile.close();
    }
    if (mut::output_var == true)
    {
        string var_name{};
        if (mut::output_returns == true)
            var_name = to_string(mut::threshold) + mut::returns + "var.dat";
        else
            var_name = mut::output_name + "var.dat";
        ofstream outfile;
        outfile.open(var_name, ios::app);
        outfile << full_dummy << '\t' << within_dummy << '\t' << between_dummy << '\t' << d1_within << '\t' << d1_between << endl;
        outfile.close();
    }

    if (mut::output_type_diff == true)
    {
        string var_name{};
        if (mut::output_returns == true)
            var_name = to_string(mut::threshold) + mut::returns + "type_diff.dat";
        else
            var_name = mut::output_name + "type_diff.dat";
        
        ofstream outfile;
        outfile.open(var_name, ios::app);
        outfile << a_total_average / (a_total_average + b_total_average) << '\t' << transmitter_k1 << '\t' 
        << transmitter_k2 << '\t' << utiliser_k1 << '\t' << utiliser_k2 << endl;
    }

    cout << largest_org <<  "     population size: " << total_cell_count << "   number of orgs: " << m_population.size() <<
    "   approx number of cells before death: " << mut::current_cells + mut::new_death_rate << endl;

}



void SaveState(vector<organism> &orgs, long ts)
{
    string s = mut::output_name + "saves";
    const char* dirname = s.c_str();
    mkdir(dirname, 0777);
    ofstream outfile;
    string f_name = s + "/" + mut::output_name + "ss" + to_string(ts) + ".dat";
    outfile.open(f_name, ios::app);
    for (organism i : orgs)
    {
        vector<cell> org = i.get_organism();
        for (cell c : org)
        {
            outfile << c.get_type() << " " << c.get_acatalysis() << " " << c.get_aconversion() << " " << c.get_bcatalysis() << " " << c.get_bconversion() << " " << c.get_dummy() << " ";
        }
        outfile << endl;
    }
    outfile.close();
}

vector<organism> LoadState()
{
    ifstream infile(mut::SaveName);
    vector<organism> popl;
    string line;
    if (infile.is_open())
    {
        while (getline(infile, line))
        {
            stringstream ss(line);
            vector<string> strings;
            string str;
            while (ss >> str)
            {
                strings.push_back(str);
            }
            vector<cell> cells;
            int x = strings.size();
            if (x % 6 > 0)
                cout << "error with state!" << endl;
            // int cells = x / 6;
            int i=0;
            while (i < x)
            {   
                char type = strings[i][0];
                double kp = stod(strings[i+1]);
                double dp = stod(strings[i+2]);
                double kq = stod(strings[i+3]);
                double dq = stod(strings[i+4]);
                double dummy = stod(strings[i+5]);

                cell c;
                c.set_cell(type, 1-kp, kp, dp, 1-kq, kq, dq, dummy);
                cells.push_back(c);
                i += 6;
            }
            organism org;
            org.set_organism(cells);
            uint64_t s = mersenne();
            org.seed_org(s);
            popl.push_back(org);
        }
    }
    else
    {
        cout << "No READ-IN FILE IN DIRECTORY" << endl;
        abort();

    }
    return popl;

}



///contains the population of organisms within a vector
class population
{
    vector<organism> m_population;

public:

    population()
    {
        m_population.reserve( mut::capacity / (mut::threshold / 4) );
    }

    void set_population(vector<organism>& total)
    {
        m_population = total;
        if (mut::invasion)
        {
            mut::out_freq = 10;
        }
    }
    

    vector<organism>& get_population()
    {
        return m_population;
    }

    void process()
    {
        long time_steps = mut::time_steps;
        ///loop for total time steps
        for (int steps = 0; steps <= time_steps; ++steps)
        {
            ///reset the death calculators which changes each time step
            mut::new_current_cells = 0;
            mut::new_death_rate = 0;

            int x = m_population.size();
            vector<int> check_splits{};
            check_splits.resize(x);
            // iterate through all organisms to find total cells in system
            for (int i = 0; i < x;++i)
            {
                mut::new_current_cells += m_population[i].get_organism_population();
            }

            #pragma omp parallel for num_threads(mut::n_threads)
            for (int i = 0; i < x; ++i)
            {
                check_splits[i] = m_population[i].growth();
            }
            
            // auto start = chrono::high_resolution_clock::now();
            // auto stop = chrono::high_resolution_clock::now();
            // auto duration = chrono::duration_cast<std::chrono::milliseconds>(stop-start);
            // cout << "growth death loop: " << duration.count() << endl;

            for (int i = 0; i < x; ++i)
            {
                mut::new_death_rate += m_population[i].get_rN();
                mut::new_current_cells += m_population[i].get_rN();
            }            
            for (int i = 0; i < x; ++i)
            {
                ///if check split returns 1, split the organism and push the new one to the back of m_popuation
                ///if organism size is 1 or less, erase the organism from the vector
                switch (check_splits[i])
                {
                    case 2:
                    {
                        // m_population.erase(m_population.begin() + i);
                        // efficiently kill the organism. 
                        // auto it = m_population.begin() + i;
                        // *it = move(m_population.back());
                        // m_population.pop_back();

                        m_population.erase(m_population.begin() + i);

                        check_splits.erase(check_splits.begin() + i);
                        --x;
                        --i;

                        break;
                    }
                    case 1:
                    {
                        vector<cell> new_vector = m_population.at(i).split();
                        organism new_organism;
                        new_organism.set_organism(new_vector);
                        uint64_t rand = mersenne();
                        new_organism.seed_org(rand);
                        ///push new organisms to the back
                        m_population.push_back(new_organism);
                        break;
                        

                    }
                }
            }
            // moran process for case with no selection. 
            if (mut::new_current_cells > (mut::capacity * 0.8) && mut::moran_on == true)
            {
                uniform_int_distribution<int> moran_choice(0, m_population.size() -1);
                int num = moran_choice(mersenne);
                // efficiently kill the organism. 
                auto it = m_population.begin() + num;
                *it = move(m_population.back());
                m_population.pop_back();

                auto new_it = check_splits.begin() + num;
                *new_it = move(check_splits.back());
                check_splits.pop_back();

            }

            //mut::new_death_rate is equivalent to rN, therefore rN * N / K. 
            mut::death_rate = (mut::new_death_rate * mut::new_current_cells) / mut::capacity;
            ///calculate next time step death rate ((<r> * N^2) / K)
            mut::current_cells = mut::new_current_cells;
            

            double simtime = double(steps) * mut::delta_factor;
            if (floor(simtime) == simtime)
            {


                if (static_cast<int>(simtime) % mut::out_freq == 0)
                {
                    get_averages(m_population, simtime);
                }

                if (mut::savestates && static_cast<int>(simtime) % mut::savet == 0)
                {
                    SaveState(m_population, simtime);
                }
            }
        }

    }

    // void init_variance()
    // {
    //     int x = m_population.size();
    //     for (int i=0;i<x;++i)
    //     {
    //         vector<cell>& org = m_population[i].get_organism();
    //         int n = org.size();
    //         for (int j=0;j<n;++j)
    //         {
    //             org[j].spread();
    //         }
    //     }
    // }

};




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

    vector<cell> fit_vector;
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
    vector<cell> suppress_vector;
    suppress_vector.resize(mut::start_size);
    for (unsigned int i = 0; i < mut::start_size; ++i)
    {
        ///ensure unequal frequencies for max fitness
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
    vector<organism> org_list;


    int n_orgs = static_cast<int>(ceil((mut::capacity * 0.85) / mut::start_size));

    vector<organism> total;
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




int main(int argc, char *argv[])
{
    if (argc > 3)
    {
        mut::threshold = stoi(argv[1]);
        mut::mut_rate = stod(argv[2]);
        mut::a_pg = stod(argv[3]);
        mut::a_switch = stod(argv[4]);
        mut::b_pg = stod(argv[5]);
        mut::b_switch = stod(argv[6]);
        string name;
        for (int i=1;i<argc;++i)
        {
            name += argv[i];
            if (i < argc - 1)
                name += "-";
        }
        mut::output_name = name;
        mut::start_size = mut::threshold / 2;
        mut::capacity = mut::array_tc * mut::threshold;
        mut::choose_start == true;
    }
    else if (argc > 2)
    {
        if (!mut::alpha_sweep)
        {
            cout << argv[1] << "  " << argv[2] << endl;
            mut::threshold = stoi(argv[1]);
            mut::mut_rate = stod(argv[2]);
            cout << mut::threshold << "  " << mut::mut_rate << endl;
            string name;
            for (int i=1;i<argc;++i)
            {
                name += argv[i];
                if (i < argc - 1)
                    name += "-";
            }
            cout << name << endl;    
            mut::output_name = name;
            mut::start_size = mut::threshold / 2;
            mut::capacity = mut::array_tc * mut::threshold;
        }
        else
        {
            cout << argv[1] << "  " << argv[2] << endl;
            mut::threshold = stoi(argv[1]);
            mut::alpha = stod(argv[2]);
            mut::mut_rate = 0.01;
            cout << mut::threshold << "  " << mut::alpha << endl;
            string name;
            for (int i=1;i<argc;++i)
            {
                name += argv[i];
                if (i < argc - 1)
                    name += "-";
            }
            cout << name << endl;    
            mut::output_name = name;
            mut::start_size = mut::threshold / 2;
            mut::capacity = mut::array_tc * mut::threshold;            
        }

    }
    else if (argc > 1)
    {
        mut::threshold = stoi(argv[1]);
        mut::start_size = mut::threshold / 2;
        mut::capacity = mut::array_tc * mut::threshold;
        mut::output_name = argv[1];
    } 


    

    population environment;
    if (mut::invasion)
    {
        // CURRENTLY DEPRECATED AS ORGANISMS ARE NOT SEEDED
        environment = invasion();
    }
    else if (mut::loadstate)
    {
        vector<organism> pop = LoadState();
        environment.set_population(pop);
    }
    else
    {
        cell one;
        cell two;
        if (mut::choose_start == true)
        {
            ///DO NOT PUT VALUES IN HERE!!!! This is to test if equilibrium is stable
            one.set_cell('A', 1-mut::a_pg, mut::a_pg, mut::a_switch, 1 - mut::b_pg, mut::b_pg, mut::b_switch, 0.5);
            two.set_cell('B', 1-mut::a_pg, mut::a_pg, mut::a_switch, 1 - mut::b_pg, mut::b_pg, mut::b_switch, 0.5);
        }
        else
        {
            one.set_cell('A', 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5);
            two.set_cell('B', 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5);
        }

        organism new_organism;
        vector<cell> new_vector;
        new_vector.resize(mut::start_size);
        if (mut::start_broken_frequencies == true)
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
            organism cp = new_organism;
            uint64_t rand = mersenne();
            cp.seed_org(rand);
            total.at(i) = cp;
        }
        environment.set_population(total);
        
    }
    // if (mut::start_broken_frequencies)
    //     environment.init_variance();
    
    environment.process();
        
    return 0;
}


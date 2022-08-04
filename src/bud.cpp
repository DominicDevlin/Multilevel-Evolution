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

using namespace std;

namespace mut
{
    constexpr int start_size{100};
    constexpr int start_pop{100};

    /// 1 / mut_rate is the chane of mutation
    constexpr double mut_rate{100};


    ///threshold before an organisms begins spitting
    int threshold{10000};

    ///total cells in the system
    long total_cells{threshold * 80};

    constexpr long time_steps{100000000};

    ///death rate is calculated based on current cells from the previous time step and an average of the "r" values
    double death_rate{};
    double new_death_rate{};
    long current_cells{};
    long new_current_cells{};

    ///the time step divisor to approximate a continuous model
    constexpr double delta_factor{0.4};


    ///standard deviation for mutations (deprecated)
    constexpr double funct_var{0.05};
    constexpr double conv_var{0.05};

    ///tallies total propagules in system
    long propagules{};

    ///data files
    bool output_values{true};
    bool output_var{true};
    bool output_rel{true};


    bool start_broken{false};
    constexpr double a_pg{0};
    constexpr double a_switch{0.5};
    constexpr double b_pg{0.8};
    constexpr double b_switch{0.00005};   

    // what values to check for keys
    constexpr int min_val{60};
    constexpr int max_val{80};

    //max propagule size
    constexpr int max_prop_size{50};
}

///RNG
auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
mt19937 mersenne( static_cast<mt19937::result_type>(seed) );

///Global distributions so they don't have to be remade every time they are called
uniform_real_distribution<double> double_num(0.0, 1.0);
uniform_int_distribution<> mut_dist(1, mut::mut_rate);
uniform_int_distribution<> splitter(0, 1);
normal_distribution<double> curve(0, 0.02);


///integer RNG
int generate_random_int(const int min, const int max)
{
    uniform_int_distribution<> num(min, max);
    int prob = num(mersenne);
    return prob;
}


void normal_mutation(double gene)
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
}


///mutate the chosen genotype log normal, accounting for boundaries
void mutation(double &gene)
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
}

class cell
{
    /// cell type can either be 'A' or 'B', each has 3 values (replication = 'r', catalysis = '1 - r', conversion = 'c'
    char m_type;
    double ma_replication;
    double ma_catalysis;
    double ma_conversion;
    double mb_replication;
    double mb_catalysis;
    double mb_conversion;
    double m_dummy;
    vector<double> codes;

public:
    void set_cell(char type, double replication1, double catalysis1, double conversion1, double replication2, double catalysis2, double conversion2, double dummy)
    {
        m_type = type;
        ma_replication = replication1;
        ma_conversion = conversion1;
        ma_catalysis = catalysis1;
        mb_replication = replication2;
        mb_conversion = conversion2;
        mb_catalysis = catalysis2;
        m_dummy = dummy;
    }
    void set_code_vector()
    {
        codes.assign(80, 0.5);
    }

    void set_code()     ///set a random <double> code at the front of the cell, and remove the one at the back. Used to determine close relatives.
    {
        double code = double_num(mersenne);
        codes.insert(codes.begin(), code);
        codes.pop_back();
    }

    vector<double>& get_codes()
    {
        return codes;
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

    void mutate()
    {
        ///if mutation is hit, mutate whole genotype.
        double prob = mut_dist(mersenne);
        if (prob == 1)
        {
            mutation(ma_replication);
            ma_catalysis = 1.0 - ma_replication;

            mutation(ma_conversion);

            mutation(mb_replication);
            mb_catalysis = 1.0 - mb_replication;

            mutation(mb_conversion);

            normal_mutation(m_dummy);
        }
    }
};

///calculate poisson distributed death rate based on the previous generation cell average, adapted to the size of the organism
int death_calculator(int organism_size)
{
    double death_calc = mut::death_rate * (static_cast<double>(organism_size) / static_cast<double>(mut::current_cells));
    poisson_distribution<int> p_distribution(death_calc);
    int deaths = p_distribution(mersenne);
    return deaths;
}

///each organism has its cell contained within 'm_current' and its propagules contained within 'propagules', which is a vector of organisms defined within the organism..
class organism
{
    vector<cell> m_current;
    vector<organism> propagules;

public:

    // organism()
    // {
    //     m_current.reserve(mut::threshold*1.1);
    // }

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


    int growth()
    {
        mut::new_current_cells += m_current.size();
        ///this will become the (1-<r>) used to calculate organism growth rate

        ///guess how many cells will be over the threshold to calculate number of propagules
        int guess = death_calculator(m_current.size());
        int prop_calc = guess + mut::threshold;

        ///poisson distribution to determine cellular replication, 'rate' = (r/<r> * (1-<r>) * delta_t)
        
        ///iterate through all cells and add their (1-r), then average
        double functionality{};
        double replication{};
        ///create vector 'fitness' for discrete disribution
        vector<double> fitness;
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
                    cout << "Switch Error" << endl;
                    break;
            }
        }

        functionality = functionality / m_current.size();
        replication = replication / m_current.size();

        discrete_distribution<int> total_fitness(fitness.begin(), fitness.end());
        double divisions{};
        divisions = m_current.size() * functionality * replication * mut::delta_factor;
        poisson_distribution<int> p_distribution(divisions);
        int replicants = p_distribution(mersenne);
        mut::new_death_rate += replicants;

        for (int i = 0; i < replicants; ++i)
        {
            int number = total_fitness(mersenne);
            if (static_cast<int>(m_current.size()) > prop_calc)
            {
                ///if the size of the organism is greater than the threshold that is predicted after deaths, cells replicate into propagules instead.
                ///create a new organism and push it onto the propagule vector, mutate and convert back cell
                vector<cell> new_organism;
                new_organism.push_back(m_current.at(number));
                new_organism.push_back(m_current.at(number));
                double conv = double_num(mersenne);
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
                        cout << "Switch Error" << endl;
                        break;
                }
                new_organism.back().mutate();
                
                // erase the cell from the host organism as it moves into a propagule

                // replace with back, pop back
                auto it = m_current.begin() + number;
                *it = move(m_current.back());
                m_current.pop_back();
                // m_current.erase(m_current.begin() + number);
            
                // auto fit_it = fitness.begin() + number; 
                // *fit_it = move(fitness.back());
                // fitness.pop_back();
                //fitness.erase(fitness.begin() + number);

                vector<double> checks = new_organism.at(0).get_codes();
                int sizer = m_current.size();   ///iterate through organism to find related individuals within 60-80 generations. THis number may be too large.
                // iterate through cells in organism backwards. 
                for (int x = (sizer - 1); x > -1; --x)
                {

                    vector<double> codes = m_current.at(x).get_codes();
                    for (int z = mut::min_val; z < mut::max_val; ++z)
                    {
                        for (int a = mut::min_val; a < mut::max_val; ++a)
                        {
                            if (codes[a] == checks[z])   ///look for identical codes.
                            {
                                // push back cell if match onto propagule
                                new_organism.push_back(m_current.at(x));
                                // m_current.erase(m_current.begin() + x);
                                // efficient removal from mother
                                auto remover = m_current.begin() + x;
                                *remover = move(m_current.back());  
                                m_current.pop_back();
                                // end code loops
                                a = mut::max_val;
                                z = mut::max_val;
                            }
                            else
                                continue;
                        }
                    }
                    if (new_organism.size() > mut::max_prop_size)
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
                discrete_distribution<int> redo(fitness.begin(), fitness.end());
                total_fitness = redo;
                organism prop;
                ///set the new organism with x related cells
                prop.set_organism(new_organism);
                propagules.push_back(prop);
            }
            else    ///if the organism is not dividing into propagules
            {
                m_current.at(number).set_code();
                ///copy the cell and push it to the back of the vector, then mutate last.
                m_current.push_back(m_current.at(number));
                m_current.back().set_code();
                //m_current.at(number).mutate();

                ///convert the cell state of the push back cell at its conversion probability.
                double conv = double_num(mersenne);
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
            int number = generate_random_int(0, m_current.size() - 1);
            // m_current.erase(m_current.begin() + number);
            auto killed = m_current.begin() + number;
            *killed = move(m_current.back());
            m_current.pop_back();
        }

        if (propagules.size() > 0)
        {
            return propagules.size();
        }
        else
            return 0;

    }

    ///spit the propagules from the host organism.
    vector<cell> split()
    {
        int x = propagules.size();
        for (int i = 0; i < x; ++ i)
        {
            vector<cell> org = propagules.at(i).get_organism();
            propagules.erase(propagules.begin() + i);
            return org;
        }
        cout << "SPITTING ERROR" << endl;
        vector <cell> error;
        return error;
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




///calculate averages to send to output
void get_averages(vector<organism>& m_population, long generations)
{
    double a_repl_average{};
    double a_catal_average{};
    double b_repl_average{};
    double b_catal_average{};
    double dummy_average{};

    double a_total_average{};
    double b_total_average{};
    double a_conv_average{};
    double b_conv_average{};

    double c_total_average{};
    long total_cell_count{};

    int total_organisms{static_cast<int>(m_population.size())};
    long largest_org = {};
    for (unsigned int i = 0; i < m_population.size(); ++i)
    {
        vector<cell> test = m_population.at(i).get_organism();

        if (static_cast<int>(test.size()) > largest_org)
            largest_org = test.size();

        for (cell &c : test)
        {
            switch (c.get_type())
            {
                case 'A':
                {
                    ++ a_total_average;
                    ++ total_cell_count;
                    break;
                }
                case 'B':
                {
                    ++ b_total_average;
                    ++ total_cell_count;
                    break;
                }
                case 'C':
                {
                    ++ c_total_average;
                    ++ total_cell_count;
                    break;
                }
            }
            a_repl_average += c.get_areplication();
            a_catal_average += c.get_acatalysis();
            a_conv_average += c.get_aconversion();
            b_repl_average += c.get_breplication();
            b_catal_average += c.get_bcatalysis();
            b_conv_average += c.get_bconversion();
            dummy_average += c.get_dummy();
        }
    }
    a_repl_average = a_repl_average / total_cell_count;
    a_catal_average = a_catal_average / total_cell_count;
    b_repl_average = b_repl_average / total_cell_count;
    b_catal_average = b_catal_average / total_cell_count;
    a_conv_average = a_conv_average / total_cell_count;
    b_conv_average = b_conv_average / total_cell_count;
    a_total_average = a_total_average / total_organisms;
    b_total_average = b_total_average / total_organisms;
    c_total_average = c_total_average / total_organisms;
    dummy_average = dummy_average / total_cell_count;



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
        string file_name = to_string(mut::threshold) + ".dat";
        ofstream outfile;
        outfile.open(file_name, ios::app);
        outfile << generations << '\t' << a_total_average / (a_total_average + b_total_average) << '\t' << a_total_average << '\t' << a_catal_average << '\t' << a_conv_average << '\t'
        << b_total_average << '\t' << b_catal_average << '\t' << b_conv_average << endl;
        outfile.close();
    }


    if (mut::output_rel == true)
    {
        string var_name = to_string(mut::threshold) + "rel.dat";
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

    cout << largest_org <<  "     population size: " << total_cell_count << "   number of orgs: " << m_population.size() <<
    "   approx number of cells before death: " << mut::current_cells + mut::new_death_rate << endl;
}


///contains the population of organisms within a vector
class population
{
    vector<organism> m_population;

public:

    void set_population(vector<organism>& total)
    {
        m_population = total;
    }

    vector<organism> get_population()
    {
        return m_population;
    }

    void process(long time_steps)
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
                int check_split = m_population.at(i).growth();
                ///if -1 is returned, erase organism, if more than 1 is returned, spit the propagules which have reached greater than 10 in size
                switch (check_split)
                {
                    case -1:
                    {
                        // m_population.at(i) = m_population.at(x - 1);
                        // m_population.erase(m_population.begin() + (x - 1));
                        // remove organism. 
                        auto it = m_population.begin() + i;
                        *it = move(m_population.back());
                        m_population.pop_back();
                        -- x;
                        -- i;
                        break;
                    }
                    case 0:
                        break;
                    default:
                        ///split propagules, iterate for how many in the organism have reached 10+
                        for (int y = 0; y < check_split; ++y)
                        {
                            vector<cell> new_vector = m_population.at(i).split();
                            organism new_organism;
                            new_organism.set_organism(new_vector);
                            ///push new organisms to the back
                            m_population.push_back(new_organism);
                            // cout << "Propagule has " << new_vector.size() << " cells" << endl;
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
};

int main()
{
    cell one;
    cell two;
    one.set_code_vector();
    two.set_code_vector();


    if (mut::start_broken == true)
    {
        // Start using specific vals
        one.set_cell('A', 1-mut::a_pg, mut::a_pg, mut::a_switch, 1 - mut::b_pg, mut::b_pg, mut::b_switch, 0.5);
        two.set_cell('B', 1-mut::a_pg, mut::a_pg, mut::a_switch, 1 - mut::b_pg, mut::b_pg, mut::b_switch, 0.5);
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
    environment.process(mut::time_steps);

    return 0;
}

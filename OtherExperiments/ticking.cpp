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

namespace mut
{


    /// 1 / mut_rate is the chance of mutation
    constexpr double mut_rate{0.05};


    ///the maximum number of cells in an organism before it splits
    constexpr int threshold{422};

    //starting parameters
    constexpr int start_size{mut::threshold / 2};
    constexpr int start_pop{150};

    ///total cells in the system
    long total_cells{threshold * 100};

    ///run steps
    constexpr long time_steps{100000000};

    ///returns on cooperation
    double alpha{1.0};

    double trait_e{1.0};

    ///the time step divisor to approximate a continuous model
    constexpr double delta_factor{0.4};


    ///standard deviation for mutations when using standard dist
    constexpr double funct_var{0.03};
    constexpr double conv_var{0.03};


    ///maximum organisms in population - moran process, turn on or off
    constexpr bool moran_on{0};
    constexpr int pop_size{200};



    ///death rate is calculated based on current cells from the previous time step and an average of the "r" values
    double death_rate{};
    double new_death_rate{};
    long current_cells{};
    long new_current_cells{};

    bool output_var{true};
    bool output_values{true};
    bool start_broken{false};


}

///RNG
auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
std::mt19937 mersenne( static_cast<std::mt19937::result_type>(seed) );

///Global distributions so they don't have to be remade every time they are called
std::uniform_real_distribution<double> double_num(0.0, 1.0);
std::uniform_int_distribution<> splitter(0, 1);
std::normal_distribution<double> curve(0, 0.02);

///integer RNG
long generate_random_int(const int min, const int max)
{
    std::uniform_int_distribution<> num(min, max);
    long prob = num(mersenne);
    return prob;
}

///mutate the chosen genotype, accounting for boundaries
double flat_mutation(double gene, double variation)
{
    std::uniform_real_distribution<double> curve(gene - (variation / 2), gene + (variation / 2));
    gene = curve(mersenne);
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

///mutate the chosen genotype, accounting for boundaries
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
    std::poisson_distribution<int> p_distribution(death_calc); ///poisson distribution for organism deaths
    int deaths = p_distribution(mersenne);
    return deaths;
}

int death_calculator_deterministic(unsigned int organism_size)
{
    int deaths = ceil( mut::death_rate * (static_cast<double>(organism_size) / static_cast<double>(mut::current_cells)));  ///calculate death rate based on the previous generation cels, adapted to the size of the organism
    return deaths;
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
        ma_replication = mutation(ma_replication);
        ma_catalysis = 1.0 - ma_replication;

        ma_conversion = mutation(ma_conversion);

        mb_replication = mutation(mb_replication);
        mb_catalysis = 1.0 - mb_replication;

        mb_conversion = mutation(mb_conversion);

        m_dummy = mutation(m_dummy);

    }
};



///cells in organism are contained within 'm_current'
class organism
{
    std::vector<cell> m_current;
    /// check how many mutations are happening per life cycle
    long mut_num{};

public:

    void set_organism(std::vector<cell>& total)
    {
        m_current = total;
    }

    std::vector<cell>& get_organism()
    {
        return m_current;
    }

    int get_organism_population()
    {
        return static_cast<int>(m_current.size());
    }

    int growth()
    {
        ///this is the (1-<r>) used to calculate organism growth rate
        double functionality{};
        double replication{};
        ///iterate through all cells and add their (1-r), then average
        for (cell &i : m_current)
        {
            switch (i.get_type())
            {
                case 'A':
                    functionality += pow(i.get_acatalysis(), mut::alpha);
                    replication += i.get_areplication();
                    ///add to the number of total cells in the system to calculate the next generation death rate
                    ++ mut::new_current_cells;
                    break;
                case 'B':
                    functionality += pow(i.get_bcatalysis(), mut::alpha);
                    replication += i.get_breplication();
                    ++ mut::new_current_cells;
                    break;
                default:
                    std::cout << "Switch Error" << std::endl;
                    break;
            }
        }
        functionality = functionality / m_current.size();
        replication = replication / m_current.size();

        ///iterate through all cells to determine cellular replication 'rate' = (r/<r> * (1-<r>) * N)


        ///create vector to be used for determining which cells will replicate
        std::vector<double> fitness;
        fitness.resize(m_current.size());
        for (unsigned int i = 0; i < m_current.size(); ++i)
        {
            switch (m_current.at(i).get_type())
            {
                case 'A':
                    fitness.at(i) = m_current.at(i).get_areplication();
                    break;
                case 'B':
                    fitness.at(i) = m_current.at(i).get_breplication();
                    break;
            }
        }
        ///discrete distribution based on vector
        double divisions{};
        std::discrete_distribution<int> total_fitness(fitness.begin(), fitness.end());

        divisions = m_current.size() * (1 - mut::trait_e + mut::trait_e * functionality) * replication * mut::delta_factor;
        mut::new_death_rate += divisions;
        std::poisson_distribution<int> p_distribution(divisions);
        int replicants = p_distribution(mersenne);

        for (int i = 0; i < replicants; ++i)
        {
            int number = total_fitness(mersenne);
            m_current.push_back(m_current.at(number));
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
                    std::cout << "Switch Error" << std::endl;
                    break;
            }
        }

        ///mutations
        int n_mutations = ceil(mut::mut_rate * mut::delta_factor * m_current.size());
        for (int i = 0; i<n_mutations; ++i)
        {
            int n = generate_random_int(0, m_current.size()-1);
            m_current.at(i).mutate();
        }




        ///proportionate killing of cells in organism

        int deaths = death_calculator_deterministic(m_current.size());
        for (int i = 0; i < deaths; ++i)
        {
            ///if the organism size is 1, kill the organism.
            if (m_current.size() < 3)
            {
                return 2;
            }
            int number = generate_random_int(0, m_current.size() - 1);
            m_current.erase(m_current.begin() + number);
        }



        if (static_cast<int>(m_current.size()) >= mut::threshold)
        {
            return 1;
        }
        return 0;

    }

    ///split the organism into two. return the new organism
    std::vector<cell> split()
    {
        std::vector<cell> new_organism{};
        std::vector<cell> new_copy{};
        ///randomly allocate each cell to old organism or new organism
        for (cell &i : m_current)
        {
            double number = splitter(mersenne);
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
};


double get_within_variance(std::vector<organism>& m_population, double average, double total_cell_count, int pick)
{
    double variance{};

    for (unsigned int i = 0; i < m_population.size(); ++i)
    {
        std::vector<cell> test = m_population.at(i).get_organism();
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

double get_between_variance(std::vector<organism>& m_population, double mean, double total_cell_count, int pick)
{
    double variance{};

    for (unsigned int i = 0; i < m_population.size(); ++i)
    {
        double x_x{};
        double between_mean{};
        std::vector<cell> test = m_population.at(i).get_organism();
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


double get_total_variance(std::vector<organism>& m_population, double mean, double total_cell_count, int pick)
{
    double variance{};
    double xi_xi{};

    for (unsigned int i = 0; i < m_population.size(); ++i)
    {
        std::vector<cell> test = m_population.at(i).get_organism();
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
void get_averages(std::vector<organism>& m_population, long generations)
{
    double a_repl_average{};
    double a_catal_average{};
    double b_repl_average{};
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
    for (unsigned int i = 0; i < m_population.size(); ++i)
    {
        std::vector<cell> test = m_population.at(i).get_organism();
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
                    break;
                }
                case 'B':
                {
                    ++ b_total_average;
                    break;
                }
                case 'C':
                {
                    ++ c_total_average;
                    break;
                }

            }
            ++ total_cell_count;
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
    dummy_average = dummy_average / total_cell_count;

    a_total_average = a_total_average / total_organisms;
    b_total_average = b_total_average / total_organisms;
    c_total_average = c_total_average / total_organisms;

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

    ///std::cout << "         AVERAGES      \n\n"
    ///<< "Total cell count   Total Organsims   MaxOrg Size  " << std::endl;
    ///std::cout << mut::current_cells << "               " << total_organisms << "                " << largest_org << "          "
    ///<< '\n' << " within dummy   between dummy    total dummy\n" << within_dummy << "      " << between_dummy << "       " << full_dummy << std::endl;



    if (mut::output_values == true)
    {
        std::string file_name = std::to_string(mut::threshold) + ".dat";
        std::ofstream outfile;
        outfile.open(file_name, std::ios::app);
        outfile << generations << "\t" << a_total_average << "\t" << a_repl_average << '\t' << a_catal_average << '\t' << a_conv_average << '\t'
        << b_total_average << '\t' << b_repl_average << '\t' << b_catal_average << '\t' << b_conv_average << std::endl;
        outfile.close();
    }


    if (mut::output_var == true)
    {
        std::string var_name = std::to_string(mut::threshold) + "var.dat";
        std::ofstream outfile;
        outfile.open(var_name, std::ios::app);
        outfile << global_average_r << "\t" << relatedness << "\t" << within_dummy << '\t' << between_dummy << '\t' << k1_relatedness << "\t" << d1_relatedness 
		<< "\t" << k2_relatedness << "\t" << d2_relatedness << std::endl;
        outfile.close();
    }


}


///contains the population of organisms within a vector
class population
{
    std::vector<organism> m_population;

public:

    void set_population(std::vector<organism>& total)
    {
        m_population = total;
    }

    std::vector<organism> get_population()
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
            int x = m_population.size();
            for (int i = 0; i < x; ++i)
            {
                ///if check split returns 1, split the organism and push the new one to the back of m_popuation
                ///if organism size is 5 or less, erase the organism from the vector
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
                        std::vector<cell> new_vector = m_population.at(i).split();
                        organism new_organism;
                        new_organism.set_organism(new_vector);
                        ///push new organisms to the back
                        m_population.push_back(new_organism);
                        break;
                    }
                }
            }
            /// if total organisms is greater than max population size (200), randomly kill organisms (moran process), TURN ON AND OFF
//            if (mut::moran_on == true)
//                int kill = m_population.size() - mut::pop_size;
//                while (kill > 0)
//                {
//                    int num = generate_random_int(0, m_population.size() - 1);
//                    m_population.erase(m_population.begin() + num);
//                    -- kill;
//                }

            ///calculate next time step death rate ((<r> * N^2) / K)
            mut::current_cells = mut::new_current_cells;
            mut::new_death_rate = mut::new_death_rate / mut::current_cells;

            mut::death_rate =  (mut::new_death_rate * pow(mut::current_cells, 2)) / mut::total_cells;

            // if (time_steps % 200 == 0)
            // {
            //     std::cout << total_steps - time_steps << " time steps have passed" << std::endl;
            // }
            
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
    ///set cells to initial values so that states 'A' and 'B' are equal, commented out starting from a symmetry broken state (change starting values here)
//    two.set_cell('A', 0.8, 1.0, 0.0, 0.0, 0.0, 1.0);
//    one.set_cell('B', 0.8, 1.0, 0.0, 0.0, 0.0, 1.0);
    if (mut::start_broken == true)
    {
        ///DO NOT PUT VALUES IN HERE!!!! This is to test if equilibrium is stable
        one.set_cell('A', 0.99, 0.01, 0.42, 0.4, 0.6, 0.0005, 0.5);
        two.set_cell('B', 0.99, 0.01, 0.42, 0.4, 0.6, 0.0005, 0.5);
    }
    else
    {
        one.set_cell('A', 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5);
        two.set_cell('B', 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5);
    }


    organism new_organism;
    std::vector<cell> new_vector;
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

    std::vector<organism> total;
    total.resize(mut::start_pop);
    for (unsigned int i = 0; i < total.size(); ++i)
    {
        total.at(i) = new_organism;
    }
    population environment;


    environment.set_population(total);
    environment.process(mut::time_steps);



    ///send start info to data file
//    std::ofstream outfile;
//    outfile.open("split_data.dat", std::ios::app);
//    outfile << "Mutation rate: " << mut::mut_rate << std::endl;
//    outfile << "Mutations occur by Normal Distribution. Variance of repl/catal mutation: " << mut::funct_var << "   Variance of conversion mutation: " << mut::conv_var << std::endl;
//    outfile << "Max Cell Number " << mut::total_cells << "   Max Organism Size: " << mut::threshold << "    Total time steps ran for: " << mut::time_steps << std::endl;
//    outfile.close();

    return 0;
}

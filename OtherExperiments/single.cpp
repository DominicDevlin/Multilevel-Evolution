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
    constexpr double mut_rate{178};
    ///the maximum number of cells in an organism before it splits
    constexpr int threshold{300};


    constexpr int start_size{100};
    constexpr int start_pop{100};
    

    constexpr long time_steps{600000};

    

    ///total cells in the system
    long total_cells{threshold * 100};


    ///standard deviation for mutations
    constexpr double funct_var{0.05};
    constexpr double conv_var{0.05};


    ///maximum organisms in population - moran process, turn on or off
    constexpr bool moran_on{0};
    constexpr int pop_size{200};

    ///the time step divisor to approximate a continuous model
    constexpr double delta_factor{0.4};

    ///death rate is calculated based on current cells from the previous time step and an average of the "r" values
    double death_rate{};
    double new_death_rate{};
    long current_cells{};
    long new_current_cells{};


    ///turn write to file on or off
    bool output_var{true};
}

///RNG
auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
std::mt19937 mersenne( static_cast<std::mt19937::result_type>(seed) );

///Global distributions so they don't have to be remade every time they are called
std::uniform_real_distribution<double> double_num(0.0, 1.0);
std::uniform_int_distribution<> mut_dist(1, mut::mut_rate);
std::uniform_int_distribution<> splitter(0, 1);
std::normal_distribution<double> curve(0, 0.02);
//std::uniform_int_distribution<> death_chance(1, 50 / mut::delta_factor);


///integer RNG
long generate_random_int(const int min, const int max)
{
    std::uniform_int_distribution<> num(min, max);
    long prob = num(mersenne);
    return prob;
}


///mutate the chosen genotype, accounting for boundaries
double mutation(double gene)
{
    //std::uniform_real_distribution<double> curve(gene - (variation / 2), gene + (variation / 2));
    gene = gene * exp(-curve(mersenne));
    //gene = curve(mersenne);
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



///calculate poisson distributed death rate based on the previous generation cell average, adapted to the size of the organism
int death_calculator(unsigned int organism_size)
{
    double death_calc = mut::death_rate * (static_cast<double>(organism_size) / static_cast<double>(mut::current_cells));  //calculate death rate based on the previous generation cels, adapted to the size of the organism
    std::poisson_distribution<int> p_distribution(death_calc); //poisson distributino for organism deaths
    int deaths = p_distribution(mersenne);
    return deaths;
}

class cell
{
    /// cell type can either be 'A' or 'B', each has 3 values (replication = 'r', catalysis = '1 - r', conversion = 'c'

    double m_replication;
    double m_catalysis;
    double m_dummy;


public:
    void set_cell(double replication1, double catalysis1, double dummy)
    {

        m_replication = replication1;
        m_catalysis = catalysis1;
        m_dummy = dummy;
    }

    double get_dummy()
    {
        return m_dummy;
    }

    double get_replication()
    {
        return m_replication;
    }

    double get_catalysis()
    {
        return m_catalysis;
    }


    void mutate()
    {
        ///if mutation is hit, mutate whole genotype.
        double prob = mut_dist(mersenne);
        if (prob == 1)
        {
            m_replication = mutation(m_replication);
            m_catalysis = 1.0 - m_replication;
            m_dummy = normal_mutation(m_dummy);
        }
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
            functionality += i.get_catalysis();
            replication += i.get_replication();
            ///add to the number of total cells in the system to calculate the next generation death rate
            ++ mut::new_current_cells;
        }
        functionality = functionality / m_current.size();
        replication = replication / m_current.size();

        ///iterate through all cells to determine cellular replication 'rate' = (r/<r> * (1-<r>) * N)


        ///create vector to be used for determining which cells will replicate based on fitness
        std::vector<double> fitness;
        fitness.resize(m_current.size());
        for (unsigned int i = 0; i < m_current.size(); ++i)
        {
            fitness.at(i) = m_current.at(i).get_replication();
        }
        ///discrete distribution based on vector
        std::discrete_distribution<int> total_fitness(fitness.begin(), fitness.end());

        double divisions{m_current.size() * functionality * replication * mut::delta_factor};
        mut::new_death_rate += divisions;
        std::poisson_distribution<int> p_distribution(divisions);
        int replicants = p_distribution(mersenne);
        for (int i = 0; i < replicants; ++i)
        {
            int number = total_fitness(mersenne);
            m_current.push_back(m_current.at(number));
            m_current.at(number).mutate();
            m_current.back().mutate();
            ///push back cell and mutate.
        }
        ///proportionate killing of cells in organism
        int deaths = death_calculator(m_current.size());
        for (int i = 0; i < deaths; ++i)
        {
            ///if the organism size is 5 or less, kill the organism.
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


double get_full_variance(std::vector<organism>& population, double total_cell_count, double mean, int pick)
{
    double variance{};
    double x_x{};

    for (unsigned int i = 0; i < population.size(); ++i)
    {
        std::vector<cell> test = population.at(i).get_organism();
        for (cell &c : test)
        {
            switch (pick)
            {
            case 1:
                x_x += pow(c.get_catalysis() - mean, 2);
                break;
            case 2:
                x_x += pow(c.get_dummy() - mean, 2);
                break;
            }
        }
    }
    variance = x_x / total_cell_count;
    return variance;
}

double get_within_variance(std::vector<organism>& population, double total_cell_count, double average, int pick)
{
    double variance{};

    for (unsigned int i = 0; i < population.size(); ++i)
    {
        std::vector<cell> test = population.at(i).get_organism();
        double within_mean{};
        double x_x{};
        for (cell &c : test)
        {
            switch (pick)
            {
            case 1:
                within_mean += c.get_catalysis();
                break;
            case 2:
                within_mean += c.get_dummy();
                break;
            }
        }
        within_mean = within_mean / test.size();
        for (cell &c : test)
        {
            switch (pick)
            {
            case 1:
                x_x += pow(c.get_catalysis() - within_mean, 2);
                break;
            case 2:
                x_x += pow(c.get_dummy() - within_mean, 2);
                break;
            }
        }
        variance += x_x;
    }
    variance = variance / total_cell_count;
    return variance;
}

double get_between_variance(std::vector<organism>& population, double total_cell_count, double mean, int pick)
{
    double variance{};
    for (unsigned int i = 0; i < population.size(); ++i)
    {
        double x_x{};
        double between_mean{};
        std::vector<cell> test = population.at(i).get_organism();
        for (cell &c : test)
        {
            switch (pick)
            {
            case 1:
                between_mean += c.get_catalysis();
                break;
            case 2:
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


///calculate averages to send to output. Averages are of all cell values & cell numbers
void get_averages(std::vector<organism>& m_population, long generations)
{
    double repl_average{};
    double catal_average{};
    long total_average{};

    double dummy_average{};

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
            ++ total_cell_count;
            repl_average += c.get_replication();
            catal_average += c.get_catalysis();
            dummy_average += c.get_dummy();
        }
    }
    repl_average = repl_average / total_cell_count;
    catal_average = catal_average / total_cell_count;
    total_average = total_cell_count / total_organisms;
    dummy_average = dummy_average / total_cell_count;



    //double full_repl_variance = get_full_variance(m_population, total_cell_count, catal_average, 1);
    double within_repl_variance = get_within_variance(m_population, total_cell_count, catal_average, 1);
    double between_repl_variance = get_between_variance(m_population, total_cell_count, catal_average, 1);

    //double full_dummy_var = get_full_variance(m_population, total_cell_count, dummy_average, 2);
    double within_dummy_var = get_within_variance(m_population, total_cell_count, dummy_average, 2);
    double between_dummy_var = get_between_variance(m_population, total_cell_count, dummy_average, 2);


    double k_relatedness = between_repl_variance / (between_repl_variance + within_repl_variance);
    double relatedness = between_dummy_var / (between_dummy_var + within_dummy_var);

    double global_average_r = catal_average * (1 - catal_average);



    if (mut::output_var == true)
    {

        std::string var_name = std::to_string(mut::threshold) + "s.dat";
        std::ofstream outfile;
        outfile.open(var_name, std::ios::app);
        outfile << catal_average << "\t" << global_average_r << "\t" << relatedness << "\t" << within_dummy_var 
        << "\t" << between_dummy_var << "\t" << k_relatedness << std::endl;
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
            mut::death_rate = (mut::new_death_rate * pow(mut::current_cells, 2)) / mut::total_cells;
            if (time_steps % 200 == 0)
            {
                std::cout << total_steps - time_steps << " time steps have passed" << std::endl;
            }
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
    one.set_cell(0.5, 0.5, 0.5);
    organism new_organism;
    std::vector<cell> new_vector;
    new_vector.resize(mut::start_size);

    for (unsigned int i = 0; i < new_vector.size(); ++i)
    {
        new_vector.at(i) = one;
    }
    new_organism.set_organism(new_vector);




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

    return 0;
}

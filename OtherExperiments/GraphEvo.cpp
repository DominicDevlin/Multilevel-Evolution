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


namespace vars
{
    constexpr int mut_rate{3};

    constexpr int cells_per_org{7};

    constexpr int population_size{200};

    constexpr int moran_steps{3000};

    constexpr int graph_mut_rate{10};

    constexpr int n_graph_mutations{3};

    constexpr long sim_steps{10000};

    constexpr bool start_moran{true};

    constexpr bool output{false};

    constexpr bool all_mutations_on{true};
}



///RNG
auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
std::mt19937 mersenne( static_cast<std::mt19937::result_type>(seed) );

///Global distributions
std::uniform_real_distribution<double> double_num(0.0, 1.0);
std::uniform_int_distribution<> mut_dist(1, vars::mut_rate);
std::normal_distribution<double> curve(0, 0.02);
std::uniform_int_distribution<> splitter(0, 1);
std::uniform_int_distribution<> new_edge(0, vars::cells_per_org - 1);
std::uniform_int_distribution<> graph_mut_dist(1, vars::graph_mut_rate);



///integer RNG
long generate_random_int(const int min, const int max)
{
    std::uniform_int_distribution<> num(min, max);
    long prob = num(mersenne);
    return prob;
}


double normal_mutation(double gene, double max_val)
{
    gene = gene + curve(mersenne);
    if (gene > max_val)
    {
        gene = max_val;
    }
    if (gene < 0)
    {
        gene = 0;
    }
    return gene;
}



class cell
{
    double m_public_good;
    double m_dummy;

public:
 
    void set_cell(double public_good, double dummy)
    {
        m_public_good = public_good;
        m_dummy = dummy;
    }

    double get_pg()
    {
        return m_public_good;
    }

    double get_dummy()
    {
        return m_dummy;
    }

    void mutate()
    {
        int mut_prob = mut_dist(mersenne);
        if (mut_prob == 1)
        {
            m_public_good = normal_mutation(m_public_good, 1);
            m_dummy = normal_mutation(m_dummy, 100);
        }
    }

    void initialise_cell()
    {
        m_public_good = 0.5;
        m_dummy = 50;
    }


};


/// Each organism holds a cell_n * cell_n matrix that determines the weight of the edges (which cells are connected to which cells)
/// every cell has to have at least one input or output. That is sum i + sum j > 0 for all i = j
class organism
{
    std::vector<std::vector<int>> weights;
    std::vector<cell> cell_list;
public:

    std::vector<cell>& get_cell_list()
    {
        return cell_list;
    }

    std::vector<std::vector<int>>& get_graph()
    {
        return weights;
    }


    double get_fitness()
    {
        double fitness = 0;
        for (cell &c : cell_list)
        {
            fitness += c.get_pg() * (1 - c.get_pg());
        }
        return fitness / vars::cells_per_org;
    }

    double get_d_average()
    {
        double average = 0;
        for (cell &c : cell_list)
        {
            average += c.get_dummy();
        }
        return average / vars::cells_per_org;
    }

    void init_weights()
    {
        if (vars::start_moran)
        {
            weights.resize(vars::cells_per_org);
            for (int i = 0; i < vars::cells_per_org; ++i)
            {
                weights.at(i).resize(vars::cells_per_org);
                for (int j = 0; j < vars::cells_per_org; ++j)
                {
                    if (i == j)
                        weights.at(i).at(j) = -1;
                    else 
                    weights.at(i).at(j) = 1;
                }
                    
            }
        }
        else
        {
            for (int i = 0; i < vars::cells_per_org; ++i)
            {
                
                std::vector<int> random_weights{};
                for (int j = 0; j < vars::cells_per_org; ++j)
                {
                    /// initial spatial arrangement is p=0.5 node from one cell to the next, p=0.5 for no node, except -1 for i=j where cell connects to itself.
                    if (i != j)
                        random_weights.push_back(splitter(mersenne));
                    else
                        random_weights.push_back(-1);
                }
                weights.push_back(random_weights);

            }
            for (int i = 0; i < vars::cells_per_org; ++i)
            {
                for (int j = 0; j < vars::cells_per_org; ++j)
                {
                    std::cout << weights.at(i).at(j) << " ";
                }
                std::cout << std::endl;
            }        
        }
        
        

    }

    void init_cells()
    {
        for (int i = 0; i < vars::cells_per_org; ++i)
        {
            cell new_cell;
            new_cell.initialise_cell();
            cell_list.push_back(new_cell);
        }

    }

    void refresh_cells()
    {
        for (cell &i : cell_list)
        {
            i.initialise_cell();
        }
    }

    void mutate_graph()
    {
        int num = graph_mut_dist(mersenne);
        if (num == 1)
        {
            if (vars::all_mutations_on)
            {
                int input = new_edge(mersenne);
                int output = new_edge(mersenne);
                while (input == output)
                    {
                        output = new_edge(mersenne);
                    }
                weights.at(input).at(output) = splitter(mersenne);
            }
            else
            {
                for (int i = 0; i < vars::n_graph_mutations; ++i)
                {
                    int input = new_edge(mersenne);
                    int output = new_edge(mersenne);
                    while (input == output)
                    {
                        output = new_edge(mersenne);
                    }
                    /// initally i will make the probability of a connection being equally 0 ir 1. 
                    ///This means that, with drift, we expect each cell to be connected to half of the other cells by output alone.
                    ///I should try favouring there being no edge eventually. 
                    int sum_in = 0;
                    int sum_out = 0;
                    for (int i = 0; i < vars::cells_per_org; ++i)
                    {
                        for (int j = 0; j < vars::cells_per_org; ++j)
                        if (i != j)
                        {
                            if (i == input || j == input)
                                sum_in += 1;
                            
                            if (i == output || j == output)
                                sum_out += 1;
                        }
                    }
                    if (sum_in < 2 && sum_out < 2)
                    {
                        if (weights.at(input).at(output))
                            continue;
                        else
                            weights.at(input).at(output) = splitter(mersenne);
                    }
                    else
                    {
                        weights.at(input).at(output) = splitter(mersenne);
                        // std::cout << "mutating" << std::endl;                       
                    }
                }
            }
        }


    }

    /// this is quite inefficient because only do a single reproduction event every step.
    void moran_step()
    {
        std::vector<double> fitness{};
        fitness.resize(cell_list.size());
        double k_average{};
        for (cell &i : cell_list)
        {
            k_average += i.get_pg();
        }
        k_average = k_average / vars::cells_per_org;

        for (unsigned int i = 0; i < cell_list.size(); ++i)
        {
            fitness.at(i) = k_average * (1 - cell_list.at(i).get_pg());
        }
        std::discrete_distribution<int> total_fitness(fitness.begin(), fitness.end());
        int c_num = total_fitness(mersenne);

        std::vector<int>weight_row = weights.at(c_num);
        int total{0};
        for (int i = 0; i < vars::cells_per_org; ++i)
        {
            if (i == c_num)
                continue;
            else
                total += weight_row.at(i);
        }
        ///this is probably slower than initialising them all globally, but will take like 100 lines. Redo if necessary
        /// NEED TO FIX THIS ERROR
        // if (total == 0)
        //     std::cout << "TOTAL IS 0" << std::endl;
        std::uniform_int_distribution<> node_dist(1, total + 1);
        int choose_node = node_dist(mersenne);
        int row_num{0};
        for (int i = 0; i < vars::cells_per_org; ++i)
        {
            if (weight_row.at(i) == 1)
            {
                if (choose_node == 1)
                {
                    row_num = i;
                    break;
                }
                else
                    --choose_node;
            }
        }
        ///get reproducing cell and replace the cell chosen based on node
        ///mutate cell before it reproduces (both cells mutate) -- can change this to see the affect of asymmetry later. 
        cell_list.at(c_num).mutate();
        cell chosen_cell = cell_list.at(c_num);
        cell_list.at(row_num) = chosen_cell;
        // std::cout << c_num << "  " << row_num << std::endl;
        // for (int i = 0; i < vars::cells_per_org; ++i)
        // {
        //     for (int j = 0; j < vars::cells_per_org; ++j)
        //     {
        //         std::cout << weights.at(i).at(j) << " ";
        //     }
        //     std::cout << std::endl;
        // }
    }

};


class population
{
    ///do wright-fisher process after calling moran_process 10^4 times for each organism.
    std::vector<organism> m_population;
public:

    void init_population()
    {
        for (int i = 0; i < vars::population_size; ++i)
        {
            organism new_org;
            new_org.init_cells();
            new_org.init_weights();
            m_population.push_back(new_org);
        }

    }

    void growth_phase()
    {
        for (organism &i : m_population)
        {
            int count = 0;
            while (count < vars::moran_steps)
            {
                i.moran_step();
                ++count;
            }
        }
    }


    void fisher_process()
    {
        std::vector<double> organism_fitness (vars::cells_per_org,0);
        organism_fitness.resize(vars::population_size);
        for (int i = 0; i < vars::population_size; ++i)
        {
            organism_fitness.at(i) = m_population.at(i).get_fitness();
        }
        std::discrete_distribution<int> pop_fitness(organism_fitness.begin(), organism_fitness.end());
        for (int i = 0; i < vars::population_size; ++i)
        {
            int num = pop_fitness(mersenne);
            organism next_org = m_population.at(num);
            m_population.at(i) = next_org;
            m_population.at(i).mutate_graph();
        }

        /// reset organisms (germline vs soma?) for next generation

    }

    void new_gen_reset()
    {
    for (organism &i : m_population)
        {
            i.refresh_cells();
        }
    }

    void output_vals()
    {
        double global_k{};
        double global_d{};
        for (organism &i : m_population)
        {
            double org_average = 0;
            std::vector<cell> c_list = i.get_cell_list();
            for (cell &c : c_list)
            {   
                global_k += c.get_pg();
                global_d += c.get_dummy();
                org_average += c.get_pg();
            }
            // std::cout << org_average / 10 << "\n" << std::endl;
        }
        global_k = global_k / (vars::population_size * vars::cells_per_org);
        global_d = global_d / (vars::population_size * vars::cells_per_org);
        std::cout << global_k << std::endl;
        std::cout << global_d << std::endl;
    }

    void output_graph()
    {
        int most_fit = 0;
        double fit_val = 0;
        for (int i = 0; i < vars::population_size; ++i)
        {
            double z = m_population.at(i).get_fitness();
            if (z > fit_val)
            {
                most_fit = i;
                fit_val = z;
            }
        }
        std::vector<std::vector<int>> graph = m_population.at(most_fit).get_graph();
        for (int i = 0; i < vars::cells_per_org; ++i)
        {
            std::vector<int> c_row = graph.at(i);
            for (int j = 0; j < vars::cells_per_org; ++j)
            {
                std::cout << c_row.at(j) << " ";
            }
            std::cout << std::endl;
        }
    }

    

    void to_file(long sim_steps)
    {

        double global_k{};
        double global_d{};
        for (organism &i : m_population)
        {
            global_k += i.get_fitness();
            global_d += i.get_d_average();
            // std::cout << org_average / 10 << "\n" << std::endl;
        }
        global_k = global_k / vars::population_size;
        global_d = global_d / vars::population_size;

        double k_within_var{};
        double d_within_var{};

        double k_between_var{};
        double d_between_var{};

        double k_total_var{};
        double d_total_var{};

        for (int i = 0; i < vars::population_size; ++i)
        {
            std::vector<cell> test = m_population.at(i).get_cell_list();

            double k_mean = m_population.at(i).get_fitness();
            double d_mean = m_population.at(i).get_d_average();

            k_between_var += pow(k_mean - global_k, 2);
            d_between_var += pow(d_mean - global_d, 2);

            for (cell &c : test)
            {
                k_within_var += pow(c.get_pg() - k_mean, 2);
                d_within_var += pow(c.get_dummy() - d_mean, 2);

                k_total_var += pow(c.get_pg() - global_k, 2);
                d_total_var += pow(c.get_dummy() - global_d, 2);
            }
        }
        k_within_var = k_within_var / (vars::cells_per_org * vars::population_size);
        d_within_var = d_within_var / (vars::cells_per_org * vars::population_size);

        k_between_var = k_between_var / vars::population_size;
        d_between_var = d_between_var / vars::population_size;

        k_total_var = k_total_var / (vars::cells_per_org * vars::population_size);
        d_total_var = d_total_var / (vars::cells_per_org * vars::population_size);


        if (vars::output == true)
        {  
            std::string file_name = std::to_string(vars::cells_per_org) + ".dat";
            std::ofstream outfile;
            outfile.open(file_name, std::ios::app);
            outfile << sim_steps << '\t' << global_k << '\t' << k_within_var << '\t' << k_between_var << '\t'
            << k_total_var << d_within_var << '\t' << d_between_var << '\t' << d_total_var << '\t' <<
            k_between_var / k_total_var << '\t' << d_between_var / d_total_var << '\t' << std::endl;
            outfile.close();
        }
    }

};

///simulation will have moran process at cell level, and wright-fisher at organism level. All cell traits are rest upon wright-fisher process (essentially new organisms)


int main()
{
    population new_pop;
    new_pop.init_population();
    long n_steps = vars::sim_steps;

    while (n_steps)
    {
        new_pop.growth_phase();
        --n_steps;

        if (n_steps % 1 == 0)
        {
            new_pop.output_vals();
            new_pop.output_graph();
            new_pop.to_file(vars::sim_steps - n_steps);
            std::cout << vars::sim_steps - n_steps << " wright-fisher processes complete" << "\n" << std::endl;

        }
        new_pop.fisher_process();
        new_pop.new_gen_reset();
        
        
    }   

    return 0;
}
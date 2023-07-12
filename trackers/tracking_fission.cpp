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


    /// 1 / mut_rate is the chance of mutation upon reproduction
    constexpr double mut_rate{100};


    ///the number of cells in an organism required for fission
    constexpr int threshold{1000};

    //starting parameters
    constexpr int start_size{mut::threshold / 2};
    constexpr int start_pop{150};

    ///total cells in the system
    long total_cells{threshold * 100};

    ///run steps
    constexpr long time_steps{100000000};

    ///returns on cooperation
    constexpr float alpha{1.0};

    ///returns on personal reproduction
    constexpr float beta{1.0};

    constexpr bool output_returns{false};
    
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

    /// 0 = all traits mutate, 1 = d mutates only, 2 = no traits mutate
    constexpr int selection_type{0};

    constexpr bool start_broken{true};
    constexpr double a_pg{0};
    constexpr double a_switch{0.46};
    constexpr double b_pg{0.56};
    constexpr double b_switch{0.000005};   


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

    int out_freq{};

    // for ancestor tracking
    constexpr int ancestor_length{10000};
    constexpr int tracking_freq{50};
    vector<int> gen_times{};


    constexpr bool track_codes{true};
    constexpr int n_codes{1000};
    constexpr int code_freq{50};


}

///RNG
auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
mt19937 mersenne( static_cast<mt19937::result_type>(seed) );

///Global distributions
uniform_real_distribution<double> double_num(0.0, 1.0);
uniform_int_distribution<> mut_dist(1, mut::mut_rate);
uniform_int_distribution<> splitter(0, 1);
normal_distribution<double> curve(0, 0.02);

///integer RNG
long generate_random_int(const int min, const int max)
{
    uniform_int_distribution<> num(min, max);
    long prob = num(mersenne);
    return prob;
}

///mutate the chosen genotype, accounting for boundaries
void normal_mutation(double &gene)
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

///calculate poisson distributed death rate based on the previous generation cell average, adapted to the size of the organism
int death_calculator(unsigned int organism_size)
{
    double death_calc = mut::death_rate * (static_cast<double>(organism_size) / static_cast<double>(mut::current_cells));  ///calculate death rate based on the previous generation cels, adapted to the size of the organism
    poisson_distribution<int> p_distribution(death_calc); ///poisson distribution for organism deaths
    int deaths = p_distribution(mersenne);
    return deaths;
}

///probably don't use this, Need variation through poisson distribution. 
int death_calculator_deterministic(unsigned int organism_size)
{
    int deaths = ceil( mut::death_rate * (static_cast<double>(organism_size) / static_cast<double>(mut::current_cells)));  ///calculate death rate based on the previous generation cels, adapted to the size of the organism
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

    vector<char> ancestors;



public:
    void set_cell(char type, double replication1, double catalysis1, double conversion1, double replication2, double catalysis2, double conversion2, double dummy, bool tag=false)
    {
        m_type = type;
        ma_replication = replication1;
        ma_conversion = conversion1;
        ma_catalysis = catalysis1;
        mb_replication = replication2;
        mb_conversion = conversion2;
        mb_catalysis = catalysis2;
        m_dummy = dummy;
        m_tag = tag;
    }

    double get_tag()
    {
        return m_tag;
    }

    double get_dummy()
    {
        return m_dummy;
    }

    void set_ancestors_vector()
    {
        ancestors.assign(mut::ancestor_length, 'B');
    }

    void set_ancestor()
    {
        ancestors.insert(ancestors.begin(), m_type);
        ancestors.pop_back();
    }

    vector<char>& get_ancestors()
    {
        return ancestors;
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
            switch (mut::selection_type)
            {
                case 0:
                {
                    normal_mutation(m_dummy);

                    mutation(ma_replication);
                    ma_catalysis = 1.0 - ma_replication;

                    mutation(ma_conversion);

                    mutation(mb_replication);
                    mb_catalysis = 1.0 - mb_replication;

                    mutation(mb_conversion);
                    break;
                }
                case 1:
                {
                    normal_mutation(m_dummy);
                    mutation(ma_conversion);
                    break;                
                }
                case 2:
                {
                    normal_mutation(m_dummy);
                    break;
                
                }
            }
        }
    }
};



///cells in organism are contained within 'm_current'
class organism
{
    vector<cell> m_current;
    /// check how many mutations are happening per life cycle
    long mut_num{};
    int gen_time{};

public:

    organism()
    {
        m_current.reserve(mut::threshold);
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

    int growth()
    {
        ++gen_time;
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
            m_current[number].set_ancestor();
            m_current.back().set_ancestor();
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
            int number = generate_random_int(0, m_current.size() - 1);
            auto it = m_current.begin() + number;
            *it = move(m_current.back());
            m_current.pop_back();
            
            // m_current.erase(m_current.begin() + number);
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
        mut::gen_times.push_back(gen_time);
        gen_time = 0;

        vector<cell> new_organism{};
        vector<cell> new_copy{};
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
    }

    cout << largest_org <<  "     population size: " << total_cell_count << "   number of orgs: " << m_population.size() <<
    "   approx number of cells before death: " << mut::current_cells + mut::new_death_rate << endl;

}








void ancestor_tracking(vector<organism>& m_population)
{
    long int p_ancestor{};
    long int q_ancestor{};
    long int p_count{};
    long int q_count{};

    for (unsigned int i = 0; i < m_population.size(); ++i)
    {
        std::vector<cell> test = m_population.at(i).get_organism();
        for (cell &c : test)
        {
            char current_state = c.get_type();
            if (current_state == 'A')
                ++ p_count;
            else
                ++ q_count;

            std::vector<char> ancestors = c.get_ancestors();
            char last = ancestors.back();
            if (last == 'A')
                ++ p_ancestor;
            else
                ++ q_ancestor;
        }
    }
    double mean = 0;
    if (mut::gen_times.size() > 0)
    {
        for (int i : mut::gen_times)
            mean += i;
        mean = mean / (double)(mut::gen_times.size());
        cout << mean << "!!" << endl;
    }

    std::ofstream outfile;
    string name = to_string(mut::threshold) + "tracking.dat";
    outfile.open(name, std::ios::app);
    outfile << p_count << "\t" << q_count << "\t" << p_ancestor << "\t" << q_ancestor << "\t" << mean << "\t" << mut::gen_times.size() << std::endl;
    outfile.close();
}







///contains the population of organisms within a vector
class population
{
    vector<organism> m_population;

public:

    population()
    {
        m_population.reserve( mut::total_cells / (mut::threshold / 4) );
    }

    void set_population(vector<organism>& total)
    {
        m_population = total;
        if (mut::invasion)
        {
            mut::out_freq = 10;
        }
        else
        {
            mut::out_freq = 1000;
        }
    }

    vector<organism> get_population()
    {
        return m_population;
    }

    void process()
    {
        long time_steps = mut::time_steps;
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
                        // m_population.erase(m_population.begin() + i);
                        // efficiently kill the organism. 
                        auto it = m_population.begin() + i;
                        *it = move(m_population.back());
                        m_population.pop_back();
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
                int num = generate_random_int(0, m_population.size() - 1);
                // efficiently kill the organism. 
                auto it = m_population.begin() + num;
                *it = move(m_population.back());
                m_population.pop_back();
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
            
            if (time_steps % mut::tracking_freq == 0)
            {
                ancestor_tracking(m_population);
            }


            if (time_steps % 50000 == 0)
            {
                mut::gen_times.clear();
            }

            --time_steps;
        }

    }
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


    int n_orgs = static_cast<int>(ceil((mut::total_cells * 0.85) / mut::start_size));

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




int main()
{
    population environment;
    if (mut::invasion)
    {
        environment = invasion();
    }
    else
    {
        cell one;
        cell two;

        one.set_ancestors_vector();
        two.set_ancestors_vector();

        if (mut::start_broken == true)
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
    environment.process();
        
    return 0;
}


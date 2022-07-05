#include <string>

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)


double RANDOM();
int Seed(int seeder);
long RandomNumber(long max);
int Randomize(void);

using namespace std;

namespace mut
{
    /// prob of mutation
    extern const double mut_rate;

    ///the number of cells in an organism required for fission
    extern const int threshold;

    //starting parameters
    extern const int start_size;
    extern const int start_pop;

    ///total cells in the system
    extern long total_cells;

    extern int out_freq;

    ///run steps
    extern const long time_steps;;

    ///returns on cooperation
    extern const float alpha;

    ///returns on personal reproduction
    extern const float beta;

        
    /// trait essentiality
    extern const double trait_e;

    ///the time step divisor to approximate a continuous model
    extern const double delta_factor;

    ///maximum organisms in population - moran process, turn on or off, used for case where there is no selection so there is collective dynamics.
    extern const bool moran_on;

    /// 0 = all traits mutate, 1 = d mutates only, 2 = all traits mutate
    extern const int selection_type;

    //initial conditions if needed. 
    extern const bool start_broken;
    extern const double a_pg;
    extern const double a_switch;;
    extern const double b_pg;
    extern const double b_switch;
    
    ///what data files to send out
    extern const bool output_rel;
    extern const bool output_var;
    extern const bool output_values;
    extern const bool output_type_diff;
    
    // change file names to alpha and beta (returns)
    extern const bool output_returns;
    extern string al;
    extern string be;
    extern string a;
    extern string b;
    extern string returns;


    ///death rate is calculated based on current cells from the previous time step and an average of the "r" values
    extern double death_rate;
    extern double new_death_rate;
    extern long current_cells;
    extern long new_current_cells;

    extern const int s;

    extern const  bool invasion;
    /// 0 = 50/50, 1 = suppressor invades, 2 = fittest invades, 3 = choose frequency
    extern const int inv_type;
    extern const int inv_freq;

    // for budding simulation
    extern long propagules;    
}


class cell
{
/// cell type can either be 'A' or 'B', each has 3 values (replication = 'r', catalysis = '1 - r', conversion = 'c'
private:
    char m_type;
    double ma_replication;
    double ma_catalysis;
    double ma_conversion;
    double mb_replication;
    double mb_catalysis;
    double mb_conversion;
    double m_dummy;
    bool m_tag;
    std::vector<double> codes;
public:

    inline void set_cell(char type, double replication1, double catalysis1, double conversion1, double replication2, double catalysis2, double conversion2, double dummy, bool tag=false)
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

    inline void set_code_vector()
    {
        codes.assign(80, 50.0);
    }


    void set_code(); ///set a random <double> code at the front of the cell, and remove the one at the back. Used to determine close relatives.

    inline std::vector<double>& get_codes()
    {
        return codes;
    }

    inline double get_tag()
    {
        return m_tag;
    }

    inline double get_dummy()
    {
        return m_dummy;
    }

    inline void set_type(char type)
    {
        m_type = type;
    }

    inline double get_aconversion()
    {
        return ma_conversion;
    }

    inline double get_areplication()
    {
        return ma_replication;
    }

    inline double get_acatalysis()
    {
        return ma_catalysis;
    }


    inline double get_bconversion()
    {
        return mb_conversion;
    }

    inline double get_breplication()
    {
        return mb_replication;
    }

    inline double get_bcatalysis()
    {
        return mb_catalysis;
    }

    inline char get_type()
    {
        return m_type;
    }

    void mutate();
};

class organism
{
private:
    vector<cell> m_current;
    /// check how many mutations are happening per life cycle
    long mut_num{};
    // propagules used for budding
    vector<organism> propagules;
public:
    void set_organism(vector<cell>& total);
    vector<cell>& get_organism();
    int get_organism_population();
    int growth();
    
    ///split the organism into two. return the new organism
    vector<cell> split();


    //budding methods (renditions on growth + split)
    vector<cell> budding_split();
    int budding_growth();

};

///contains the population of organisms within a vector
class population
{
private:
    std::vector<organism> m_population;

public:
    void set_population(std::vector<organism>& total);
    std::vector<organism> get_population();
    void process(long time_steps);

    // budding rendition on method "process"
    void budding_process(long time_steps);
};

population invasion(void);
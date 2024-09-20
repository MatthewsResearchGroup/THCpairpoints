#ifndef _ENERGYANALYSIS_JOBINFO_HPP_
#define _ENERGYANALYSIS_JOBINFO_HPP_

#include <string>
#include <ctime>
#include <filesystem>

/********************************************
 * Jobinfo
 *
 * struct that contains information for an 
 * energyAnalysis job. 
 *
********************************************/
struct Jobinfo
{
    int num_points;
    int num_pairs;
    int seed = (int) time(nullptr);
    std::string molecule;
    std::string method;
    std::string basis;
    std::string curr_working_path;
    std::string path_to_orb;
    std::string path_to_grid;
    std::string path_to_fa;
    std::string path_to_fi;
    std::string path_to_xa;
    std::string path_to_xi;
    std::string path_to_v;
    std::string path_to_t;
    std::string path_to_csv;

    std::string print_str() const;
    int validate() const;

    void set_molecule(const std::string& mol);
    void set_method(const std::string& cal_method);
    void set_basis(const std::string& bas);
    void set_work_path();
    void set_orb_path();
    void set_grid_path();
    void set_t_path();
    void set_csv_path();
    //void set_points_set_path(const std::string& path);
};


/********************************************
 *
 * print_str()
 *
 * Generates a std::string to be used in printing 
 * information about the job
 *
********************************************/
std::string Jobinfo::print_str() const
{
    std::string str;

    str += "Fa data path      : " + path_to_fa + "\n";
    str += "Fi data path      : " + path_to_fi + "\n";
    str += "V  data path      : " + path_to_v + "\n";
    str += "t  data path      : " + path_to_t + "\n";
    str += "Xa data path      : " + path_to_xa + "\n";
    str += "Xi data path      : " + path_to_xi + "\n";
    str += "CSV data path     : " + path_to_csv + "\n";

    return str;
}

/*******************************************
 * validate()
 *
 * validates the jobinfo. Should be used after 
 * input parsing. Returns 0 only if there is no error
 * detected 
*******************************************/
int Jobinfo::validate() const 
{
    if (num_points < 0) {return 1;}
    if (num_pairs <= 0) {return 1;}

    return 0;
}

/*******************************************
 * get the current working directory
 *
 *
******************************************/

void Jobinfo::set_work_path()
{
    if (basis != "dz")
        curr_working_path = std::filesystem::current_path().c_str() + std::string("/") + basis;
    else
        curr_working_path = std::filesystem::current_path();
}


/*******************************************
 * get the calculation method
 *
******************************************/
void Jobinfo::set_method(const std::string& cal_method)
{
    method = cal_method;
}

void Jobinfo::set_molecule(const std::string& mol)
{
    molecule = mol;
}

void Jobinfo::set_basis(const std::string& bas)
{
    basis = bas;
}

/*******************************************
 * set_orb_path(std::string&)
 *
 * Sets data based on the orbital directory path
*******************************************/
void Jobinfo::set_orb_path() 
{
    if (method == "MP2" or method == "MP3")
        path_to_orb = curr_working_path;
    else if (method == "CCSD")
        path_to_orb = curr_working_path + std::string{"/ccsd"};

    path_to_fa = path_to_orb + std::string{"/fa.dat"};
    path_to_fi = path_to_orb + std::string{"/fi.dat"};
    path_to_v  = path_to_orb + std::string{"/v.dat"};
}

/*******************************************
 * set_grid_path(std::string&)
 *
 * Sets data based on the grid directory path
*******************************************/
void Jobinfo::set_grid_path() 
{
    path_to_grid = curr_working_path + std::string{"/grid"};
    path_to_xa = path_to_grid + std::string{"/XA.dat"};
    path_to_xi = path_to_grid + std::string{"/XI.dat"};
}


/*******************************************
 * set_csv_path(std::string&)
 *
 * Sets csv file based on the project root directory path
*******************************************/
void Jobinfo::set_csv_path()
{
    
    path_to_csv = curr_working_path + "/partialenergyanalysis" + "_" + method + "_" + basis  + ".csv";
    // if (basis == "dz")
    //     path_to_csv = curr_working_path + "/../../partialenergyanalysis.csv";
    // else
    //     path_to_csv = curr_working_path + "/../../../partialenergyanalysis.csv";
}

/*******************************************
 * set_amp_path(std::string&)
 *
 * Sets data based on the project root directory path
*******************************************/
void Jobinfo::set_t_path()
{
    if (method == "MP3")
        path_to_t = curr_working_path + std::string{"/t.dat"};
    else if (method == "CCSD")
        path_to_t = curr_working_path + std::string{"/ccsd/t.dat"};
}
#endif

/************************************************* 
 * input.hpp 
 *
 * Input parser for partialGridenergyAnalysis.cxx
 *
 *  Uses the docopt library found at:
 *  https://github.com/docopt/docopt.cpp
*************************************************/

#ifndef _ENERGYANALYSIS_INPUT_HPP_
#define _ENERGYANALYSIS_INPUT_HPP_

#include "jobinfo.hpp"
#include "docopt.h"

#include <string>
#include <map>
#include <iostream>

/************************************************* 
 *  String for docopt
*************************************************/ 
static const char USAGE[] =
R"(EnergyAnalysis.

  Usage:
    partialGridenergyAnalysis --molecule=<molecule>  --method=<method>  --basis=<basis>  
    partialGridenergyAnalysis (-h | --help)
    partialGridenergyAnalysis --version

  Options:
    -h --help    Show this screen.
     --version   Show version.
)";


/************************************************* 
 *  Parse and check input 
*************************************************/ 
int parse_input(const std::map<std::string,docopt::value> args, 
                Jobinfo& jobinfo)
{
    for (auto const& arg : args) 
    {

//Turn this on if you want debugging
#if 1
       std::cout << arg.first << ": " << arg.second << std::endl; 
#endif

       // moleucle options
       if (arg.first == "--molecule") {jobinfo.set_molecule(arg.second.asString());}

       // Method options
       if (arg.first == "--method") {jobinfo.set_method(arg.second.asString());}

       // basis options
       if (arg.first == "--basis") {jobinfo.set_basis(arg.second.asString());}


       // Path options, set relative path.
       jobinfo.set_work_path();
       jobinfo.set_orb_path();
       jobinfo.set_grid_path();
       jobinfo.set_t_path();
       jobinfo.set_csv_path();

    }//end loop over arguements

    printf("Job Information\n%s\n",jobinfo.print_str().c_str());

    return jobinfo.validate();
}


#endif

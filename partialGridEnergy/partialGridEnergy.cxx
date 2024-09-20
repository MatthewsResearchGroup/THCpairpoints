#include "marray.hpp"
#include "docopt.h"
#include "jobinfo.hpp"
#include "input.hpp"
#include "io.hpp"
#include "qc_utility.hpp"
#include "tensor_ops.hpp"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <stdio.h>
#include <iomanip>
#include <vector>
#include <chrono>
#include <random>

using namespace std;

int main(int argc, char **argv)
{
    std::map<std::string, docopt::value> args = docopt::docopt(USAGE, {argv+1, argv + argc}, true, "patialGridEnergy 1.0");

    // Parse the input and setup the jobinfo struct
    //
    Jobinfo jobinfo;
    if (parse_input(args, jobinfo) != 0)
    {
        printf("Bad input to partialGridenergyanalysis\n");
        exit(1);
    } 

    //read the dimensions we need
    // no 	     -> number of occupied orbitals
    // nv	     -> number of virtual orbitals
    // nps	     -> number of TOTAL grid points
    int nv, no, nps;
    read_dimensions(nv,  jobinfo.path_to_fa,
                    no,  jobinfo.path_to_fi,
                    nps, jobinfo.path_to_xa);
    
    // int nvo = nv * no;
    printf("nvrt: %d \nnooc : %d \nngrd : %d \n", nv, no, nps);

    // the rotio of # of grid points over # of orbitals, which define the number of grid points to keep
    double ratio = nps / (no + nv) ;
    int norb = no + nv ;

    // if (nps < 10 * ( nv + no))
    // {
    //     printf("# of grid points is less than 10 times of # of (no + nv)\n");
    //     exit(1);
    // }
    

    // srand(time(NULL));
    // int direc = 8;
    // char buffer [50];

    // cout << setprecision(15);

    // read xa, xi, fa, fi from files
    auto xpa = tensor_from_file(jobinfo.path_to_xa, sizeof(Int), nps, nv);
    auto xpi = tensor_from_file(jobinfo.path_to_xi, sizeof(Int), nps, no);

    auto faa = tensor_from_file(jobinfo.path_to_fa, sizeof(Int), nv, nv); 
    auto fii = tensor_from_file(jobinfo.path_to_fi, sizeof(Int), no, no); 
    auto FA = to_diagonal(faa);
    auto FI = to_diagonal(fii);
    
    //Read in vc,vx
    auto Vaibj = aibj_from_file(jobinfo.path_to_v, 0, nv, no); 
    auto Vajbi = aibj_to_ajbi(Vaibj);
    auto VCmm = Vaibj.lowered(2);
    auto VXmm = Vajbi.lowered(2);

    // We need to calculate the t from v in MP2 method, however just read it from t.dat file for MP3 and CCSD method.
    tensor<4> Taibj{nv, no, nv, no};
    if (jobinfo.method == "MP2")
    {
        // printf("Calculation method is MP2, making T from V...\n");
        make_c(Vaibj, FA, FI, Taibj);
    }
    else if (jobinfo.method == "MP3")
    {
        // first order of amplitude from MP2, second order of amplitude from the file
        make_c(Vaibj, FA, FI, Taibj);
        auto Taibj_2 = aibj_from_file(jobinfo.path_to_t, 0, nv, no);
        Taibj += Taibj_2;
    }
    else if (jobinfo.method == "CCSD")
    {
        // printf("Calculation method is MP3 or CCSD. Reading t from file: %s\n", jobinfo.path_to_t.c_str());
        Taibj = aibj_from_file(jobinfo.path_to_t, 0, nv, no);
    }
    auto Tmm = Taibj.lowered(2);

    /*************************************
     *
     * Energy clalculations start!
     *
     * **********************************/

    double E_exact_c, E_exact_x;
    E_exact_c = 2 * E1(Vaibj, Taibj);
    E_exact_x = -E1(Vajbi, Taibj);
    // if (jobinfo.method == "MP2")
    // {
    //     E_exact_c = 2 * E1(Vaibj, Vaibj, FA, FI);
    //     E_exact_x = -E1(Vaibj, Vajbi, FA, FI);
    // }
    // else if (jobinfo.method == "MP3")
    // {
    //     E_exact_c = 2 * E1(Vaibj, Vaibj, FA, FI, Taibj)
    //     E_exact_c = -E1(Vaibj, Vajbi, FA, FI, Taibj)
    // }
    // else if (jobinfo.method == "CCSD")
    // {    
    //     E_exact_c = 2 * E1(Vaibj, Taibj);
    //     E_exact_x = -E1(Vajbi, Taibj);
    // }

    // printf("E_exact_c : %f\nE_exact_x : %f\n", E_exact_c, E_exact_x);
    //******
    tensor<2> VCMM = Vaibj.lowered(2);
    tensor<2> VXMM = Vajbi.lowered(2);
    
    for (auto i = 1; i < 10 * ratio; i++)
    // for (auto num_grid_keep = 50; num_grid_keep < nps; num_grid_keep += 50)
    {
        // Num of grid points need to be kept, here we start from 0.1(nv+no) to 10(nv+no)
        auto num_grid_keep = int(i * (nv + no) / 10);
        double ratio_grid_orb = num_grid_keep / (double)(nv + no); // the ratio of # of grid over # of orbitals
        printf("num_of_keep : %d\n", num_grid_keep);

        /***************************************
         * Form Y   , krp of occ and virt gridpoints  
         *       MP' 
         *
         *  M  = a x i (outer product), MO basis
         *  P  = THC gridpoints
         *  P' = (potentially) reduced set of THC gridpoints
         *            
         *        
         *  Y    =  X    *  X
         *   p'ai    p'a     p'i
         *
         *          
         *  Y    =   Y     -> lower(ai) -> transpose
         *   MP'      P'ai
         *
         *  L2 normalize so that sum   Y      = 1
         *                          ai  aiP'
         ***************************************/

        tensor<2> Ymp = krp(xpa[range(num_grid_keep)][all],
                            xpi[range(num_grid_keep)][all]).lowered(1).T();   

    
        tensor<2> Ymp_nonorm = Ymp;
        /***************************************** 
         * Form Spp
         *          T
         * S   =   Y   .  Y
         *  PP      PM     MP
         *
        *****************************************/
        tensor<2> Spp = gemm(Ymp.T(), Ymp);


        /***************************************** 
         * Cholesky Decompose S 
         *                T
         * S   = SL   . SL
         *  PP     PP     PP
        *****************************************/
        tensor<2> SLpp = Spp;
        if (potrf('L', num_grid_keep, SLpp.data(), num_grid_keep) != 0) {printf("Bad exit on potrf\n"); exit(1);}

        /***************************************** 
         * Form B from the solution of linear eqs:
         *         T
         * B   . SL      =  Y
         *  MP     PP        MP
        *****************************************/
        tensor<2> Bmp = Ymp_nonorm;
        solve_tri(RIGHT, UPPER, SLpp.T(), Bmp); 

                 
        
        /***************************************** 
        * Form MC and MX
        *           (           )        (   T        )
        * M    =   (  T   . B    )  .   (   B   . V    )
        *  MM       (  MM    MP )        (   PM     MM)
        *****************************************/ 
        auto MCmm = gemm(gemm(Tmm, Bmp), gemm(Bmp.T(), VCmm));
        auto MXmm = gemm(gemm(Tmm, Bmp), gemm(Bmp.T(), VXmm));

        /***************************************** 
        *  Create E(2) contributions for each 
        *  pair (I think these are the non-orthogonalized
        *  pairs?)
        *                 T
        *  E2   =  diag( B    .  M     . B   ) 
        *    P            PM      MM      MP
        *****************************************/ 
        //How it is originally done in justenergy.cc.sb 
        // tensor<1> E2Cp =  2.e0 * gemm3diag(Bmp.T(), MCmm, Bmp); 
        // tensor<1> E2Xp = -1.e0 * gemm3diag(Bmp.T(), MXmm, Bmp); 
        auto E2Cp =  2.e0 * sum(gemm3diag(Bmp.T(), MCmm, Bmp)); 
        auto E2Xp = -1.e0 * sum(gemm3diag(Bmp.T(), MXmm, Bmp)); 


        // Percentage Error of Ex and Ec
        // auto PErr_c = abs( (E2Cp -  E_exact_c) /  E_exact_c);
        // auto PErr_x = abs( (E2Xp -  E_exact_x) /  E_exact_x);
        // auto PErrxminsc = PErr_x - PErr_c;

        output_to_csv(jobinfo,
                      num_grid_keep,
                      nps,       norb, 
                      E_exact_c, E_exact_x,
                      E2Cp,      E2Xp);
                      // PErr_c,    PErr_x,
                      // PErrxminsc);
    }
    return 0;
}

#include "energyAnalysis.hpp"
#include "marray.hpp"
#include "input.hpp"
#include "jobinfo.hpp"
#include "io.hpp"
#include "tensor_ops.hpp"
#include "pair_points.hpp"
#include "qc_utility.hpp"

#include "docopt.h"

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
#include <sstream>
#include <tuple>
#include <map>
#include <utility>

/********************************************
 * main function
 ********************************************/

int main(int argc, char **argv) {

    std::map<std::string, docopt::value> args = docopt::docopt(USAGE,
                                                 { argv+1, argv+argc },
                                                 true,
                                                 "Energy Analsyis 1.0");

    //Parse the input and setup the jobinfo struct
    Jobinfo jobinfo;
    if (parse_input(args,jobinfo) != 0) 
    {
        printf("Bad input to energyAnalysis\n");
        exit(1);
    }

    //read the dimensions we need
    // no 	     -> number of occupied orbitals
    // nv	     -> number of virtual orbitals
    // nps	     -> number of TOTAL grid points
    // num_pairs     -> number of pairs of subset of gridpoints to keep
    int nv, no, nps;
    std::vector<int> psvec;
    read_dimensions(nv,  jobinfo.path_to_fa,
                    no,  jobinfo.path_to_fi,
                    nps, jobinfo.path_to_xa);
                    // psvec, jobinfo.path_to_points_set);
    int nvo = nv * no;
    printf("nvrt : %d \nnocc : %d \nngrd : %d \n", nv, no, nps);
    printf("points_set : \n");
    for (auto i : range(31)) // we enlarge the scope of num_grid_keep to 6*(nv+no)
    {
        psvec.push_back(int(i * (nv+no) / 5));
    }
    // for (auto i : psvec)
    // {
    //     printf("%d\n",i);
    // }
    printf("\n");
    // num of pairs
    auto num_pairs     = jobinfo.num_pairs;
    double threshold = 1.0e-10;

    //read Fock matrix elements
    //Note that the fa, fi, Xa, and Xi files all have offsets of Int*1,
    //as the first entry in the file is the relevant dimension
    auto faa = tensor_from_file(jobinfo.path_to_fa, sizeof(Int), nv, nv); 
    auto fii = tensor_from_file(jobinfo.path_to_fi, sizeof(Int), no, no); 
    auto FA = to_diagonal(faa);
    auto FI = to_diagonal(fii);

    //Read in vc,vx
    auto Vaibj = aibj_from_file(jobinfo.path_to_v, 0, nv, no); 
    auto Vajbi = aibj_to_ajbi(Vaibj);
    auto VCmm = Vaibj.lowered(2);
    auto VXmm = Vajbi.lowered(2);
    // printf("checking the symmetry of VCmm...\n");
    sym_check(VCmm, threshold);
    // printf("checking the symmetry of VXmm...\n");
    sym_check(VXmm, threshold);



    //Generate T
    // if method is MP2, we will generate amplitude from V.
    // Otherwise, we will read amplitude from files.
    //
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

    //read the THC matrix elements (note offset from start of file, 
    //which skips the number of gridpoints 
    auto xpa = tensor_from_file(jobinfo.path_to_xa, sizeof(Int), nps, nv);
    auto xpi = tensor_from_file(jobinfo.path_to_xi, sizeof(Int), nps, no);

/*
    printf("FAA tensor is...\n"); jht_print(faa); 
    printf("FII tensor is...\n"); jht_print(fii); 
    printf("FA vector is...\n"); jht_print(FA);
    printf("FI vector is...\n"); jht_print(FI);
    printf("XA tensor is...\n"); jht_print(xpa); 
    printf("XI tensor is...\n"); jht_print(xpi); 
    printf("Vaibj matirx in C++ is...\n");jht_print(v.lowered(2));
    printf("T2 matrix is...\n"); jht_print(Tmm);
*/

    /*************************************
     * Form orthonormal basis of T
     *
     *  M = a x i (outer product), MO basis
     *  D = a x i (outer proudct), diagonal T basis
     *                       T
     *  T    =  C     S     C
     *   MM      MD    DD    DM
     *
     *  W   = |S|  * C    (column wise weighting)
     *   MD      D    MD
     *
     *  And we return the transpose of W
     *************************************/

    double Snorm; 
    auto res = weighted_evd(Tmm,Snorm);
    tensor<2> WT = std::get<0>(res);
    tensor<1> S = std::get<1>(res);
//    printf("Eigenvalues of T2 are...\n"); jht_print(S);
//    printf("Reweighted eigenvectors of T2^T are ...\n"); jht_print(WT);

    //Number of pair points
    printf("psvec size is %ld.\n", psvec.size());


    for (auto num_grid_keep : psvec)
    {

        // printf("the number of grid point : %d\n", num_grid_keep);

        if (num_grid_keep != 0)
        {
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
            //Constuct Ymp matrix in the reduced space:
            //Ymp = X[a][P]*X[i][p]/sqrt(sum_ai X[a][P]*X[i][P])
            tensor<2> Ymp = krp(xpa[range(num_grid_keep)][all],
                                xpi[range(num_grid_keep)][all]).lowered(1).T();   

            tensor<2> Ymp_nonorm = Ymp;

            L2_normalize(Ymp);
            // printf("JHT here 1\n");
               
//            printf("Ymp is...\n");jht_print(Ymp);

            /***************************************
             *  Project Y out of WT
             *
             *  M  = a x i (outer product), MO basis
             *  P  = THC gridpoints
             *  P' = (potentially) reduced set of THC gridpoints
             *  P" = orthogonal basis vectors
             *
             *  QR decompose Y:
             *
             *  Y    = Q    R
             *   MP'    MP"  P"P"
             *
             *  Project out
             *                              T
             *  WTP   = WTP   -  WTP   Q   Q
             *     DM      DM       DM  MP" P"M   
             *
             ****************************************/
            auto Q = qrdecomp(Ymp);
            tensor<2> QQT = gemm(Q, Q.T());
            // printf("Chao here we measure the demision of QQT: %ld, %ld\n", QQT.length(0), QQT.length(1));
            //for (auto i : range(QQT.length(0)))
                // printf("%ld : %e\n", i, QQT[i][i]);
            // tensor<2> WTP = WT - gemm3(WT,Q,Q.T()); 
            // printf("JHT here 2\n");

//            printf("Q is ...\n"); jht_print(Q);
//            printf("Q.QT is...\n"); jht_print(gemm(Q,Q.T()));
//            printf("WTP is... \n"); jht_print(WTP);
                
            std::vector<std::string> ftype_string = {"all", "cross", "nocross"};
            // std::vector<std::string> ftype_string = {"all"};
            for (auto ftype : ftype_string)
            {    
                printf("working on ....\n");   
                printf("method:                         %s\n", jobinfo.method.c_str());
                printf("ftype:                          %s\n", ftype.c_str());
                printf("num_grid_keep:                  %d\n", num_grid_keep);
                printf("number of pairs of grid points: %d\n", num_pairs);
                if (ftype ==  "all") {jobinfo.ftype = Jobinfo::Ftypes::all;}
                else if (ftype == "cross") {jobinfo.ftype = Jobinfo::Ftypes::cross;}
                else if (ftype == "nocross") {jobinfo.ftype = Jobinfo::Ftypes::nocross;}

                /* 
                * Generate the pair points
                * Here we use the full grid points to generate the pairs of grid points
                * which is not affected by the starting grid points, and benefits for 
                * the following analysis. 
                */
                printf("Generating %d pairs of grid points...\n", num_pairs);
                auto pairs = make_random_pairs(num_pairs, nps, jobinfo.seed);

                // printf("Pair list is...\n");
                // for (auto pair : pairs) 
                // {
                //     printf("[%d,%d] ", pair.first, pair.second);
                // }
                // printf("\n");

                printf("\n----------------------------------\n");



                /*****************************************   
                 * Form the pair points we want to test 
                 *
                 *  M  = a x i (outer product), MO basis
                 *  O  = pairs of THC gridpoints 
                 *  P  = THC gridpoints
                 *
                 * XXA       = X      +   X
                 *    O a       P a        Q a
                 *
                 * XXI       = X      +   X
                 *    O i       P i        Q i 
                 *
                 *
                 * Y2    = krp(XXA    , XXI   ) -> lower(1) -> transpose -> L2 norm
                 *   M O          O a      O i
                 *
                 *****************************************/   
                auto Y2mo = make_Y2(xpa,
                                    xpi,
                                    pairs,
                                    jobinfo);

                tensor<2> Y2mo_nonorm = Y2mo;
                tensor<2> Y2P = Y2mo - gemm3(Q, Q.T(), Y2mo);
                // for (auto i : range(Y2P.length(0)))
                // {
                //     printf("Y2P[%ld][0] = %e\n", i, Y2P[i][0]);
                // }
                

                tensor<1> Svals{num_pairs};
                for (auto pair : range(num_pairs)) {
                    double sum = 0;
                    for (auto ai : range(nvo)) {
                        sum += Y2mo[ai][pair] * Y2mo[ai][pair];
                    }
                    Svals[pair] = sqrt(sum);
                }

                // L2_normalize(Y2mo);
                L2_normalize(Y2P);

//                printf("Y2 matrix is...\n"); jht_print(Y2mo);

                /*****************************************   
                 * Form WPT.Y2 temp
                 *
                 * TMP   = WTP   Y2
                 *    DO      DM   MO
                 *
                 *           |     |^2
                 * F  = sum  |TMP  |    /  SNORM 
                 *  O      D |   DO|          
                 *****************************************/   
                 auto tmp = gemm(WT, Y2P); 
                 // for (auto i : range(tmp.length(0)))
                 // {
                 //     printf("%e\n", tmp[i][0]);
                 // }
                 auto Fo = sum_col_squares(tmp);
                 Fo *= 1e0/Snorm;
                 // printf("Snorm: %e\n", Snorm);
                 // printf("Fo size: %ld\n", Fo.length(0));
                 // for (auto i : range(Fo.length(0)))
                     // printf("%ld, %e\n", i, Fo[i]);

// E            ND F-VALUE ANALYSIS
//
// B            EGIN ENERGY ANALYSIS

                /* justenergy.cc code from Sara Beth */ 
                // It seems that "D" is equiv. to Y2mo, here,
                //   and "b" is 
                //
                //D = dU - B B^T dU
                //mu = dU^T dU - dU^T B B^T dU = dU^T D
                /***************************************** 
                 * Form Spp
                 *          T
                 * S   =   Y   .  Y
                 *  PP      PM     MP
                 *
                *****************************************/
                tensor<2> Spp = gemm(Ymp_nonorm.T(), Ymp_nonorm);


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
                 * Form Dmo via gemm3
                 *
                 *                        T
                 * D    =  Y2    - B   . B   . Y2
                 *  MO       MO     MP    PM     MO
                *****************************************/ 
                tensor<2> Dmo = Y2mo_nonorm;
                gemm3(-1.0e0, Bmp, Bmp.T(), Y2mo_nonorm, 1.0e0, Dmo);
//                 for (auto i : range(Dmo.length(0)))
//                 for (auto j : range(Dmo.length(1)))
//                     printf("%ld, %ld, %e\n", i, j, Dmo[i][j]);
// 

                /***************************************** 
                 * Form MUo 
                 *                  T
                 * mu    =  diag( Y2   . D   )
                 *   o              OM    MO
                *****************************************/ 
                auto MUo = gemmdiag(Y2mo_nonorm.T(), Dmo);

                
                /***************************************** 
                * Form MC and MX
                *           (           )        (   T        )
                * M    =   (  T   . B    )  .   (   B   . V    )
                *  MM       (  MM    MP )        (   PM     MM)
                *****************************************/ 
                auto MCmm = gemm(gemm(Tmm, Bmp), gemm(Bmp.T(), VCmm));
                auto MXmm = gemm(gemm(Tmm, Bmp), gemm(Bmp.T(), VXmm));
                // auto VXT_trace = sum(gemmdiag(VXmm, Tmm));
                // auto TVX_trace = sum(gemmdiag(Tmm, VXmm));
                // auto VCT_trace = sum(gemmdiag(VCmm, Tmm));
                // auto TVC_trace = sum(gemmdiag(Tmm, VCmm));
                // printf("VXT_trace, TVX_trace = %e, %e\n", VXT_trace, TVX_trace);
                // printf("VCT_trace, TVC_trace = %e, %e\n", VCT_trace, TVC_trace);


                // printf("Chao is checking the symmetry of MCmm\n");
                // sym_check(MCmm, threshold);

                // printf("Chao is checking the symmetry of MXmm\n");
                // sym_check(MXmm, threshold);

                // printf("MCmm size = %ld,%ld\n\n", MCmm.length(0), MCmm.length(1));
                // auto MCmm_VBBT = gemm(gemm(VCmm, Bmp), gemm(Bmp.T(), Tmm));
                // auto MXmm_VBBT = gemm(gemm(VXmm, Bmp), gemm(Bmp.T(), Tmm));
                 // printf("Here we check MCmm and MCmm_VBBT\n");
                 // mat_check(MCmm, MCmm_VBBT, threshold);
                // for (auto i : range(MCmm.length(0)))
                // for (auto j : range(MCmm.length(1)))
                // {
                //     if (MCmm[i][j] - MCmm_VBBT[i][j] > 1.e-5)
                //         printf("i,j, MCmm, MCmm_VBBT = %ld, %ld, %e, %e \n ", i, j, MCmm[i][j], MCmm_VBBT[i][j]);
                //     if (MXmm[i][j] - MXmm_VBBT[i][j] > 1.e-5)
                //         printf("i,j, MXmm, MCmm_VBBT = %ld, %ld, %e, %e \n ", i, j, MXmm[i][j], MXmm_VBBT[i][j]);
                // }


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
                // auto E2Cp_VBBT =  2.e0 * sum(gemm3diag(Bmp.T(), MCmm_VBBT, Bmp)); 
                // auto E2Xp_VBBT = -1.e0 * sum(gemm3diag(Bmp.T(), MXmm_VBBT, Bmp)); 

                /***************************************** 
                *  Create E(4) contributions for each 
                *  pair in O
                *                T
                *  E4   = diag( D   .  M     . D    ) / MU
                *    O           OM     MM     MO        O
                *****************************************/ 
                tensor<1> E4Co =  2.e0 * 2.e0 * gemm3diag(Dmo.T(), MCmm, Dmo) / MUo;
                tensor<1> E4Xo = -1.e0 * 2.e0 * gemm3diag(Dmo.T(), MXmm, Dmo) / MUo;
                // tensor<1> E4Co_VBBT =  2.e0 * 2.e0 * gemm3diag(Dmo.T(), MCmm_VBBT, Dmo) / MUo;
                // tensor<1> E4Xo_VBBT = -1.e0 * 2.e0 * gemm3diag(Dmo.T(), MXmm_VBBT, Dmo) / MUo;

                // auto Dmo_init = Dmo;
                // auto Dmo_T = Dmo.T();
                // auto M4Co = sum(gemm3diag(Dmo_T, MCmm, Dmo_init));
                // auto M4Co_VBBT = sum(gemm3diag(Dmo_T, MCmm_VBBT, Dmo_init));
                // auto M4Xo = sum(gemm3diag(Dmo_T, MXmm, Dmo_init));
                // auto M4Xo_VBBT = sum(gemm3diag(Dmo_T, MXmm_VBBT, Dmo_init));

                // printf("Here we check MCmm and MCmm_VBBT\n");
                // mat_check(MXmm, MXmm_VBBT, threshold);
                // printf("\n");

                // printf("Here we check M4Co and M4Co_VBBT\n");
                // printf("M4Co and M4Co_VBBT are %e, %e\n\n", M4Co, M4Co_VBBT);
                // printf("Here we check M4Xo and M4Xo_VBBT\n");
                // printf("M4Xo and M4Xo_VBBT are %e, %e\n\n", M4Xo, M4Xo_VBBT);

                // printf("Here we check E4Co and E4Co_VBBT\n");
                // printf("M4Co and M4Co_VBBT are %e, %e\n\n", E4Co[0], E4Co_VBBT[0]);
                // printf("E4Co and E4Co_VBBT size = %ld, %ld\n",  E4Co.length(0),  E4Co_VBBT.length(0));
                // printf("E4Co and E4Co_VBBT = %e, %e\n",  E4Co[0],  E4Co_VBBT[0]);
                // printf("E4Xo and E4Xo_VBBT size = %ld, %ld\n",  E4Xo.length(0),  E4Xo_VBBT.length(0));
                // printf("E4Xo and E4Xo_VBBT = %e, %e\n",  E4Xo[0],  E4Xo_VBBT[0]);
                // tensor<1> E4Co = 2.e0 * 2.e0 * gemm
                

                /***************************************** 
                *  Create E(8) contributions for each 
                *  pair in O
                *                  T
                *  E8   = diag(   D   .  V    . D    
                *    O             OM     MM     MO        
                *
                *               *
                *
                *                  T 
                *                 D    . T    . D    )
                *                  OM     MM     MO
                *****************************************/ 
                tensor<1> E8Co =  2.e0 * gemm3diag(Dmo.T(), VCmm, Dmo) 
                                       * gemm3diag(Dmo.T(), Tmm,  Dmo)
                                       / (MUo * MUo);
                tensor<1> E8Xo = -1.e0 * gemm3diag(Dmo.T(), VXmm, Dmo) 
                                       * gemm3diag(Dmo.T(), Tmm,  Dmo)
                                       / (MUo * MUo);
                //tensor<1> E8Co =  2.e0 * gemmdiag(gemm3(Dmo.T(), VCmm, Dmo), gemm3(Dmo.T(), Tmm,  Dmo)) 
                //                       / (MUo * MUo);
                //tensor<1> E8Xo = -1.e0 * gemmdiag(gemm3(Dmo.T(), VXmm, Dmo), gemm3(Dmo.T(), Tmm,  Dmo)) 
                //                       / (MUo * MUo);

                printf("jht guess right before output\n");
                output_to_csv(jobinfo,
                               pairs,
                               num_grid_keep,
                               Fo,
                               Svals,
                               MUo,
                               E2Cp, E2Xp,
                               // E2Cp_VBBT, E2Xp_VBBT,
                               E4Co, E4Xo,
                               // E4Co_VBBT, E4Xo_VBBT,
                               E8Co, E8Xo);
            }
        }
        else
        {
            std::vector<std::string> ftype_string = {"all", "cross", "nocross"};
            // std::vector<std::string> ftype_string = {"all"};
            for (auto ftype : ftype_string)
            {       
                printf("working on ....\n");   
                printf("method:                         %s\n", jobinfo.method.c_str());
                printf("ftype:                          %s\n", ftype.c_str());
                printf("num_grid_keep:                  %d\n", num_grid_keep);
                printf("number of pairs of grid points: %d\n", num_pairs);
                if (ftype ==  "all") {jobinfo.ftype = Jobinfo::Ftypes::all;}
                else if (ftype == "cross") {jobinfo.ftype = Jobinfo::Ftypes::cross;}
                else if (ftype == "nocross") {jobinfo.ftype = Jobinfo::Ftypes::nocross;}

                
                /* 
                * Generate the pair points
                * Here we use the full grid points to generate the pairs of grid points
                * which is not affected by the starting grid points, and benefits for 
                * the following analysis. 
                */
                printf("Generating %d pairs of grid points...\n", num_pairs);
                auto pairs = make_random_pairs(num_pairs, nps, jobinfo.seed);

                // printf("Pair list is...\n");
                // for (auto pair : pairs) 
                // {
                //     printf("[%d,%d] ", pair.first, pair.second);
                // }
                // printf("\n");

                // printf("\n----------------------------------\n");


                /*****************************************   
                 * Form the pair points we want to test 
                 *
                 *  M  = a x i (outer product), MO basis
                 *  O  = pairs of THC gridpoints 
                 *  P  = THC gridpoints
                 *
                 * XXA       = X      +   X
                 *    O a       P a        Q a
                 *
                 * XXI       = X      +   X
                 *    O i       P i        Q i 
                 *
                 *
                 * Y2    = krp(XXA    , XXI   ) -> lower(1) -> transpose -> L2 norm
                 *   M O          O a      O i
                 *
                 *****************************************/   
                auto Y2mo = make_Y2(xpa,
                                    xpi,
                                    pairs,
                                    jobinfo);

                tensor<2> Y2mo_nonorm = Y2mo;
                tensor<2> Y2P = Y2mo; 
                // for (auto i : range(Y2P.length(0)))
                // {
                //     printf("Y2P[%ld][0] = %e\n", i, Y2P[i][0]);
                // }
                

                tensor<1> Svals{num_pairs};
                for (auto pair : range(num_pairs)) {
                    double sum = 0;
                    for (auto ai : range(nvo)) {
                        sum += Y2mo[ai][pair] * Y2mo[ai][pair];
                    }
                    Svals[pair] = sqrt(sum);
                }

                // L2_normalize(Y2mo);
                L2_normalize(Y2P);

//                printf("Y2 matrix is...\n"); jht_print(Y2mo);

                /*****************************************   
                 * Form WPT.Y2 temp
                 *
                 * TMP   = WTP   Y2
                 *    DO      DM   MO
                 *
                 *           |     |^2
                 * F  = sum  |TMP  |    /  SNORM 
                 *  O      D |   DO|          
                 *****************************************/   
                 auto tmp = gemm(WT, Y2P); 
                 // for (auto i : range(tmp.length(0)))
                 // {
                 //     printf("%e\n", tmp[i][0]);
                 // }
                 auto Fo = sum_col_squares(tmp);
                 Fo *= 1e0/Snorm;

                /***************************************** 
                 * Form Dmo via gemm3
                 *
                 *                        T
                 * D    =  Y2    - B   . B   . Y2
                 *  MO       MO     MP    PM     MO
                *****************************************/ 
                tensor<2> Dmo = Y2mo_nonorm;

                /***************************************** 
                 * Form MUo 
                 *                  T
                 * mu    =  diag( Y2   . D   )
                 *   o              OM    MO
                *****************************************/ 
                auto MUo = gemmdiag(Y2mo_nonorm.T(), Dmo);

                /***************************************** 
                *  Create E(8) contributions for each 
                *  pair in O
                *                  T
                *  E8   = diag(   D   .  V    . D    
                *    O             OM     MM     MO        
                *
                *               *
                *
                *                  T 
                *                 D    . T    . D    )
                *                  OM     MM     MO
                *****************************************/ 
                tensor<1> E8Co =  2.e0 * gemm3diag(Dmo.T(), VCmm, Dmo) 
                                       * gemm3diag(Dmo.T(), Tmm,  Dmo)
                                       / (MUo * MUo);
                tensor<1> E8Xo = -1.e0 * gemm3diag(Dmo.T(), VXmm, Dmo) 
                                       * gemm3diag(Dmo.T(), Tmm,  Dmo)
                                       / (MUo * MUo);
                //tensor<1> E8Co =  2.e0 * gemmdiag(gemm3(Dmo.T(), VCmm, Dmo), gemm3(Dmo.T(), Tmm,  Dmo)) 
                //                       / (MUo * MUo);
                //tensor<1> E8Xo = -1.e0 * gemmdiag(gemm3(Dmo.T(), VXmm, Dmo), gemm3(Dmo.T(), Tmm,  Dmo)) 
                //                       / (MUo * MUo);

                // printf("jht guess right before output\n");

                output_to_csv(jobinfo,
                               pairs,
                               num_grid_keep,
                               Fo,
                               Svals,
                               MUo,
                               E8Co, E8Xo);
        }
        }
    }
}

/***************************************************
 * io.hpp
 *
 * Header file for io operators in energyAnalysis.cxx
 *
 * Current list of functions:
 *     read_dimensions
 *     tensor_from_file
 *     aibj_from_file 
 *     output_to_csv
****************************************************/
#ifndef _ENERGYANALYSIS_IO_HPP_
#define _ENERGYANALYSIS_IO_HPP_

#include "partialGridEnergy.hpp"
#include "marray.hpp"

#include <string>
#include <stdio.h>
#include <fstream>
#include <filesystem>


/***************************************************
 * read_dimensions
 *
 * reads the dimensions of the neccesary tensors  
***************************************************/
//Read the relevant dimensions from the data files
void read_dimensions(Int& nv,  const std::string& nv_file_path,
                     Int& no,  const std::string& no_file_path,
                     Int& nps, const std::string& nps_file_path) {
                     // std::vector<int>& psvec, const std::string& ps_file_path) {

    //Read number of virtuals 
    auto nv_file = fopen(nv_file_path.c_str(), "rb");
    if (!nv_file) {
        printf("Could not open %s \n",nv_file_path.c_str());
        exit(1);
    }
    fread(&nv, sizeof(Int), 1, nv_file);
    if (fclose(nv_file) != 0) {
        printf("Error closing %s\n", nv_file_path.c_str());
        exit(1);
    }

    //Read number of occupied
    auto no_file = fopen(no_file_path.c_str(), "rb");
    if (!no_file) {
        printf("Could not open %s \n", no_file_path.c_str());
        exit(1);
    }
    fread(&no, sizeof(Int), 1, no_file);
    if (fclose(no_file) != 0) {
        printf("Error closing %s\n", no_file_path.c_str());
        exit(1);
    }

    //Read number of gridpoInts
    auto nps_file = fopen(nps_file_path.c_str(),"rb");
    if (!nps_file) {
        printf("Could not open %s \n", nps_file_path.c_str());
        exit(1);
    }
    fread(&nps, sizeof(Int), 1, nps_file);
    if (fclose(nps_file) != 0) {
        printf("Error closing %s\n", nps_file_path.c_str());
        exit(1);
    }

    // Read the points sets
    // std::string s;
    // std::fstream ps_file;
    // ps_file.open(ps_file_path.c_str());
    // if (!ps_file.is_open())
    // {
    //     printf("Could not open %s \n", ps_file_path.c_str());
    //     exit(1);
    // }
    // std::getline(ps_file, s); // Get rid of the first row, which may be the name of data
    // while(std::getline(ps_file, s))
    // {
    //     printf("%s\n",s.c_str());
    //     psvec.push_back(stoi(s));
    // }
}

/***************************************************
 * tensor_from_file
 * 
 * Read a (double) tensor from a file, with an offset 
***************************************************/
tensor<2> tensor_from_file(const std::string& file_name,
                           const size_t offset,
                           const Int n0,
                           const Int n1) {
    tensor<2> data{{n0, n1},COLUMN_MAJOR};

    auto data_file = fopen(file_name.c_str(),"rb");
    if (!data_file) {
        printf("Error opening %s \n",file_name.c_str());
    }

    if (fseek(data_file, offset, SEEK_SET) != 0) {
        printf("Error seeking to %zu in %s \n",offset,file_name.c_str());
    } 

    for (auto j : range(data.length(1))) {
        fread(&data[0][j], sizeof(double), data.length(0), data_file);
    }

    if (fclose(data_file) != 0) {
        printf("Error closing %s \n",file_name.c_str());
    }

    return data;
    
}

tensor<4> tensor_from_file(const std::string& file_name,
                           const size_t offset,
                           const Int n0,
                           const Int n1,
                           const Int n2,
                           const Int n3) {

    tensor<4> data{{n0, n1, n2, n3},COLUMN_MAJOR};

    auto data_file = fopen(file_name.c_str(),"rb");
    if (!data_file) {
        printf("Error opening %s \n",file_name.c_str());
    }

    if (fseek(data_file, offset, SEEK_SET) != 0) {
        printf("Error seeking to %zu in %s \n",offset,file_name.c_str());
    } 

    for (auto l : range(data.length(3)))
    for (auto k : range(data.length(2)))
    for (auto j : range(data.length(1)))
    {
        fread(&data[0][j][k][l], sizeof(double), data.length(0), data_file);
    }

    if (fclose(data_file) != 0) {
        printf("Error closing %s \n",file_name.c_str());
    }

    return data;
}

/***************************************************
 * aibj_cx_from_file
 *
 * Reads a tensor stored (vrt, occ, vrt, occ) from a 
 * path. 
****************************************************/
tensor<4> aibj_from_file(const std::string& file_name,
                         const size_t offset,
                         const Int nv,
                         const Int no) {
                         
    tensor<4> data{{nv, no, nv, no},COLUMN_MAJOR};
    
    auto data_file = fopen(file_name.c_str(),"rb");
    if (data_file == NULL) {
        printf("Error opening %s \n",file_name.c_str());
    }   
    
    if (fseek(data_file, offset, SEEK_SET) != 0) {
        printf("Error seeking to %zu in %s \n",offset, file_name.c_str());
    } 

    for (auto j : range(no))
    for (auto i : range(no))
    for (auto b : range(nv)) 
    {
        fread(&data[0][i][b][j], sizeof(double), nv, data_file);
    }

    if (fclose(data_file) != 0) 
    {
        printf("Error closing %s \n",file_name.c_str());
    }

    return data;
    
}

/***************************************************
 * output_to_csv
 *
 * writes the output of the calculation to a files called
 * partialGridEnergyAnalysis.csv 
***************************************************/
void output_to_csv(const Jobinfo& jobinfo,
              int& num_grid_keep,
              int& ngrid,        int& norb,
              double& Ec_exact,  double& Ex_exact,
              double& Ec_THC,    double& Ex_THC)
              // double& PErr_c,    double& PErr_x,
              // double& PErrxminsc)
 {
     std::string filename = jobinfo.path_to_csv;

     auto out_csv = fopen(filename.c_str(), "a+");

    if (!out_csv) 
    {
        printf("\nERROR: could not open %s for output\n",filename.c_str());
        exit(1);
    }

    // Generate formatting for output
    std::string header   = "Molecule,Method";
    std::string values_f = "%s, %s";

    header    += ",Basis,num_grid_keep";
    values_f  += ",%s,%6d";

    header    += ",ngrid,norb";
    values_f  += ",%6d,%6d";

    header    += ",Ec_exact,Ex_exact";
    // values_f  += ",%20.17g,%20.17g";
    values_f  += ",%e,%e";

    header    += ",Ec_THC,Ex_THC";
    values_f  += ",%e,%e";
    // values_f  += ",%20.17g,%20.17g";

    // header    += ",PErr_c,PErr_c";
    // values_f  += ",%e,%e";

    // header    += ",PErrxminsc";
    // values_f  += ",%e";


    values_f  += "\n";


    if (std::filesystem::is_empty(filename))
        fprintf(out_csv,"%s\n",header.c_str());

    fprintf(out_csv, values_f.c_str(),
            jobinfo.molecule.c_str(),  jobinfo.method.c_str(),
            jobinfo.basis.c_str(),     num_grid_keep,
            ngrid,                     norb,
            Ec_exact,                  Ex_exact,
            Ec_THC,                    Ex_THC);
            // PErr_c,                    PErr_x,
            // PErrxminsc);  

    fclose(out_csv);
 }
#endif

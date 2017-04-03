#include "main.h"

void mol_selection (int& type_mol_1, int& type_mol_2, int& apm_mol_1, int& apm_mol_2, vector<double>& M1, vector<double>& M2, double& totM1, double& totM2)
{
        if (type_mol_1 == 1) // Type "Ethylene Oxide"
        {
            apm_mol_1 = 3;
            M1.push_back(14.026);
            M1.push_back(14.026);
            M1.push_back(15.999);
            totM1 = M1[0] + M1[1] + M1[2];
        }
        else if (type_mol_1 == 2) // Type "Methanol"
        {
            apm_mol_1 = 3;
            M1.push_back(15.034);
            M1.push_back(15.999);
            M1.push_back(1.008);
            totM1 = M[0] + M[1] + M[2];
        }
        else if (type_mol_1 == 3) // Type "CO2"
        {
            apm_mol_1 = 3;
            M1.push_back(12.01);
            M1.push_back(15.999);
            M1.push_back(15.999);
            totM1 = M1[0] + M1[1] + M1[2];
        }
        else
        {
            cerr << "Error in selection of molecule type: molecule not accepted." << endl;
        }

        if (type_mol_2 == 1) // Type "Ethylene Oxide"
        {
            apm_mol_2 = 3;
            M2.push_back(14.026);
            M2.push_back(14.026);
            M2.push_back(15.999);
            totM2 = M2[0] + M2[1] + M2[2];
        }
        else if (type_mol_2 == 2) // Type "Methanol"
        {
            apm_mol_2 = 3;
            M2.push_back(15.034);
            M2.push_back(15.999);
            M2.push_back(1.008);
            totM2 = M2[0] + M2[1] + M2[2];
        }
        else if (type_mol_2 == 3) // Type "CO2"
        {
            apm_mol_2 = 3;
            M2.push_back(12.01);
            M2.push_back(15.999);
            M2.push_back(15.999);
            totM2 = M2[0] + M2[1] + M2[2];
        }
        else
        {
            cerr << "Error in selection of molecule type: molecule not accepted." << endl;
        }
}
double pbc_unwrap(float& boxlength, double& coord, double& prevcoord)
{
    double stepdiff = coord - prevcoord;
    double halfbox = boxlength / 2.0;
    while ( abs(stepdiff) >= halfbox )
    {
        if ( stepdiff > halfbox)
        {
            coord = coord - boxlength;     
        }
        else if ( stepdiff < -halfbox )
        {
            coord = coord + boxlength; 
        }
        else
        {
            // Do Nothing
        }
    stepdiff = coord - prevcoord;
    }
}

int main(int argc, char *argv[])
    int starttime = 0; endtime = 0, dumpfreq = 0, n_bin = 0, ntop = 0, nbeh = 0, n_mol_1 = 0, n_mol_2 = 0, type_mol_1 = 0; type_mol_2 = 0, selection_1 = 0, selection_2 = 0, pbc_unwrap_flag = 0; // Parameter File integers
    int apm_mol_1 =0, apm_mol_2 = 0, t_start = 0, t_end = 0, tot_time = 0, bin_time = 0;
    float boxlength = 0.0; // Parameter File Floats
    double totM1 = 0.0, totM2 = 0.0;
    vector<double> X, Y, Z, M1, M2, g_r;
    string bufferstring = "tmp", filename = "file.dat", delim = "tmp_"; // Parameter file strings

    if ( argc < 2 )
    {
        // argc should be 2 for correct execution
        cout<<"usage: "<< argv[0] <<" <filename>\n";
    }
    else 
    {
        ifstream param_file ( argv[1] );
        if (param_file.fail())
        {
            cerr << "Error: Parameter File Open Failed" << endl;
        }
        
        getline(param_file, bufferstring);
        param_file >> boxlength >> filename >> delim;
        getline(param_file, bufferstring);
        getline(param_file, bufferstring);
        param_file >> starttime >> endtime >> dumpfreq;
        getline(param_file, bufferstring);
        getline(param_file, bufferstring);
        param_file >> n_bin;
        getline(param_file, bufferstring);
        getline(param_file, bufferstring);
        param_file >> ntop >> nbeh;
        getline(param_file, bufferstring);
        getline(param_file, bufferstring);
        param_file >> n_mol_1 >> type_mol_1;
        getline(param_file, bufferstring);
        getline(param_file, bufferstring);
        param_file >> n_mol_2 >> type_mol_2;
        getline(param_file, bufferstring);
        getline(param_file, bufferstring);
        param_file >> selection_1 >> selection_2;
        getline(param_file, bufferstring);
        getline(param_file, bufferstring);
        param_file >> pbc_unwrap_flag;
        
        // Handles Molecular Identity
        mol_selection(type_mol_1, type_mol_2, apm_mol_1, apm_mol_2, M1, M2, totM1, totM2);
        // Calculates various time things.
        t_start = starttime / dumpfreq;
        t_end = endtime / dumpfreq;
        tot_time = t_end - t_start;
        bin_time = tot_time / n_bin;


        // Open the File
        ifstream input_file ( argv[1] );
        if (input_file.fail())
        { 
            cerr << "Error: Trajectory File Open Failed." << endl;
        }


        for (int t = 0; t < t_end; t++)
        {
            if (ntop != 0)
            {
                for (int m = 0; m < ntop; m++)
                {
                    getline(input_file, bufferstring);
                }
            }
            for (int m = 0; m < natoms; m++)
            {
                ID.push_back(0);
                X.push_back(0);
                Y.push_back(0);
                Z.push_back(0);
                input >> ID[index] >> X[index] >> Y[index] >> Z[index];
                if ( pbcwrapflag == 1 && t != 0 )
                {
                    pbc_unwrap(boxlength, X[index], X[index-natoms]);
                    pbc_unwrap(boxlength, Y[index], Y[index-natoms]);
                    pbc_unwrap(boxlength, Z[index], Z[index-natoms]);
                }
                index++;
            }
            getline(input, bufferstring);
            if (nbeh != 0)
            {
                for (int m = 0; m < nbeh; m++)
                {
                    getline(input, bufferstring);
                }
            }
        }
    }
  } 

    // Loop over blocks
    for (int b = 0; b < n_bin; b++)
    {
        // Loops over first atom
        for (int i = 0; i < n_mol_1; i++)
        {
            // Loops over second atom
            if (selection_2 > apm_mol_1)
            {
                for (int j = 0; j < n_mol_2; j++)
                {
                     
                }
            }
            else
            {
                int j = 0;
                while (j<i)
                {

                }
            }
        }
    }


    return 0;
}

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

double separation(double& val1, double& val2, float& boxlength)
{
    double DR = 0.0;
    DR = val1 - val2;
    //Added handling for PBC
    if (DR > boxlength/2)
    {
        DR = boxlength - DR;
    }
    else if (DR < -boxlength/2)
    {
        DR = boxlength + DR;
    }
    else
    {
        // Do Nothing
    }
    return DR;
}

int distance(double& dx, double& dy, double& dz, double& dr)
{
    double r = 0;
    int r_bin = 0;
    r = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));
    r_bin = floor(r/dr);
    return r_bin;
}

int main(int argc, char *argv[])
{
    int starttime = 0; endtime = 0, dumpfreq = 0, n_bin = 0, ntop = 0, nbeh = 0, n_mol_1 = 0, n_mol_2 = 0, type_mol_1 = 0; type_mol_2 = 0, selection_1 = 0, selection_2 = 0, pbc_unwrap_flag = 0; // Parameter File integers
    int apm_mol_1 =0, apm_mol_2 = 0, t_start = 0, t_end = 0, tot_time = 0, bin_time = 0, tot_r = 0, n_ref_mol = 0, n_sel_mol = 0;
    float boxlength = 0.0, dr = 0.0, r_end = 0.0; // Parameter File Floats
    double totM1 = 0.0, totM2 = 0.0;
    double dx = 0.0, dy = 0.0, dz = 0.0;
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
        param_file >> n_bin >> dr >> r_end;
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
        tot_r = r_end/dr;

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
        // Loop over blocks
        for (int b = 0; b < n_bin; b++)
        {
            int tmp_start_t =t_start + b*bin_time;
            int tmp_end_t = tmp_start_t + bin_time;
            vector<double> GofR(tot_r, 0);
            // loop over time
            for (int t = tmp_start_t; t < tmp_end_t; t++)
            {
                vector<int> rdf(tot_r, 0);
                //IF STATEMENT: Checks whether it is the rdf between same molecules, or different.
                if (selection_1 <= apm_mol_1 && selection_2 <= apm_mol_1) //CASE1 - First Mol Type - First Mol Type
                {
                    n_ref_mol = n_mol_1;
                    n_sel_mol = n_mol_1;
                    // Loops over first atom
                    for (int i = 0; i < n_mol_1; i++)
                    {
                        int j = 0;
                        while (j < i)
                        {
                            int atom1 = i*apm_mol_1+(selection_1-1);
                            int atom2 = j*apm_mol_1+(selection_2-1);
                            dx=separation(X[atom1],X[atom2],boxlength);
                            dy=separation(Y[atom1],Y[atom2],boxlength);
                            dz=separation(Z[atom1],Z[atom2],boxlength);
                            int tmp_rbin = distance(dx,dy,dz,dr);
                            rdf[tmp_rbin]+=1;
                            j+=1;
                        }
                    }
                }
                else if (selection_1 <= apm_mol_1 && selection_2 > apm_mol_1) //CASE2 - First Mol Type - Second Mol Type
                {
                    n_ref_mol = n_mol_1;
                    n_sel_mol = n_mol_2;
                    // Loops over first atom
                    for (int i = 0; i < n_mol_1; i++)
                    {      
                        for (int j = 0; j < n_mol_2; j++)
                        {
                            int atom1 = i*apm_mol_1 + (selection_1-1);
                            int atom2 = j*apm_mol_2 + (selection_2 - apm_mol_1 - 1);
                            dx=separation(X[atom1],X[atom2],boxlength);
                            dy=separation(Y[atom1],Y[atom2],boxlength);
                            dz=separation(Z[atom1],Z[atom2],boxlength);
                            int tmp_rbin = distance(dx,dy,dz,dr);
                            rdf[tmp_rbin]+=1;
                        }
                    }
                }
                else if (selection_1 > apm_mol_1 && selection_2 > apm_mol_1) //CASE3 - Second Mol Type - Second Mol Type
                {
                    n_ref_mol = n_mol_2;
                    n_sel_mol = n_mol_2;
                    // Loops over first atom
                    for (int i = 0; i < n_mol_1; i++)
                    {
                        for (int j = 0; j < n_mol_2; j++)
                        {
                            int atom1 = i*apm_mol_2 + (selection_1-apm_mol_1-1);
                            int atom2 = j*apm_mol_2 + (selection_2 - apm_mol_1 - 1);
                            dx=separation(X[atom1],X[atom2],boxlength);
                            dy=separation(Y[atom1],Y[atom2],boxlength);
                            dz=separation(Z[atom1],Z[atom2],boxlength);
                            int tmp_rbin = distance(dx,dy,dz,dr);
                            rdf[tmp_rbin]+=1;
                        }
                    }
                }
                else
                {
                    cerr << "There was an invalid selection" << endl;
                    exit(0);
                }
                for (int r = 0; r < tot_r; r++)
                {
                    GofR[r] += (rdf[r]/((double) n_ref_mol))*(1/(4*PI*pow(r*dr,2)*dr))*(pow(boxlength,3)/((double) n_sel_mol));
                }
            }
        }
        string gofr_out_fname = "gofr_" + delim + "_bin_" + to_string(b) + ".dat";
        ofstream gofrout;
        gofrout.open(gofr_out_fname);
        for (int r = 0; r < tot_r; r++)
        {
            GofR[r] = GofR[r]/((double) bin_time);
            gofrout << ((float) (r*dr+0.5*dr)) << GofR[r] << endl;
        }
    }
    return 0;
}

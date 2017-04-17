/*************************************************/
/*     Local Liquid Structure                    */
/*                                               */
/* Ezekiel Piskulich, Brian Laird, Ward Thompson */
/* University of Kanas, Lawrence                 */
/* Copyright 2017, All Rights Reserved           */
/*************************************************/

//This code requires an parameter file in the same directory as which this code is being run.
//An example of this parameter file is included as sample_input.txt
//Current Capabilities:
//- Pair Correlation Function
//- Coordination Number
//- Potential of Mean Force
//
//This calculation outputs a file for each block based on the paramter file.
//Usage: ./localliquidstructure.exe [parameterfile]

#include "main.h"

void mol_selection (int& type_mol_1, int& type_mol_2, int& apm_mol_1, int& apm_mol_2, vector<double>& M1, vector<double>& M2, double& totM1, double& totM2)
//This subroutine takes care of the molecule type selection.
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
            totM1 = M1[0] + M1[1] + M1[2];
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

/*
double pbc_unwrap(float& boxlength, double& coord, double& prevcoord)
//This subroutine takes care of all of the unwrapping.
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
*/

double separation(double& val1, double& val2, float& boxlength)
//This calculates the difference between two values, and then modifies based on periodic boundaries. 
{
    double DR = 0.0;
    DR = val1 - val2;
    //Added handling for PBC
    if (DR > boxlength/2.0)
    {
        DR = boxlength - DR;
    }
    else if (DR < -boxlength/2.0)
    {
        DR = boxlength + DR;
    }
    else
    {
        // Do Nothing
    }
    return DR;
}

int distance(double& dx, double& dy, double& dz, float& dr)
{
    double r = 0;
    int r_bin = 0;
    r = sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));
    r_bin = floor(r/dr);
    return r_bin;
}

double trapezoidal_rule(int& a, int& b, int& N, vector<double>& f_x)
{
    double trap_out = 0.0;
    double h=(b-a)/(double) N;
    for (int i=a; i<=b; i++)
    {
        if (i == a)
        {
            trap_out+=0.5*f_x[i];
        }
        else if (i == b)
        {
            trap_out+=0.5*f_x[i];
        }
        else
        {
            trap_out+=f_x[i];
        }
    }
    trap_out*=h;
    return trap_out;
}


int main(int argc, char *argv[])
{
    int starttime = 0, endtime = 0, dumpfreq = 0, n_bin = 0, ntop = 0, nbeh = 0, n_mol_1 = 0, n_mol_2 = 0, type_mol_1 = 0, type_mol_2 = 0, selection_1 = 0, selection_2 = 0, pbcwrapflag = 0; // Parameter File integers
    int apm_mol_1 =0, apm_mol_2 = 0, t_start = 0, t_end = 0, tot_time = 0, bin_time = 0, tot_r = 0, n_ref_mol = 0, n_sel_mol = 0, natoms = 0, CASE = 0;
    float boxlength = 0.0, dr = 0.0, r_end = 0.0; // Parameter File Floats
    double totM1 = 0.0, totM2 = 0.0, number_rho = 0.0, temperature;
    double dx = 0.0, dy = 0.0, dz = 0.0;
    double kb = 0.0019872041; //kcal/mol K
    double norm_const = 0.0;
    vector<double> ID, X, Y, Z, M1, M2;
    string bufferstring = "tmp", filename = "file.dat", delim = "tmp_"; // Parameter file strings

    if ( argc < 2 )
    {
        // argc should be 2 for correct execution
        cout<<"Proper Usage: "<< argv[0] <<" [filename]\n";
    }
    else 
    {
        // Reads in the parameter file
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
        param_file >> pbcwrapflag >> temperature;
            
        // Handles Molecular Identity
        mol_selection(type_mol_1, type_mol_2, apm_mol_1, apm_mol_2, M1, M2, totM1, totM2);
        // Calculates various time things.
        t_start = starttime / dumpfreq;
        t_end = endtime / dumpfreq;
        tot_time = t_end - t_start;
        bin_time = tot_time / n_bin;
        tot_r = r_end/dr;
        natoms = n_mol_1*apm_mol_1 + n_mol_2*apm_mol_2;
        
        //Handles Which Case is the correct one:
        if (selection_1 <= apm_mol_1 && selection_2 <= apm_mol_1) // Mol Type 1 - Mol Type 1
        {
            CASE=1;
            n_ref_mol = n_mol_1;
            n_sel_mol = n_mol_1-1;
        }
        else if (selection_1 <= apm_mol_1 && selection_2 > apm_mol_1) // Mol Type 1 - Mol Type 2
        {
            CASE=2;
            n_ref_mol = n_mol_1;
            n_sel_mol = n_mol_2;
        }
        else if (selection_1 > apm_mol_1 && selection_2 > apm_mol_1) //Mol Type 2 - Mol Type 2
        {
            CASE=3;
            n_ref_mol = n_mol_2;
            n_sel_mol = n_mol_2;
        }
        else
        {
            cerr << "Error: Atom Selection out of Proper Scope" << endl;
            throw exception();
        }
        //Calculates some necessary RDF constants (see Allen and Tildesley)
        number_rho=n_ref_mol/pow(boxlength,3);
        norm_const=4.0*M_PI*number_rho/3.0;

        //Outputs a check of the parameter file.
        cout << "boxlength filename delim" << endl;
        cout << boxlength << " " << filename << " " << delim << endl;
        cout << "starttime endtime dumpfreq" << endl;
        cout << starttime << " " << endtime << " " << dumpfreq << endl;
        cout << "n_bin dr r_end" << endl;
        cout << n_bin << " " << dr << " " << r_end << endl;
        cout << "ntop nbeh" << endl;
        cout << ntop << " " << nbeh << endl;
        cout << "n_mol_1 type_mol_1" << endl;
        cout << n_mol_1 << " " << type_mol_1 << endl;
        cout << "n_mol_2 type_mol_2" << endl;
        cout << n_mol_2 << " " << type_mol_2 << endl;
        cout << "selection_1 selection_2" << endl;
        cout << selection_1 << " " << selection_2 << endl;
        cout << "pbcwrapflag temperature" << endl;
        cout << pbcwrapflag << " " << temperature << endl;
        cout << "--------Computed Quantities---------" << endl;
        cout << "t_start t_end" << endl;
        cout << t_start << " " << t_end << endl;
        cout << "tot_time bin_time" << endl;
        cout << tot_time << " " << bin_time << endl;
        cout << "tot_r natoms" << endl;
        cout << tot_r << " " << natoms << endl;

        // Opens the Trajectory File
        ifstream input_file ( filename );
        if (input_file.fail())
        { 
            cerr << "Error: Trajectory File Open Failed." << endl;
            throw exception();
        }
        
        // Reads the trajectory into the ID, X, Y, Z vectors!
        int index=0;
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
                input_file >> ID[index] >> X[index] >> Y[index] >> Z[index];
                index++;
            }
            getline(input_file, bufferstring);
            if (nbeh != 0)
            {
                for (int m = 0; m < nbeh; m++)
                {
                    getline(input_file, bufferstring);
                }
            }
        }
        int count = 0;
        // Loop over blocks
        for (int b = 0; b < n_bin; b++)
        {
            cout << "Calculating block " << b << endl;
            // Calculates Starting and Ending times for the block
            int tmp_start_t =t_start + b*bin_time;
            int tmp_end_t = tmp_start_t + bin_time;
            // Defines the Vectors that Store the GOFR, NOFR, and AOFR data.
            vector<double> GofR(tot_r, 0);
            vector<double> NofR(tot_r, 0);
            vector<double> Ntmp(tot_r, 0);
            vector<double> AofR(tot_r, 0);
            vector<int> rdf(tot_r, 0);

            // TIME LOOP
            for (int t = tmp_start_t; t < tmp_end_t; t++)
            {
                //IF STATEMENT: Checks whether it is the rdf between same molecules, or different.
                if (CASE==1) //CASE1 - First Mol Type - First Mol Type
                {
                    // Loops over first atom
                    for (int i = 0; i < n_mol_1-1; i++)
                    {
                        // Loops over the second atom
                        for (int j = i+1; j<n_mol_1; j++)
                        {
                            // Calculates the atom indices of both atom based on the frame
                            int atom1 = i*apm_mol_1+(selection_1-1)+t*natoms;
                            int atom2 = j*apm_mol_1+(selection_2-1)+t*natoms;
                            // Calculates separation components and bins distance
                            dx=separation(X[atom1],X[atom2],boxlength);
                            dy=separation(Y[atom1],Y[atom2],boxlength);
                            dz=separation(Z[atom1],Z[atom2],boxlength);
                            int tmp_rbin = distance(dx,dy,dz,dr);
                            // Increments the count by two.
                            if (tmp_rbin<=tot_r)
                            {
                                rdf[tmp_rbin]+=2;
                            }
                        }
                    }
                }
                else if (CASE==2) //CASE2 - First Mol Type - Second Mol Type
                {
                    // Loops over first atom
                    for (int i = 0; i < n_mol_1; i++)
                    {      
                        for (int j = 0; j < n_mol_2; j++)
                        {
                            int atom1 = i*apm_mol_1 + (selection_1-1)+t*natoms;
                            int atom2 = j*apm_mol_2 + (selection_2 - apm_mol_1 - 1)+t*natoms;
                            dx=separation(X[atom1],X[atom2],boxlength);
                            dy=separation(Y[atom1],Y[atom2],boxlength);
                            dz=separation(Z[atom1],Z[atom2],boxlength);
                            int tmp_rbin = distance(dx,dy,dz,dr);
                            if (tmp_rbin<=tot_r)
                            {
                                rdf[tmp_rbin]+=1;
                        
                            }
                        }
                    }
                }
                else if (CASE==3) //CASE3 - Second Mol Type - Second Mol Type
                {
                    // Loops over first atom
                    for (int i = 0; i < n_mol_1; i++)
                    {
                        for (int j = 0; j < n_mol_2; j++)
                        {
                            int atom1 = i*apm_mol_2 + (selection_1-apm_mol_1-1)+t*natoms;
                            int atom2 = j*apm_mol_2 + (selection_2 - apm_mol_1 - 1)+t*natoms;
                            dx=separation(X[atom1],X[atom2],boxlength);
                            dy=separation(Y[atom1],Y[atom2],boxlength);
                            dz=separation(Z[atom1],Z[atom2],boxlength);
                            int tmp_rbin = distance(dx,dy,dz,dr);
                            if (tmp_rbin<=tot_r)
                            {
                                rdf[tmp_rbin]+=1;
                            }
                        }
                    }
                }
                else
                {
                    cerr << "Error: There was an invalid selection that somehow got missed earlier." << endl;
                    throw exception();
                }
                /*
                // Calculates the Normalized Cumulative Distribution aka NofR
                int SUMRDF = 0;
                for (int r = 0; r < tot_r; r++)
                {
                    SUMRDF += rdf[r]; 
                    NofR[r] += SUMRDF;
                }
                */
            }

            //Defines the filenames for the three output files (per block)
            string gofr_out_fname = "gofr_" + delim + "_bin_" + to_string(b) + ".dat";
            string nofr_out_fname = "nofr_" + delim + "_bin_" + to_string(b) + ".dat";
            string aofr_out_fname = "aofr_" + delim + "_bin_" + to_string(b) + ".dat";
            //Opens the streams, and then opens the files
            ofstream gofrout,nofrout,aofrout;
            gofrout.open(gofr_out_fname), nofrout.open(nofr_out_fname), aofrout.open(aofr_out_fname);
            //Defines the number density of the simulation cell.
            for (int r = 0; r < tot_r; r++)
            {
                double inner_r=r*dr;
                double outer_r=inner_r+dr;
                double ideal_n=norm_const * ( pow(outer_r,3) - pow(inner_r,3) );
                GofR[r] = rdf[r] / ((double) bin_time * (double) n_ref_mol * ideal_n);
                
                // Calculates the Potential of Mean Force
                if (GofR[r] > 0)
                {
                    AofR[r] = -kb*temperature*log(GofR[r]);
                }
                else
                {
                    AofR[r]=0;
                }
                Ntmp[r]=GofR[r]*pow(outer_r,2);
            }
            int a=0;
            for (int b=a+1; b<tot_r; b++)
            {
                int npoints=b; //Only works since a=0;
                NofR[b]=dr*4.0*M_PI*number_rho*trapezoidal_rule(a,b,npoints,Ntmp);
            }

            for (int r = 0; r < tot_r; r++)
            {
                gofrout << ((float) (r*dr+0.5*dr)) << " " << GofR[r] << endl;
                nofrout << ((float) (r*dr+0.5*dr)) << " " << NofR[r] << endl;
                aofrout << ((float) (r*dr+0.5*dr)) << " " << AofR[r] << endl;
            } 
        }
        
    }
    cout << "Exit was a success!" << endl;
    return 0;
}

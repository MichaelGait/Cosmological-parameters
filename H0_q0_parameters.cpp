//Analysis of the SNIa UNION21 catalogue  (Suzuki et al. 2012)
//Estimation of the present Hubble's parameter H0 and deceleration parameter q0
//Samples used: 580 supernovae Ia

#include<iostream>
#include<iomanip>
#include<fstream>
#include<cstdlib>
#include<cmath>
#include<ctime>
#include<vector>

using namespace std;

const float c = 3.0*pow(10,5); //light speed [km/sec]
float H1,H2,dH;
float q1,q2,dq;

//SN's effective luminosity distance for fixed parameters (H0,q0)
//Good approximation for SN with z <= 0.3
float EffectiveLuminosityDistance(float H0, float q0, float z)
{
    float d = (c/H0)*(z + 0.5*(1-q0)*z*z);
    return d;
}

//SN's theoretical distance modulus for fixed parameters (H0,q0)
float DistanceModulus(float d)
{
    float mu = 5*log10(d) + 25.0; //d [Mpc]
    return mu;
}

//i-term of chi_squared sum
float chi_squared(float muExp, float muTheo, float sigma)
{
    return ((muExp - muTheo)*(muExp - muTheo))/(sigma*sigma);
}

int main()
{
    cout << fixed << setprecision(4);
    cout << "Analysis of the SNIa UNION21 catalogue  (Suzuki et al. 2012)\n";
    ifstream fp1;
    fp1.open("Union21.txt");
    if (!fp1.is_open())
    {
        cout << "File not found. Aborting...\n";
        exit(EXIT_FAILURE);
    }

    vector<string> SN_name; string name;
    vector<float>  z;       float redshift;
    vector<float>  muExp;   float dModulusExp;
    vector<float>  sigma;   float s;
                            float omit;

    while (fp1 >> name >> redshift >> dModulusExp >> s >> omit)
    {
        SN_name.push_back(name);
        z.push_back(redshift);
        muExp.push_back(dModulusExp);
        sigma.push_back(s);
    }
    int fileRows = SN_name.size();
    cout << "SN samples used : " << fileRows << endl;

    ofstream fp2;
    fp2.open("H0_q0.txt");
    fp2 << fixed << setprecision(3);
    fp2 << "H0\t\tq0\t\tx^2\n\n";

    float z_min = z[0];
    float z_max = z[0];
    for (int i = 0; i < fileRows; ++i)
    {
        if (z[i] < z_min) z_min = z[i];
        if (z[i] > z_max) z_max = z[i];
    }
    cout << "z_min = " << z_min << "\t" << "z_max = " << z_max << endl;

    //H0-q0 grid's dimensions and steps
    cout << "Define H0-q0 grid\n";
    cout << "H1 = "; cin >> H1;
    cout << "H2 = "; cin >> H2;
    cout << "dH = "; cin >> dH;
    cout << "q1 = "; cin >> q1;
    cout << "q2 = "; cin >> q2;
    cout << "dq = "; cin >> dq;
    cout << "Calculating data. Please wait...\n";

    float min_chi_squared = 9999999.0;
    float result_q0 = 9999999.0, result_H0 = 9999999.0;
    float muTheo; //theoretical distance modulus

    clock_t begin, end;
    begin = clock();
    for (float H0 = H1; H0 <= H2; H0 += dH)
    {
        for (float q0 = q1; q0 <= q2; q0 += dq)
        {
            float sum = 0.0; //chi_squared sum
            bool flag = true;
            for (int i = 0; i < fileRows; i++)
            {
                if (z[i] <= 0.3)
                {
                    if (EffectiveLuminosityDistance(H0,q0,z[i]) > 0.0)
                    {
                        float d = EffectiveLuminosityDistance(H0,q0,z[i]);
                        muTheo = DistanceModulus(d);
                        sum += chi_squared(muExp[i],muTheo,sigma[i]);
                    }
                    else //DistanceModulus() returns NaN because of log10()
                    {
                        flag = false;
                        break;
                    }
                }
            }
            if (flag == true)
            {
                fp2 << H0 << "\t\t" << q0 << "\t\t" << sum << endl;
                if ((sum < min_chi_squared))
                {
                    min_chi_squared = sum;
                    result_q0 = q0;
                    result_H0 = H0;
                }
            }
        }
    }
    end = clock();
    float dt = (float)(end-begin)/CLOCKS_PER_SEC;
    cout << "Estimated calculation time >> dt = " << dt << " sec" << endl;
    fp2 << "Estimation >> q0 = " << result_q0 << "  H0 = " << result_H0
    << "  min_chi_squared = " << min_chi_squared << endl;
    cout << "Estimation >> q0 = " << result_q0 << "  H0 = " << result_H0
    << "   min_chi_squared = " << min_chi_squared << endl;
    cout << "Data were written in 'H0_q0.txt'\n";
    SN_name.clear();
    z.clear();
    muExp.clear();
    sigma.clear();
    fp1.close();
    fp2.close();
    return 0;
}

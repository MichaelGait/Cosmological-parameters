//Analysis of the SNIa UNION21 catalogue  (Suzuki et al. 2012)
//Estimation of the present density parameters (Omega_matter and Omega_Lambda)
//Samples used: 580 supernovae Ia

#include<iostream>
#include<iomanip>
#include<fstream>
#include<cstdlib>
#include<ctime>
#include<cmath>
#include<vector>

using namespace std;

const float c = 3.0*pow(10,5); //light speed [km/sec]
const float w = -1.0; //dark energy
const float H0 = 70.022; //known from previous calculations

float LuminosityDistance(float z, float Integral)
{
    float d = (1+z)*(c/H0)*Integral;
    return d;
}

//SN's theoretical distance modulus for fixed parameters (H0,q0)
float DistanceModulus(float d)
{
    float mu = 5*log10(d) + 25.0; //d [Mpc]
    return mu;
}

//Assuming flat Universe (Omega_k = 0)
//Assuming no radiation (Omega_rad = 0)
//Assuming cosmological constant (Omega_lambda > 0)
//Assuming baryonic luminous + dark matter (Omega_matter > 0)
float E(float z, float Omega_m, float Omega_L)
{
    return 1.0/sqrt(Omega_m*pow(1+z,3) + Omega_L*pow(1+z,3.*(1+w)));
}

//Simpson's integration of the function 1/E(z) for fixed Omega parameters
float Integral(const float a, const float b, int N, float Omega_m, float Omega_L)
{
    float h = abs((b-a)/((float)(N-1))); //N must be odd!!
    float I = 0.0;
    unsigned int i = 1;
    for (float z = (a+h); z < b; z += h)
    {
        if (i%2 == 1)
            I = I + 4*E(z,Omega_m,Omega_L);
        else
            I = I + 2*E(z,Omega_m,Omega_L);
        ++i;
    }
    I = (h/3)*(E(a,Omega_m,Omega_L) + I + E(b,Omega_m,Omega_L));
    return I;
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
    fp2.open("Omega_parameters.txt");
    fp2 << fixed << setprecision(3);
    fp2 << "Omega_m\t\tOmega_L\t\tx^2\n\n";

    //Omega_m - Omega_L grid's dimensions and steps
    cout << "Define Omega steps\n";
    float dOmega_m,dOmega_L;
    cout << "dOmega_m = "; cin >> dOmega_m;
    cout << "dOmega_L = "; cin >> dOmega_L;
    const float Omega_m1 = 0.0, Omega_m2 = 1.0;
    const float Omega_L1 = dOmega_L, Omega_L2 = 1.0;

    cout << "Calculating data. Please wait...\n";

    const int N = 51; //integration points for the Simpson's method
    float min_chi_squared = 9999999.0;
    float result_Omega_m0 = 999999.0, result_Omega_L0 = 999999.0;
    float muTheo = 0.0;

    clock_t begin, end;
    begin = clock();
    for (float Omega_m0 = Omega_m1; Omega_m0 <= Omega_m2; Omega_m0 += dOmega_m)
    {
        for (float Omega_L0 = Omega_L1; Omega_L0 <= Omega_L2; Omega_L0 += dOmega_L)
        {
            float sum = 0.0;
            for (int i = 0; i < fileRows; i++)
            {
                float I = Integral(0.0,z[i],N,Omega_m0,Omega_L0);
                float d = LuminosityDistance(z[i],I);
                muTheo = DistanceModulus(d);
                sum += chi_squared(muExp[i],muTheo,sigma[i]);
            }
            fp2 << Omega_m0 << "\t\t" << Omega_L0 << "\t\t" << sum << endl;
            if ((sum < min_chi_squared))
            {
                min_chi_squared = sum;
                result_Omega_m0 = Omega_m0;
                result_Omega_L0 = Omega_L0;
            }
        }
    }
    end = clock();
    float dt = (float)(end-begin)/CLOCKS_PER_SEC;
    cout << "Estimated calculation time >> dt = " << dt << " sec" << endl;
    fp2 << "Estimation >> Omega_m0 = " << result_Omega_m0
    << "  Omega_L0 = " << result_Omega_L0 << "  min_chi_squared = " << min_chi_squared << endl;
    cout << "Estimation >> Omega_m0 = " << result_Omega_m0
    << "  Omega_L0 = " << result_Omega_L0 << "  min_chi_squared = " << min_chi_squared << endl;
    cout << "Data were written in 'Omega_parameters.txt'\n";
    SN_name.clear();
    z.clear();
    muExp.clear();
    sigma.clear();
    fp1.close();
    fp2.close();
    return 0;
}

#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */


# include <cmath>
# include <cctype>
# include <string>
# include <iostream>
# include <math.h>
# include <stdio.h>
# include <fstream>
# include <vector>
using namespace std;




//Structure containing the cosmological parameters
struct s_cosmological_parameters{
  double cosmological_redshift;
  double Om_matter;     //Energy density of total matter in units of the the critical density
  double Om_cdm;  //Energy density of cold_dark_matter in units of the the critical density
  double Om_radiation;  //Energy density of radiation in units of the the critical density
  double Om_baryons;    //Energy density of baryons in units of the the critical density
  double Om_vac;        //Energy density of dark energy in units of the the critical density
  double Om_k;          //'Energy density' of curvature in units of the the critical density 
  double f_baryon;
  double Hubble;        //Hubble parameter in units of h km / s / Mpc 
  double hubble;        //dimensionless Hubble parameter
  double w_eos;         //equation of state of dark energy
  double N_eff;         //Effective number of relativistic particles
  double sigma8;        //RMS of matter fluctuations at R = 8 Mpc/h
  double A_s;            // Amplitude of primordial ower spectrum
  double n_s;           //primordial spectral index
  double alpha_s;           //primordial spectral index
  double Tcmb;          //CMB temperature
  double Mabs;          //Something related to a magnitude, used for vmax stuff
  double mlim;          //Something related to a magnitude, used for vmax stuff
  double RR;
  double Delta_SO;
  bool use_wiggles;
  double e_index_a;  // constant factor in the e-correction
  double e_index_b;  // constant factor in the e-correction
  double e_index_c;  // constant factor in the e-correction
  double e_index_d;  // constant factor in the e-correction
  double e_index_zstar;  // constant factor in the e-correction
  double beta_rsd;

  double k_index;  // constant factor in the k-correction
  double d_index; // related to the density evolution
  double GAL_BIAS;
  double alpha_BIAS;
  // aca van variables que estan asociadas a integracion en k, o
  // que vana servir para poner mas parámetros y asi pasar una sola
  // estructura a las rutinas de integracion
  double kmin_int;

  double kmax_int;
  string mass_function_fit;
  string halo_mass_bias_fit;
  string density_profile;


  double h4; // values of integrals used in the HF for intergaration wrt to z of P(k,z)
  double h2;

  double M_max_mf;
  double M_min_mf;
  double kmin_ps;
  double kmax_ps;
 
  double aux_var1; //use it for m
  double aux_var2; //use ot for z
  double aux_var3; //use ot for r
  double aux_var4; //use ot for k

  double Mnl;
  double knl;
  double rnl;
  
  double knl_hf;
  double rnl_hf;
  double kstar;
  double M_min_effective;
  double M_max_effective;

  double A_PS;
  double Q_PS;


  double Amc; //Amplitude of the non linear correctionto P(k) in PT

  // para cosas dependientes de z
  double critical_density;
  double density_contrast_top_hat;
  double mean_matter_density;
  double growth_factor;
  double pk_normalization;
  
  int n_points_mass;
  int n_points_mass_integration;

  vector<double> v_mass;
  vector<double> v_sigma_mass;
  vector<double> v_mass_function;
  vector<double> v_halo_mass_bias;
  vector<double> v_k_ps;
  vector<double> v_lin_power_spectrum;
  vector<double> v_nl_power_spectrum;
  vector<double> v_density_profile_k;
  vector<double> v_density_profile_r;
  vector<double> v_galaxy_power_spectrum_1h_ss;
  vector<double> v_galaxy_power_spectrum_1h_sc;
  vector<double> v_galaxy_power_spectrum_2h;
  bool use_K_correction;
  bool use_e_correction;


  vector<double>rv;
  vector<double>gv;
  vector<double>zv;
  vector<double>trv;




  int hod_model;
  double mmin_hod;
  double alpha_hod;
  double scatter_hod;
  double muno_hod;
  double coef_concentration;
  double coef_concentration_amp;

};





struct s_hods_par{
  int hod_model;
  double mmin;
  double alpha;
  double scatter;
  double muno;
};


struct s_astrophysical_parameters{
  double A_gas;
  double B_gas;
  double mstar;
  double sigma_red;
  double sigma_ln;
  double missing_flux;
};


struct matrices{
  vector<double>Cmed;
  vector<double>Dmed;
  vector<vector<double> >R;
  vector<vector<double> >V;  
  vector<vector<double> >N;
  vector<vector<double> >iCov;  
  vector<vector<double> >nCov;
  vector<vector<double> >Step;
  double det_matrix;
};



// Structure used in the Cl
struct params_clth{
  double r;
  double k;
  int l;
  string wtype;
  double zmax_bin;
  double zmin_bin;
  double rmax_bin;
  double rmin_bin;

  double k_min_integration;
  double k_max_integration;
  double sigma_p_errors;
  string pdf_zerrors;
  double zaux;
  vector<double>pk;
  vector<double>kv;
  vector<double>Fkernel;
  vector<double>rv;
  vector<double>zv;
  vector<double>gv;
  vector<double>bias_zv;
  vector<double>gfv;
  vector<double>Hv;
  vector<double>dn_photo;
  vector<double>dn_spect;
  vector<double>dz;
  vector<double>klnv;
  vector<double>sigma_mass;
  vector<double>M_nl;
  vector<double>HF_i4;
  vector<double>HF_i2;


  vector< vector<double > > sBessel;
  vector<double> sxBessel;
  string redshift;
  double pk_normalization;
  double n_s; // not used
  bool use_non_linear_pk;
  string type_of_nl_power;
};




// This structure is created in order to avoid missconfusing of information among threads when using
// open MP.
template<typename T>
struct s_aux{
  s_cosmological_parameters *scp_a;
  params_clth *s_clth;
  double raux;
  int laux;
  double kaux;
  double zaux;
  vector<double>v_aux;
  vector<double>XX;
  vector<double>WW;
  vector<double>XX_mu;
  vector<double>WW_mu;
  vector<double>XX_z;
  vector<double>WW_z;
  T Ps;
};


// Structure used in the passage from Cl to Ps to spped up things
struct A1{
    s_cosmological_parameters *s_cp;
    vector<double>MASS;
    vector<double>MASS_FUNCTION;
    vector<double>MASS_BIAS;
    double aux_k;
    double aux_z;
    double aux_m;

};

struct experiments{
  vector<vector<double> > acc_par1;
  vector<double>  weight_par1;
  vector<vector<double> > acc_par2;
  vector<double>  weight_par2;
  vector<vector<double> > acc_par3;
  vector<double>  weight_par3;
  vector<vector<double> > acc_par4;
  vector<double>  weight_par4;
  vector<vector<double> > acc_par5;
  vector<double>  weight_par5;
  vector<vector<double> > acc_par6;
  vector<double>  weight_par6;
};

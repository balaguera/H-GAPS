#define nbess 200
#define NUM_ITER_MAP2ALM 0   /* Number of iterations for the map2alm_iter HealPix routine*/
#define EPSILON_MK 1.0
// ######################################################################
// ######################################################################
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
// ######################################################################
// ######################################################################

# include <unistd.h>
# include <ctime>
# include <cmath>
# include <cctype>
# include <string>
# include <iostream>
# include <math.h>
# include <stdio.h>
# include <fstream>
# include <iomanip>
# include <sstream>
# include <cassert>
# include <vector>
# include <complex>
# include <stdexcept>

# include <gsl/gsl_sf_coupling.h>
# include <gsl/gsl_sf_legendre.h>
# include <gsl/gsl_sf_gamma.h>
# include <gsl/gsl_sf.h>
# include <gsl/gsl_matrix.h>
# include <gsl/gsl_linalg.h>
# include <gsl/gsl_rng.h>
# include <gsl/gsl_randist.h>
# include <gsl/gsl_vector.h>
# include <gsl/gsl_eigen.h>
# include <gsl/gsl_blas.h>
# include <gsl/gsl_errno.h>
# include <gsl/gsl_integration.h>
# include <omp.h>

# include "../../Lib/Healpix_3.31/src/cxx/generic_gcc/include/alm.h"
# include "../../Lib/Healpix_3.31/src/cxx/generic_gcc/include/alm_fitsio.h"
# include "../../Lib/Healpix_3.31/src/cxx/generic_gcc/include/healpix_map.h"
# include "../../Lib/Healpix_3.31/src/cxx/generic_gcc/include/healpix_map_fitsio.h"
# include "../../Lib/Healpix_3.31/src/cxx/generic_gcc/include/healpix_data_io.h"
# include "../../Lib/Healpix_3.31/src/cxx/generic_gcc/include/healpix_base.h"
# include "../../Lib/Healpix_3.31/src/cxx/generic_gcc/include/alm_powspec_tools.h"
# include "../../Lib/Healpix_3.31/src/cxx/generic_gcc/include/alm_healpix_tools.h"
# include "../../Lib/Healpix_3.31/src/cxx/generic_gcc/include/xcomplex.h"
# include "../../Lib/Healpix_3.31/src/cxx/generic_gcc/include/sort_utils.h"
# include "../../Lib/Healpix_3.31/src/cxx/generic_gcc/include/sharp_cxx.h"
# include "FileManager.h"
# include "Numerics.h"

using namespace std;

//#define fac M_PI/180.0
#define pi_ 3.141592653589793238462643383279502884197
#define fac 3.141592653589793238462643383279502884197/180.0
#define no_zs -999.0





// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void message(){
  time_t rawtime;
  time (&rawtime);
  std::cout<<"*******************************************************"<<endl;
  std::cout<<"*******************************************************"<<endl;
  std::cout<<"ESTIMATES OF ANGULAR POWER SPECTRUM                   *"<<endl;
  std::cout<<"v1.2                                                  *"<<endl;
  std::cout<<"abalant@R3 2014-2016                                  *"<<endl;
  std::cout<<"*******************************************************"<<endl;
  std::cout<<"*******************************************************"<<endl;
  std::cout<<"Starting time and date"<<std::endl;
  std::cout<<ctime (&rawtime)<<std::endl;
  std::cout<<"*******************************************************"<<endl;
}



void usage(string s){
  cout<<CYAN<<" ********************************************************************************"<<endl;
  cout<<CYAN<<" ********************************************************************************"<<endl;
  cout<<RED<<"  Clgal a code to measure angular power spectrum of cosmological mass tracers"<<endl;
  cout<<RED<<"  Usage "<<s<<" [-option] [argument]"<<endl;
  cout<<RED<<"  Options:     -h for more information "<<endl;
  cout<<RED<<"               -a for information on the author"<<endl;
  cout<<RED<<"               -p parameter_file.ini to execture with input parameter file "<<endl;
  cout<<CYAN<<" ********************************************************************************"<<endl;
  cout<<CYAN<<" ********************************************************************************"<<endl;}

void author(){
  cout<<CYAN<<"                                                                                 "<<endl;
  cout<<" ***************************************************************************************"<<endl;
  cout<<" ***************************************************************************************"<<endl;
  cout<<" Copyright. 2017. Andres Balaguera Antolinez"<<endl;
  cout<<" Code developed for the power spectrum analysis of the 2MPZ sample"<<endl;
  cout<<" This code is public. If you want to use it for your research,     "<<endl;
  cout<<" contact the author at balaguera@iac.es"<<endl;
  cout<<" ***************************************************************************************"<<endl;
  cout<<" ***************************************************************************************"<<endl;
  cout<<CYAN<<"                                                                                 "<<endl;
}



void message_time(time_t start_all){
    time_t end;
    time(&end);
    std::cout<<RED<<"=================================="<<std::endl;
    std::cout<<"=================================="<<RESET<<std::endl;
    double lapse=difftime(end,start_all);
    std::cout<<CYAN<<"Time elapsed: "<<lapse<<"  secs \r"<<RESET<<std::endl;
    time(&start_all);
    std::cout<<RED<<"=================================="<<std::endl;
    std::cout<<"=================================="<<RESET<<std::endl;
}




// ######################################################################
// ######################################################################
// ######################################################################
void gal2equ(double b, double l, double *ra, double *dec){
    double alpha_p=fac*192.859508;                        /*Right ascention of galactic north pole in radians*/
    double delta_p=fac*27.128336;                         /*Declination of galactic north pole in radians*/
    double l_n   =fac*122.932;                           /*Galactic longitude of the celestial pole, in radians*/
    double delta=sin(b*fac)*sin(delta_p)+cos(b*fac)*cos(delta_p)*cos(fac*l-l_n);
    *dec=asin(delta);
    double alpha1=sin(fac*l-l_n);
    double alpha2=cos(fac*l-l_n)*sin(delta_p)-tan(fac*b)*cos(delta_p);
    *ra=atan2(alpha1, alpha2)+12.25;
}



// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// Note that when we do averages over the m-modes,
// we do the following:
// (Sum_m Jlm) / (2l+1) from -l to +l is splitted
// in Jl0/(2l+1)+ sum_(m=1) Jlm/(l+0.5)
// Since we put all in the same loop, we make theinf (m==) to divide things by 2
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

#ifndef _CL_FUNCTIONS_
#define _CL_FUNCTIONS_


/**
 *  @class Cl_FUNCTIONS.h
 *  @brief The class Cl_FUNCTIONS_
 *
 */


class Cl_FUNCTIONS{
  
 protected:

 private:
  //////////////////////////////////////////////////////////
  /**
  * @brief Number of columns in galaxy cat
  */
  int n_columns_gal;
  //////////////////////////////////////////////////////////
  /**
  * @brief Number of columns in random
  */
  int n_columns_ran;

  //////////////////////////////////////////////////////////
  /**
   * @brief Number of columns in the mask
   */
  int n_columns_mask;
  
  //////////////////////////////////////////////////////////
  /**
   * @brief Object of type pointing, used in HealPix operations
  */
  pointing point;
  //////////////////////////////////////////////////////////
  /**
   * @brief RMS of galaxies in pixels
   */
  double rms_nran;
  //////////////////////////////////////////////////////////
  /**
   * @brief Mean number of random objects in pixels
   */
  double mean_number_randoms_pix;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of rings, depends on nside
   */
  int Nrings;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of redshift bins, set default zero
   */
  int nzbins;
  //////////////////////////////////////////////////////////
  /**
   * @brief Index identifying redhsift bin
   */
  int IZ;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector containing the minimum redshifts in the redshift bins
   */
  vector<double> z_min;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector containing the maximum redshifts in the redshift bins
   */
  vector<double> z_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the Wigner 3J symbols, used
   * to compute the mixing matrix Rll in it-s approximated
   * version,
   *  in the VLS
   */

  vector<vector<vector<double> > > Wigner3J;
  //////////////////////////////////////////////////////////






public:
  
  //////////////////////////////////////////////////////////
  /**
   * @brief Constructor
   */
  Cl_FUNCTIONS(){
     std::cout.precision(6);
     std::cout.setf(ios::showpoint);
     std::cout.setf(ios::scientific);
  }

  Cl_FUNCTIONS(string par_file){
     read_pars(par_file);
     std::cout.precision(6);
     std::cout.setf(ios::showpoint);
     std::cout.setf(ios::scientific);
  }


  //////////////////////////////////////////////////////////
  /**
   * @brief Default destructor
   */
  ~Cl_FUNCTIONS(){}



  //////////////////////////////////////////////////////////
  /**
  *  @brief Galaxy catalog
  */
  vector<double>   prop;
  //////////////////////////////////////////////////////////
  /**
  *  @brief Random catalog
  */
  vector<double>  prop_r;



  //////////////////////////////////////////////////////////
  /**
  *  @brief Read input cats
  */

  void read_input_cats(string cat, string file);


  //////////////////////////////////////////////////////////
  /**
  *  @brief Get redshift bins
  */
  void get_zbins_same_ngal(int, int, double, double,vector< vector<double> >&);
  //////////////////////////////////////////////////////////

  /**
  *    @brief Inpput/Output
  */
  FILE_MANAGER<int> Fmi;
  //////////////////////////////////////////////////////////
  /**
   * @brief Fmd
  */
  FILE_MANAGER<double> Fmd;
  //////////////////////////////////////////////////////////
  /**
   * @brief Random generator
  */
  gsl_rng * r;
  //////////////////////////////////////////////////////////
  /**
  * @brief  Type of code running. Do not touch it.
  */
  string code;
  //////////////////////////////////////////////////////////
  /**
   * @brief Statistics to measure, basically Cl
  */
  string statistics;
  //////////////////////////////////////////////////////////
  /**
   * @brief Name of the input catalog from wich the
   * Statistics will be measured. Overwritten by param file
  */
  string name_catalog;
  //////////////////////////////////////////////////////////
  /**
   * @brief Output directory where the Cl will be written
   * Overwritten by param file
  */
  string name_output_dir;
  //////////////////////////////////////////////////////////
  /**
   * @brief Input directory where the input cat is located
   * Overwritten by param file
  */
  string name_input_dir;
  //////////////////////////////////////////////////////////
  /**
   * @brief Name of catalog
   * Overwritten by param file
  */
  string input_file;
  //////////////////////////////////////////////////////////
  /**
   * @brief Type of file for the galaxy catalog, options are
   * ascii and fits
   * Overwritten by param file
  */
  string file_type;

  //////////////////////////////////////////////////////////
  /**
   * @brief Type of file for the Healpix mask, options are
   * ascii and fits
   * Overwritten by param file
  */
  string file_type_mask;

  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  int i_alpha;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  int i_delta;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  int i_z;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  int i_M;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  int i_w;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  int i_alpha_ran;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  int i_delta_ran;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  int i_z_ran;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  int i_M_ran;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  int i_w_ran;
  //////////////////////////////////////////////////////////
  /**


  //////////////////////////////////////////////////////////
  /**
   * @brief  We can use a random catalog.
   * In that case the map built from this catalog will
   * play the role of mask.
   * Overwritten by param file
  */
  bool use_random_cat;

  //////////////////////////////////////////////////////////
  /**
   * @brief This parameters
   *  will help us to avoid the calculation of Jlm or Ilm
   *  Options are full_sky /  masked_sky
   *  (def masked_sky)
  */
  string sky;

  //////////////////////////////////////////////////////////
/**
 * @brief Choose the coordinate system between galactic and
 * equatorial. This is only useful in case we split between
 * south and  north. When hemis=all, by default, we use galactic
 *  coordinates according to te input mask.
 *  The code internaly selects north / south if galacitic
 *  but demands an input mask if we split in equatorial
*/
  string coord;

  //////////////////////////////////////////////////////////
  /**
   * @brief Input file for the mask in HealPix format
  */
  string input_file_mask;
  //////////////////////////////////////////////////////////
  /**
   * @brief Input file for the north equatorial mask
  */
  string input_file_mask_north_equatorial;
  //////////////////////////////////////////////////////////
  /**
   * @brief Input file for the southern equatorial mask
  */
  string input_file_mask_south_equatorial;
  //////////////////////////////////////////////////////////
  /**
   * @brief Input file for the fulls sky mask
  */
  string input_file_mask_fs;

  //////////////////////////////////////////////////////////
  /**
   * @brief no (yes) if you (do not) want to generate fits files
   * containing the overdensity map. Note that if these are already
   * created, the code will stop. Delete fits files first.
  */
  bool generate_fits_files;

  //////////////////////////////////////////////////////////
  /**
   * @brief Number of REDSHIFT BINS z_bins, used if selection=fls
  */
  int n_z_bins;
  //////////////////////////////////////////////////////////
    /** @brief Select the type of redshift bins. If the variable
     *  define_z_bins is set to "delta",
     *  the code takes zmin and zmax (either given here or
     * taken from the catalog,  as has been specified in the
     * variable min_max_from_cat)  and generate n_z_bins in
     * redshift with the same width. If set  to "number",
     * the code take the same redshift range and  generate
     * n_z_bins each containing the same number of galaxies
*/
  string define_z_bins;
  //////////////////////////////////////////////////////////
  /**
   * @brief Determines hether to use Healpix or direct sum over galaxies
  */
  string sampling;
  //////////////////////////////////////////////////////////
  /**
   * @brief The kind of Peeble -lile estimator to use
  */
  string type_of_P_estimator;
  //////////////////////////////////////////////////////////
  /**
   * @brief Type of selection. This is fixed
  */
  string selection;
  //////////////////////////////////////////////////////////
  /**
   * @brief Hemisphere
  */
  string hemis;
  //////////////////////////////////////////////////////////
  /**
   * @brief Define the type of redshift
  */
  string redshift;
  //////////////////////////////////////////////////////////
  /**
   * @brief Type of Magnitude related property, off
  */
  string property_i_M_type;
  //////////////////////////////////////////////////////////
  /**
   * @brief TYpe of redshift bins
  */
  string select_z_bins;
  //////////////////////////////////////////////////////////
  /**
   * @brief Name of random ascii file
  */
  string input_file_random;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  string output_file_raw ;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  string output_file;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  string output_file_Jlm;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  string shot_noise_correction;
  ///////////////////////////////////////////////////////////
  /**
   * @brief
  */
  string output_file_window;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  string bin_type;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  string fits_map;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  bool compute_jlm;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  int Lmin;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  int Lmax;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  int nside;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  int N_L_bins;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  int ngal;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  bool compute_mixing_matrix;
  //////////////////////////////////////////////////////////

  /**
   * @brief
  */
  bool mixing_matrix_exact=false;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  bool use_weight;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  bool use_Mk_min_max_from_cat;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  bool use_z_min_max_from_cat;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  bool compute_property_weighted_cl;
  //////////////////////////////////////////////////////////
  /**
   * @brief
*/
  bool use_external_pk_file;
  //////////////////////////////////////////////////////////
   /**
    * @brief
  */
  int i_mask_pixel;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  int i_mask_alpha;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  int i_mask_delta;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  int i_mask_flag;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  int n_M_bins;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  double MKmin;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  double MKmax;

  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  double zmin;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  double zmax;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  int nran;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  double alpha;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  double zmin_bin;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  double zmax_bin;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  int number_of_realizations;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  double zmax_cat;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  double zmin_cat;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  double MKmax_cat;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  double MKmin_cat;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  double area_pixel;
  //////////////////////////////////////////////////////////
  /**
   * @brief Expected number of galaxies in one pixel
  */
  double mean_number_galaxies_pix;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  double mean_number_galaxies; //expected number of galaxies. Ngal/Area
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  double rms_ngal;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  double shot_noise;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  double area_survey;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  double sky_fraction;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<double> theta_new;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<double> phi_new;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<double> theta_new_pix;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<double> Mean_ngal_pix; // mean surface density of galaxies per pixel
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<double> Mean_ngal; //mean surface density of galaxies in sample
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<double> Shot_Noise;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<double> pixmask;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<vector<double> > Jlm;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector to the Alms in Magnitude or Redshift bins
   */
  vector<vector<vector<complex<double> > > > Blm;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector to allocate zmin and zmax of nbins in z.
   * The first 0-0 refers to the fill z-interval. Used only for
   * FLS z-interval. Used only for FLS
   */
  vector< vector<double> >zminmax;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector containing the Window function
  */
  vector<double> Wl;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector containing the l-modes
  */
  vector<int> lvec;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<double> Clvec;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<double> lbin;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<double> lbin_min;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<double> lbin_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<double> Clbin;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<double> Clbin_meas;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<double> eClbin;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of angular modes
  */
  vector<int> nmodes;
  //////////////////////////////////////////////////////////
  /**
   * @brief ANgular power of the mask in lbins
  */
  vector<double> Wlbin;
  //////////////////////////////////////////////////////////
  /**
   * @brief Mixing matrix Rll
  */
  vector<vector<double> > R;
  //////////////////////////////////////////////////////////
  /**
   * @brief MIxing matrix is bins, l-lbin
  */
  vector<vector<double> > Rll_bin;
  //////////////////////////////////////////////////////////
  /**
   * @brief NUmber of used galaxies, i.e, not masked
  */
  int ngal_used;
  //////////////////////////////////////////////////////////
  /**
   * @brief NUmber of used random objects, i.e, not masked
  */
  int nran_used;
  //////////////////////////////////////////////////////////
  /**
   * @brief NUmber of observed pixels from the mask
  */
  int n_observed_pixels;
  //////////////////////////////////////////////////////////
  /**
   * @brief NUmber of pixels in the mask
  */
  int n_total_pixels;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  int n_pixels;
  //////////////////////////////////////////////////////////
  /**
   * @brief Get the parameters from input file
  */
  void read_pars(string &);

  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  void set_vectors();
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  void set_ilm_jlm();
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  void set_mean_ngal_pix(int);

  //////////////////////////////////////////////////////////
  /**
   * @brief   Index used to pindown the redshift label
  */
  void set_IZ(int);
  //////////////////////////////////////////////////////////
  /**
   * @brief   Define the z-bins
   */
  void set_zbins();
  //////////////////////////////////////////////////////////
  /**
   * @brief   Define the type of L-bins
  */
  void set_Lbins();
  //////////////////////////////////////////////////////////
  /**
   * @brief  Function relating Healpix index and ring
  */
  int npix_ring(int, int);
  //////////////////////////////////////////////////////////
  /**
   * @brief  COMpute and set Healpix related quantities from the mask
  */
  void set_healpix_pars();
  //////////////////////////////////////////////////////////
  /**
   * @brief  Read the mask from input file
  */
  void get_mask(string);
  //////////////////////////////////////////////////////////
  /**
   * @brief   Get the map from the cat
  */
  void get_map(char,string, Healpix_Map<double>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief  Array for the ALm in redshift bins
  */
  void set_BLMZ(int);
  //////////////////////////////////////////////////////////
  /**
   * @brief  Get parameters associated to the CL estimator
  */
  void get_pars_cl_estimator();
  //////////////////////////////////////////////////////////
  /**
   * @brief   Get parmeters asociated to the mask
  */
  void get_pars_mask();
  //////////////////////////////////////////////////////////
  /**
   * @brief Get the mixing matrix
  */
  void get_mixing_matrix();
  //////////////////////////////////////////////////////////
  /**
   * @brief Write the parameters in screen
  */
  void write_pars_cl_estimator();
  //////////////////////////////////////////////////////////
  /**
   * @brief COmputes the Wigner symbols
  */
  void W3J();
  //////////////////////////////////////////////////////////
  /**
   * @brief   Get the Map from the catalog
  */
  void Cat2Map(char, Healpix_Map<double>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief  Get the map from the catalog,epressed in terms of deltas
  */
  void Cat2Map_delta(char, Healpix_Map<double>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  void Cat2Map_noz(char, Healpix_Map<double>&);

   //////////////////////////////////////////////////////////
  /**
   * @brief  Get the ALms from the map, our way
  */
  void Map2Alm(Healpix_Map<double>, Alm<xcomplex <double> >&,  Alm<xcomplex <double> >&, arr<arr<double> >& );
  //////////////////////////////////////////////////////////
  /**
   * @brief   Get the Alm from the map
  */
  void Map2Alm(Healpix_Map<double>, Alm<xcomplex <double> >&);
  //////////////////////////////////////////////////////////
  /**
   * @brief  Get the map from the Alms , using Healpix
  */
  void Alm2Map(Alm<xcomplex <double> >&, Healpix_Map<double>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief  Get Alm from the random map
  */
  void Map2Alm_ran(Healpix_Map<double>,Healpix_Map<double>, Alm<xcomplex <double> >&,  Alm<xcomplex <double> >&);
  //////////////////////////////////////////////////////////
  /**
   * @brief   Get the function Jlm used n the D-like estimator
  */
  void Map2Jlm(Healpix_Map<double>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  void Map2Ilm_Jlm(Alm<xcomplex <double> >&);
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  void Alm2Cl(Alm<xcomplex <double> >&, vector<int>&, vector<double>&);
  /**
   * @brief   .
  */
  //////////////////////////////////////////////////////////
  /**
   * @brief   COnvert Alm to Cl when partial sky coverage is prensent
  */
  void Alm2Cl(string est, Alm<xcomplex <double> >&, Alm<xcomplex <double> >&, arr<arr<double> >&, vector<int>&, vector<double>&, vector<double>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief  COnvert Alm to Cl for full sky
  */
  void Alm2Cl(int, int, Alm<xcomplex <double> >&,  vector<double> &);  // the one is
  //////////////////////////////////////////////////////////
  /**
   * @brief   Get Gaussian errors on the estimateos of Cl
  */
  void get_eCl(int, int, vector<double>,vector<double>,vector<double>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief  Get mixing matrix
  */
  void Mixing_matrix();
  //////////////////////////////////////////////////////////
  /**
   * @brief ANgular power in bins
  */
  void Cl_bins(vector<double>, vector<double>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief Gaussian errors in bins
  */
  void eCl_bins(vector<double>, vector<double>&);
   //////////////////////////////////////////////////////////
  /**
   * @brief Get the Mixing matix in bins
  */
  void Rll_bins();
  //////////////////////////////////////////////////////////
  /**
   * @brief Main member function doing all the job
  */
  void get_my_power();
  //////////////////////////////////////////////////////////
  /**
  *  @brief  Obtain minimum value of the i-th column of gal properties
  *          The column ic is passed explicitely
  */
  double get_min(char, int ic);
  //////////////////////////////////////////////////////////
   /**
   *  @brief  Obtain maximum value of the i-th column of gal properties
   *          The column is passed explicitely
   */
  double get_max(char, int);



};

#endif

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************




void Cl_FUNCTIONS::get_my_power(){


  string file_Jlm_init_fixed=this->output_file_Jlm+".dat";
  this->output_file_Jlm=file_Jlm_init_fixed;
  
  if(this->file_type=="ascii")this->read_input_cats("g",this->input_file);


  if(!this->use_random_cat){
    // If using a mask instead of a random
    // catalog, use the nside (resolution) from this point on
    string mask_f=(this->sky=="full_sky"? this->input_file_mask_fs: this->input_file_mask);
    this->get_mask(mask_f);
    this->nside=floor(sqrt((double)this->n_pixels/12.));
    cout<<"Mask with nside = "<<this->nside<<endl;
  }

  // *********************************************************

  set_healpix_pars();
  get_pars_mask();
  set_vectors();     // Allocate memory for used vectors


  // *********************************************************


  if(this->n_z_bins!=-1){
    cout<<RED<<"SELECTION INFO "<<RESET<<endl;
    cout<<CYAN;
    if(selection=="fls"){
      cout<<"Flux limited sample selected. Bins in redshift"<<endl;
      if(this->define_z_bins=="delta")cout<<"zBins with constant width"<<endl;
      if(this->define_z_bins=="number")cout<<"zBins chosen with equal number of galaxies"<<endl;
    }
    
    // Get zmax
    this->zmax_cat=this->get_max('r',this->i_z);
    this->zmin_cat=this->get_min('d',this->i_z);
    cout<<RED<<"REDSHIFT INFO "<<RESET<<endl;
    cout<<CYAN<<"Maximim redshift found in catalogue = "<<this->zmax_cat<<endl;
    cout<<CYAN<<"Minimum redshift found in catalogue = "<<this->zmin_cat<<endl;
    if(!this->use_z_min_max_from_cat){
      cout<<CYAN<<"Maximim redshift selected in parameter file = "<<this->zmax<<endl;
      cout<<CYAN<<"Minimum redshift selected in parameter file = "<<this->zmin<<endl;
    }
    else{
      this->zmin=this->zmin_cat;
      this->zmax=this->zmax_cat;
    }
    std::cout<<"*******************************************************"<<endl;
    // Get Magnitude limits
    this->MKmax_cat=this->get_max('d',this->i_M);
    this->MKmin_cat=this->get_min('d',this->i_M);
    cout<<RED<<"MAGNITUDE Mk INFO "<<RESET<<endl;
    cout<<CYAN<<"Maximim Mk found in catalogue = "<<this->MKmax_cat<<endl;
    cout<<CYAN<<"Minimum Mk found in catalogue = "<<this->MKmin_cat<<endl;
    if(!this->use_Mk_min_max_from_cat){
      cout<<CYAN<<"Maximim Mk selected in parameter file = "<<this->MKmax<<endl;
      cout<<CYAN<<"Minimum Mk selected in parameter file = "<<this->MKmin<<endl;
    }
    else{
      this->MKmin=this->MKmin_cat;
      this->MKmax=this->MKmax_cat;
    }
    std::cout<<"*******************************************************"<<endl;
    
    // SOME WARNINGS:
    std::cout<<RED;
    if(this->MKmax > this->MKmax_cat)cout<<"Warning: Maximum Mk in parameter file ("<<this->MKmax<<") greater than maximum value found in the catalog ("<<this->MKmax_cat <<")"<<endl;
    if(this->MKmin < this->MKmin_cat)cout<<"Warning: Minimum Mk in parameter file ("<<this->MKmin<<") smaller than minimum value found in the catalog ("<<this->MKmin_cat <<")"<<endl;
    if(this->zmax > this->zmax_cat){
      //cout<<"Warning: Maximum z in parameter file ("<<this->zmax<<") greater than maximum value found in the catalog ("<<this->zmax_cat <<"). Setting value from catalog"<<endl;
      cout<<"Warning: Maximum z in parameter file ("<<this->zmax<<") greater than maximum value found in the catalog ("<<this->zmax_cat <<"). Doing nothing"<<endl;
      //   this->zmax= this->zmax_cat;
    }
    if(this->zmin < this->zmin_cat){
      cout<<"Warning: Minimum z in parameter file ("<<this->zmin<<") smaller than minimum value found in the catalog ("<<this->zmin_cat <<"). Doing nothing"<<endl;
    }
  }
  
  else{
    this->zmin=-1000;
    this->zmax=+1000;
  }
  
  
  // *******************************************************
  
  int NB;
  if(this->n_z_bins!=-1)NB = (this->n_z_bins);
  else NB=0;
  
  this->set_BLMZ(NB);
  
  // *******************************************************
  // Feed class with redshift bins
  
  if(this->n_z_bins!=-1)this->set_zbins();
  
  if(this->use_random_cat)this->read_input_cats("r", this->input_file_random);

  // ******************************************************************************************************************
  // *********************************************************
  // *********************************************************
  // *********************************************************
  // OPERATIONS IN HARMONIC SPACE
  // *********************************************************
  // *********************************************************
  // ******************************************************************************************************************

  
  // Estos los dejo definidos aca en directo...
  Alm<xcomplex <double> > Ilm(this->Lmax+1,this->Lmax+1);
  // Compute Jlm, Ilm. These objects are inizialized to the
  // full sky- case, Jlm=1, Ilm=0 (l>0)
  // such that if not computed in Map2ILM_Jlm, they remain full sky.
  // but if used, these have to be initialized to zero in that funcition
  // For full sky, Ilm os not set to zero by directly not used in the estimation

  this->set_ilm_jlm(); // So far here we only set Jlm
  if(this->use_random_cat)cout<<"Nothing to do here, read cat"<<endl;
  else if(!this->use_random_cat){  //IF WE HAVE A MASK instead of a random, get the Ilm, even IF DIRECT SUM ESTIMATION
    if(this->sky!="full_sky")this->Map2Ilm_Jlm(Ilm);
  }
  if(this->sky=="full_sky"){ // Intead if we have full sky, set Ilm to 1
    for(int i=0;i<this->Lmax+1;i++)for(int im=0;im<this->Lmax+1;im++)Ilm(i,im).real(1);//Just for full sky
    for(int i=0;i<this->Lmax+1;i++)for(int im=0;im<this->Lmax+1;im++)Ilm(i,im).imag(1);
  }

  cout<<"*******************************************************"<<endl;
  cout<<BLUE<<"Computing Alm in z bins"<<RESET<<endl;

  // *******************************************************
  int izf;
  int imf;
  int icut=0;

  if(this->n_z_bins!=-1){
    imf=0;
    izf=(this->n_z_bins==1? 0 : this->n_z_bins);
    this->set_mean_ngal_pix(this->n_z_bins);
  }
  else{
    this->set_mean_ngal_pix(1);
    izf=0;
  }
  // *******************************************************
  // *******************************************************

  this->set_Lbins();
  
  // START LOOP OVER THE REDSHIFT BINS ONLY IF FLS
  for(int IZ=0;IZ<=izf;++IZ){
      
    this->set_IZ(IZ);
      
    ///Get map from cat
    Healpix_Map<double>map(log2(this->nside), RING);
    this->get_map('d', this->file_type, map);
      
    // Define weights
    arr<double>weight_data(2*map.Nside());
    weight_data.fill(1.0);
      
    Healpix_Map<double>map_ran(log2(this->nside), RING);
    if(this->use_random_cat){
      this->get_map('r',this-> file_type, map_ran);
    }
      
    // ***************************************************************************************
    this->get_pars_cl_estimator();
    this->write_pars_cl_estimator();
    // ***************************************************************************************
    this->Mean_ngal_pix[IZ]=this->mean_number_galaxies_pix;
    this->Mean_ngal[IZ]=this->ngal_used/(this->area_survey);
    this->Shot_Noise[IZ]=this->shot_noise;
      
    // ***************************************************************************************
    Alm<xcomplex <double> > Alm(this->Lmax+1,this->Lmax+1);
    for(int i=0;i<this->Lmax+1;i++)for(int im=0;im<this->Lmax+1;im++)Alm(i,im).real(0);
    for(int i=0;i<this->Lmax+1;i++)for(int im=0;im<this->Lmax+1;im++)Alm(i,im).imag(0);
    weight_data.fill(1);
      
    cout<<RED<<"Computing Alm..."<<RESET<<endl;
    map2alm(map,Alm, weight_data,false); //HealPix
    cout<<"Done"<<endl;
      
    if(this->use_random_cat){
      arr<double>weight_random(4*map.Nside()); //Do not understand why 2*map.Nside()
      for(int i=0; i < 4*map.Nside(); ++i) weight_random[i]=1.0;
      map2alm(map_ran,Ilm, weight_random,false);
      this->Map2Jlm(map_ran);
    }
    cout<<"*******************************************************"<<endl;
      
    // ***********************************************************************************************
    // Fill the vectors Blm. These are the Alm in different redshift or Magnitude bins
    for(int il=0;il<this->Lmax+1;il++)for(int im=0;im<this->Lmax+1;im++)this->Blm[il][im][IZ].real(Alm(il,im).real());
    for(int il=0;il<this->Lmax+1;il++)for(int im=0;im<this->Lmax+1;im++)this->Blm[il][im][IZ].imag(Alm(il,im).imag());
      
  } //close loop over z, opened in sample=fls
    // ***********************************************************************************************
    
  cout<<BLUE;
  // Having allocated quantities in z bin (fls)
  // we now do a double loop (fls) over the zbins, or
  // we do only one loop (vls) (the second is currently running)
  // and compute things for the current Mbins
    
  cout<<"*******************************************************"<<endl;
  
  if(this->n_z_bins>1){
    cout<<"Computing Cross Cls"<<RESET<<endl;
    for(int IZ=1;IZ<=izf;IZ++){
      vector<double> ClzbinI(this->Lmax+1,0);
      this->Alm2Cl(IZ,IZ, Ilm, ClzbinI);//Auto Cl bin IZ
      for(int JZ=IZ;JZ<=izf;JZ++){
	vector<double> Cl(this->Lmax+1,0);
	vector<double> eCl(this->Lmax+1,0);
	vector<double> ClzbinJ(this->Lmax+1,0);
	
	// ***********************************************************************************************
	// We need to compute the autos within each loop in order to compute th
	// variance for the cross power spectrum (C_i+sn_j)*(C_j+sn_j).
	// NOTE THAT FIRST WE NEED TO COMPUTE THE AUTO. DO not change ordering here.
	// ***********************************************************************************************
	//Auto Cl in bin J:
	this->Alm2Cl(JZ,JZ,Ilm,ClzbinJ);
	//Cross:
	this->Alm2Cl(IZ,JZ,Ilm,Cl);
	//get variance of cross power:
	this->get_eCl(IZ,JZ,ClzbinI, ClzbinJ, eCl);
	
	// ***********************************************************************************************
	// Write raw estimates
	string rest;
	if(this->sampling=="H"){
	  rest="_"+this->sampling+"_Nside_"+to_string(this->nside)+"_zbins_"+to_string(IZ)+"_"+to_string(JZ)+"_"+this->hemis+"_"+this->coord+".dat";
	}
	else if(this->sampling=="DS"){
	  rest="_"+this->sampling+"_zbins_"+to_string(IZ)+"_"+to_string(JZ)+"_"+this->hemis+"_"+this->coord+".dat";
	  
	}
	string cfile=this->output_file_raw+rest;
	this->Fmi.write_to_file(cfile, this->lvec, Cl, this->Wl,eCl);
	// ***********************************************************************************************
	// Cl in Bins
	this->Cl_bins(Cl, this->Clbin);
	// ***********************************************************************************************
	
	if(IZ==JZ){
	  int ll=-1;
	  do{ll++;}while(this->Clbin[ll]>this->Shot_Noise[JZ]);
	  if(ll==this->lbin.size())std::cout<<RED<<"No shot noise domination reached for estimates in redshift bin "<<JZ<<". Select lmax = 100. "<<RESET<<std::endl;
	  else std::cout<<RED<<"Maxium scale before shot-noise regime for Auto Cl in zBin "<< JZ<<" is lmax = "<<(int)this->lbin[ll]<<RESET<<std::endl;
	}
	// ***********************************************************************************************
	// Variance in Bins
	this->eCl_bins(eCl,this->eClbin);
	// Wl in Bins
	this->Cl_bins(this->Wl, this->Wlbin);
	// ***********************************************************************************************
	// Write Binned estimates
	string ofile=this->output_file+rest;
	this->Fmd.write_to_file(ofile, this->lbin, this->lbin_min, this->lbin_max,this->Clbin,this->eClbin);
      }
    }
    
    {  // Get the Cl of the full sample
      vector<double> Cl(this->Lmax+1,0);
      vector<double> eCl(this->Lmax+1,0);
      this->Alm2Cl(0,0, Ilm,Cl);
      this->get_eCl(0,0,Cl,Cl,eCl);
      // ***********************************************************************************************
      string rest;
      rest="_"+this->sampling+"_Nside_"+to_string(this->nside)+"_zbins_"+to_string(0)+"_"+to_string(0)+"_"+this->hemis+"_"+this->coord+".dat";
      	
      string cfile=this->output_file_raw+rest;
      this->Fmi.write_to_file(cfile, this->lvec, Cl, this->Wl, eCl);
      // ***********************************************************************************************
      this->Cl_bins(Cl,  this->Clbin);
      this->eCl_bins(eCl, this->eClbin);
      // ***********************************************************************************************
      this->Cl_bins(this->Wl, this->Wlbin);
	
	
      {
	int ll=-1;
	do{ll++;}while(this->Clbin[ll]>this->Shot_Noise[0]);
	if(ll==this->lbin.size())std::cout<<RED<<"No shot noise domination reached for estimates in redshift bin "<<0<<". Select lmax = 100. "<<RESET<<std::endl;
	else std::cout<<RED<<"Maxium scale before shot-noise regime for Auto Cl in zBin "<< 0<<" is lmax = "<<(int)this->lbin[ll]<<RESET<<std::endl;
      }
      // ***********************************************************************************************
      // ***********************************************************************************************
      string ofile=this->output_file+rest;
      this->Fmd.write_to_file(ofile, this->lbin, this->lbin_min, this->lbin_max,this->Clbin,this->eClbin);
    }
      
  }

  else{   // Get the Cl when no redshift info is available
    vector<double> Cl(this->Lmax+1,0);
    vector<double> eCl(this->Lmax+1,0);
    this->Alm2Cl(0,0, Ilm,Cl);
    this->get_eCl(0,0,Cl,Cl,eCl);
    // ***********************************************************************************************
    string rest;
    rest="_"+this->sampling+"_Nside_"+to_string(this->nside)+"_"+this->hemis+"_"+this->coord+".dat";
    
    string cfile=this->output_file_raw+rest;
    this->Fmi.write_to_file(cfile, this->lvec, Cl, this->Wl, eCl);
    // ***********************************************************************************************
    this->Cl_bins(Cl,  this->Clbin);
    this->eCl_bins(eCl, this->eClbin);
    // ***********************************************************************************************
    this->Cl_bins(this->Wl, this->Wlbin);
    {
      int ll=-1;
      do{ll++;}while(this->Clbin[ll]>this->Shot_Noise[0]);
      if(ll==this->lbin.size())std::cout<<RED<<"No shot noise domination reached for estimates in redshift bin "<<0<<". Select lmax = 100. "<<RESET<<std::endl;
      else std::cout<<RED<<"Maxium scale before shot-noise regime for Auto Cl in zBin "<< 0<<" is lmax = "<<(int)this->lbin[ll]<<RESET<<std::endl;
    }
    string ofile=this->output_file+rest;
    this->Fmd.write_to_file(ofile, this->lbin, this->lbin_min, this->lbin_max,this->Clbin,this->eClbin);
  }
  

  // Finally, get the mixig matrix, if desired
  if(this->compute_mixing_matrix)this->get_mixing_matrix();
    
  this->  prop.clear();
    
  return;
    
}

// *******************************************************************************************************************************************************
// *****************************************************************************************************************
// *******************************************************************************************************************************************************
// *****************************************************************************************************************



void Cl_FUNCTIONS::read_pars(string &file){

 // ***********************************************************************************************

  // Intializing parameters:
  nzbins=0;
  code="cross";
  statistics = "Cl";
  name_catalog = "Galaxy_Catalog";
  name_output_dir ="../Output/";
  name_input_dir = "../Input/";
  input_file = "catalog.dat";
  file_type="ascii";
  file_type_mask="ascii";
  use_random_cat = false;
  sky="masked";
  coord="galactic";
  generate_fits_files=false;
  n_z_bins=-1;
  define_z_bins= "delta";
  bin_type="linear";
  sampling="H";
// ***********************************************************************************************

    // Reading from input file
  
  ifstream fin_parameters (file.c_str());
  if (!fin_parameters) { cerr <<"Error in opening the parameters file "<<file<<"!"<<endl; exit(1); }

  string line_in_file;
  string par_name;
  string equality;
  string par_value;

  
  while (getline(fin_parameters,line_in_file)) {
    // ignore lines starting with hashtag
    if(line_in_file[0] != '#' && line_in_file.empty()==0) {
      stringstream line_string (line_in_file);      // read parameter name
      line_string << line_in_file;
      line_string >> par_name;
      line_string >> equality;      // check that second word is "="
      if (equality != "=") {
	cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
	cerr << "Using a default value for " << par_name << endl; exit(1);
      }
      line_string >> par_value;      // read parameter value
      if (par_value.empty()) {
	cout << "Value of " << par_name << " not specified in " <<file << endl;
	cout << "Assuming a default value for " << par_name << endl;
	continue;
      }
      if (par_name == "Lmax")Lmax = atoi(par_value.c_str());
      else if (par_name == "Lmin")Lmin = atoi(par_value.c_str());
      else if (par_name == "N_L_bins")N_L_bins = atoi(par_value.c_str());
      else if (par_name == "nside")nside = atoi(par_value.c_str());

      else if (par_name == "file_type")file_type = par_value;
      else if (par_name == "file_type_mask")file_type_mask = par_value;
      else if (par_name == "name_input_dir")name_input_dir = par_value;
      else if (par_name == "input_file")input_file = par_value;
      else if (par_name == "bin_type")bin_type = par_value;
      else if (par_name == "name_catalog")name_catalog = par_value;
      else if (par_name == "name_output_dir") name_output_dir = par_value;

      else if (par_name == "use_random_cat"){
        if(par_value=="yes")use_random_cat=true;
	else use_random_cat=false;
      }
      else if (par_name == "generate_fits_files"){
        if(par_value=="yes")generate_fits_files=true;
        else generate_fits_files=false;
      }

      else if (par_name == "compute_mixing_matrix"){
        if(par_value=="yes")compute_mixing_matrix=true;
	else compute_mixing_matrix=false;
      }
      else if (par_name == "mixing_matrix_exact"){
        if(par_value=="yes") mixing_matrix_exact=true;
	else mixing_matrix_exact=false;
      }
      else if (par_name == "use_z_min_max_from_cat"){
        if(par_value=="yes") use_z_min_max_from_cat=true;
	else use_z_min_max_from_cat=false;
      }
      else if (par_name == "use_weight"){
        if(par_value=="yes") use_weight=true;
        else use_weight=false;
      }
      else if (par_name == "compute_jlm"){
        if(par_value=="yes")compute_jlm=true;
        else compute_jlm=false;
      }
      else if (par_name == "compute_property_weighted_cl"){
        if(par_value=="yes")this->compute_property_weighted_cl=true;
        else this->compute_property_weighted_cl=false;
      }
      else if (par_name == "define_z_bins")this->define_z_bins= par_value;
      else if (par_name == "sky")this->sky= par_value;
      else if (par_name == "type_of_P_estimator")this->type_of_P_estimator= par_value;
      else if (par_name == "hemis")this->hemis = par_value;
      else if (par_name == "redshift")this->redshift = par_value;
      else if (par_name == "shot_noise_correction")this->shot_noise_correction = par_value;
      else if (par_name == "input_file_random")this->input_file_random = par_value;
      else if (par_name == "input_file_mask") this->input_file_mask = par_value;
      else if (par_name == "input_file_mask_fs") this->input_file_mask_fs = par_value;
      else if (par_name == "input_file_mask_north_equatorial") this->input_file_mask_north_equatorial = par_value;
      else if (par_name == "input_file_mask_south_equatorial") this->input_file_mask_south_equatorial = par_value;
      else if (par_name == "coord")this->coord = par_value.c_str();
      else if (par_name == "output_file_raw") this->output_file_raw = par_value;
      else if (par_name == "output_file")this->output_file = par_value;
      else if (par_name == "output_file_Jlm") this->output_file_Jlm = par_value;
      else if (par_name == "output_file_window") this->output_file_window = par_value;
      else if (par_name == "fits_map")this->fits_map = par_value;
      else if (par_name == "i_alpha") this->i_alpha= atoi(par_value.c_str());
      else if (par_name == "i_delta") this->i_delta= atoi(par_value.c_str());
      else if (par_name == "i_z") this->i_z= atoi(par_value.c_str());
      else if (par_name == "i_M") this->i_M= atoi(par_value.c_str());
      else if (par_name == "i_alpha") i_alpha_ran= atoi(par_value.c_str());
      else if (par_name == "i_w") i_w= atoi(par_value.c_str());
      else if (par_name == "i_w_ran") i_w_ran= atoi(par_value.c_str());
      else if (par_name == "i_delta_ran") this->i_delta_ran= atoi(par_value.c_str());
      else if (par_name == "i_z_ran") this->i_z_ran= atoi(par_value.c_str());
      else if (par_name == "i_M_ran") this->i_M_ran= atoi(par_value.c_str());
      else if (par_name == "i_mask_alpha") this->i_mask_alpha= atoi(par_value.c_str());
      else if (par_name == "i_mask_delta") this->i_mask_delta= atoi(par_value.c_str());
      else if (par_name == "i_mask_pixel") this->i_mask_pixel= atoi(par_value.c_str());
      else if (par_name == "i_mask_flag") this->i_mask_flag= atoi(par_value.c_str());
      else if (par_name == "zmin")this->zmin= (double)(atof(par_value.c_str()));
      else if (par_name == "zmax")this->zmax= (double)(atof(par_value.c_str()));
      else if (par_name == "zmin_bin")this->zmin_bin= (double)(atof(par_value.c_str()));
      else if (par_name == "zmax_bin")this->zmax_bin= (double)(atof(par_value.c_str()));



    }
  }

  // Ammend name of files:
  input_file=name_input_dir+input_file;
  output_file=name_output_dir+statistics+"_"+type_of_P_estimator+"_"+name_catalog+"_zbtype_"+this->define_z_bins+"_"+output_file;
  output_file_raw=name_output_dir+statistics+"_"+type_of_P_estimator+"_"+name_catalog+"_zbtype_"+this->define_z_bins+"_"+output_file_raw;
  output_file_window=name_output_dir+statistics+"_"+name_catalog+"_"+output_file_window+"_"+this->hemis+"_"+this->coord;
;
  fits_map=name_output_dir+fits_map;
  if(this->compute_property_weighted_cl)output_file=output_file+"_marked";
  if(this->compute_property_weighted_cl)output_file_raw=output_file_raw+"_marked";
  output_file_Jlm=name_output_dir+output_file_Jlm+"_2MPZ_mask_nside_"+to_string(this->nside)+"_Lmax_"+to_string(this->Lmax)+"_"+this->hemis+"_"+this->coord;


  cout<<RESET<<endl;
 }

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *****************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void Cl_FUNCTIONS::set_BLMZ(int nb){
  // Resize Blm. I have to do it in an independend method
  this->Blm.resize(this->Lmax+1);
  for(int i=0; i<Blm.size();++i)this->Blm[i].resize(this->Lmax+1);
  for(int i=0; i<this->Blm.size();++i)
    for(int j=0;j<this->Blm[i].size();++j)
      this->Blm[i][j].resize(nb+1);

}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void Cl_FUNCTIONS::set_vectors(){


    this->z_min.resize(this->n_z_bins+1,0);
    this->z_max.resize(this->n_z_bins+1,0);
    this->zminmax.resize(this->n_z_bins+1);
    for(int i=0;i<=this->n_z_bins;i++)this->zminmax[i].resize(2,0);
    this->Wl.resize(this->Lmax+1);
    this->lvec.resize(this->Lmax+1);
    this->lbin.resize(this->N_L_bins,0);
    this->lbin_min.resize(this->N_L_bins,0);
    this->lbin_max.resize(this->N_L_bins,0);
    this->Clbin.resize(this->N_L_bins,0);
    this->eClbin.resize(this->N_L_bins,0);
    this->nmodes.resize(this->N_L_bins,0);
    this->Wlbin.resize(this->N_L_bins,0);


    this->R.resize(this->Lmax+1);
    for(int i=0;i<this->R.size();i++)this->R[i].resize(this->Lmax+1,0);

    // THIS IS ALSO USED WHTN WE READ THE MIXING MATRIX!!!

    if(this->compute_mixing_matrix){
      R.resize(this->Lmax+1);
      for(int i=0;i<this->R.size();i++)R[i].resize(this->Lmax+1,0);
      Rll_bin.resize(this->N_L_bins);
      for(int i=0;i<this->Rll_bin.size();i++)Rll_bin[i].resize(this->Lmax+1,0);
    }
}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************


void Cl_FUNCTIONS::set_Lbins(){

   if(this->bin_type=="linear"){
    double deltal=((double)(this->Lmax-this->Lmin))/((double)lbin.size()); 
    //here and just here I define N_L_bins as a double instead of an integer
    cout<<CYAN<<endl;
    std::cout<<"*********************************************************"<<endl;
    cout<<BLUE<<"L-bins INFO "<<endl;
    cout<<"Bin-averaged power spectrum"<<endl;
    cout<<"Number of bins = "<<this->lbin.size()<<endl;
    cout<<"l-Bin width = "<<deltal<<RESET<<endl;
    std::cout<<"*********************************************************"<<endl;
    cout<<RESET<<endl;
    for(int i=0;i<lbin.size();i++)this->lbin[i]=this->Lmin+(i+0.5)*deltal;
    for(int i=0;i<lbin.size();i++)this->lbin_min[i]=this->Lmin+(i)*deltal;
    for(int i=0;i<lbin.size();i++)this->lbin_max[i]=this->Lmin+(i+1)*deltal;
  }
  else{
    if(this->bin_type=="log"){
      double deltal=log10(Lmax/1.0)/((double)lbin.size());
      for(int i=0;i<lbin.size();i++)lbin[i]=pow(10,log10(1)+(i+0.5)*deltal);
      for(int i=0;i<lbin.size();i++)lbin_min[i]=pow(10,log10(1)+i*deltal);
      for(int i=0;i<lbin.size();i++)lbin_max[i]=pow(10,log10(1)+(i+1)*deltal);
    }
  }
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void Cl_FUNCTIONS::set_IZ(int IZZ){IZ=IZZ;}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void Cl_FUNCTIONS::set_zbins(){

  this->z_min[0]=this->zmin;
  this->z_max[0]=this->zmax;
  double Dz=(this->zmax-this->zmin)/(double(this->n_z_bins));


   if(this->define_z_bins=="number")get_zbins_same_ngal(this->n_z_bins,this->i_z,this->zmin,this->zmax,this->zminmax);

  for(int i=1;i<=this->n_z_bins;i++){
    z_min[i]= (this->define_z_bins=="number"? this->zminmax[i][0]:this->zmin+(i-1)*Dz);
    z_max[i]= (this->define_z_bins=="number"? this->zminmax[i][1]:this->zmin+(i)*Dz);
  }
  cout<<RED<<"REDSHIFT BINS: "<<endl;
  if(this->n_z_bins>1)for(int i=0;i<=this->n_z_bins;i++)cout<<i<<"  ("<<z_min[i]<<"-"<<z_max[i]<<")"<<endl;
  else cout<<"0  ("<<z_min[0]<<"-"<<z_max[0]<<")"<<endl;
  cout<<RESET<<endl;
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void Cl_FUNCTIONS::set_healpix_pars(){
  n_pixels=12*this->nside*this->nside;
  Nrings=4*this->nside-1;
}


// *******************************************************************************************************************************************************
int Cl_FUNCTIONS::npix_ring(int nside_h, int ir){
  // Number of pixels contained in a ring labled with ir. with 0<ir< 4nside-1 rings
  int Npixels_ring;
  if(ir<nside_h-1)Npixels_ring=4*(ir+1);
  if(ir>=nside_h-1 && ir<3*nside_h)Npixels_ring=4*nside_h;
  if(ir>=3*nside_h)Npixels_ring=4*(4*nside_h-ir-1);
  return Npixels_ring;
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void Cl_FUNCTIONS::get_eCl(int iz, int jz, vector<double>cli,vector<double>clj, vector<double>&ecl){
  /// WARNING: THE CROSS POWER SPECTRUM NEEDS HAS UNDESTIMATED GAUSSIAN ERROR BARS, CHECN WHITE, SUNG, PERCIVAL

  for(int l=this->Lmin;l<=this->Lmax;++l){
    ecl[l]=sqrt(  (2./(this->sky_fraction*(2*l+1)))*(     (cli[l]+this->Shot_Noise[iz])*(clj[l]+this->Shot_Noise[jz])) );
  }

}


// *******************************************************************************************************************************************************
void Cl_FUNCTIONS::get_mask(string input_file_mask){


  if(this->file_type_mask=="ascii"){
    FILE_MANAGER<double>Fmd;

    vector<double>mask2;
    this->n_columns_mask=Fmd.read_file2(input_file_mask, mask2);

    // set the n_pixels of the mask
    this->n_pixels=mask2.size()/this->n_columns_mask;
    this->n_total_pixels=this->n_pixels;
    
  // set the nside of the mask
    this->nside=sqrt(this->n_pixels/12);

  // nr=nrings, have to compute it here for nside 
  // is going to be the size of the mask and has not been yet computed
    int nr=4*this->nside-1;
  
    this->theta_new.resize(nr);
    fill(this->theta_new.begin(), this->theta_new.end(), 0);

    this->phi_new.resize(mask2.size()/this->n_columns_mask);
    fill(phi_new.begin(), this->phi_new.end(), 0);


    this->theta_new_pix.resize(n_pixels);
    fill(theta_new_pix.begin(), this->theta_new_pix.end(), 0);


    this->pixmask.resize(mask2.size()/this->n_columns_mask);
    fill(this->pixmask.begin(), this->pixmask.end(), 0);
  
    n_observed_pixels=0;
  
  // Create an auxiliary Healpix map to
  // copy the mask and find the angles
    Healpix_Map<double>mask_aux(log2(this->nside), RING);

    
    for(int i=0;i<n_pixels;i++)mask_aux[i]=mask2[this->i_mask_flag+this->n_columns_mask*i];
    
    pointing point_rev;
    int ip=0;

    for(int ir=0;ir<nr;ir++){
      int Npixels_ring=npix_ring(nside, ir);
      for(int ipr=0;ipr<Npixels_ring;++ipr){
        point_rev=mask_aux.pix2ang(ip);
        this->theta_new[ir]=point_rev.theta;
        this->phi_new[ip]=point_rev.phi;
        this->theta_new_pix[ip]=point_rev.theta;
        int paux;

	int pmask=mask2[this->i_mask_flag+this->n_columns_mask*ip];
	if(this->hemis=="all"){
          paux=(pmask== 0 ? 0 : 1);  //what comes from a real mask
        }
        else{
          if(this->hemis=="north"){
            if(theta_new[ir]>=0.5*M_PI)paux=0; // mask south
	    //            else paux=mask[ip][this->i_mask_flag];
	    else paux=pmask;
          }
          else if(this->hemis=="south"){
            if(theta_new[ir]<0.5*M_PI)paux=0;  //mask north
	    else paux=pmask;
	    //            else paux = mask[ip][this->i_mask_flag];
          }
        }
        this->pixmask[ip]=paux;
        if(paux==1)++n_observed_pixels;
        ++ip;
      }
    }

    mask2.clear();
  }

  else if(this->file_type_mask=="fits"){
      /*
    cout<<CYAN<<"Reading mask from "<<endl;
    cout<<input_file_mask<<endl;
    cout<<RESET<<endl;
    Healpix_Map<double>mask_aux(log2(this->nside), RING);
    read_Healpix_map_from_fits(this->input_file_mask,mask_aux, 1);
    this->nside=mask_aux.Nside();
    gal.mask.resize(12*nside*nside);
    for(int i=0;i<mask.size();++i)mask[i].resize(3,0);
    for(int i=0;i<n_pixels;i++)gal.mask[i][this->i_mask_flag]=mask_aux[i];
    set_healpix_pars();
*/
  }

}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void Cl_FUNCTIONS::get_map(char c, string file_type, Healpix_Map<double>&map){

  if(file_type=="ascii"){
    if(c=='d')cout<<CYAN<<"Building map for real catalogue..."<<endl;
    if(c=='r')cout<<"Building map for random catalogue..."<<RESET<<endl;

    if(this->n_z_bins!=-1)Cat2Map(c,map);
    else Cat2Map_noz(c,map);

    string map_fits;
    map_fits=fits_map+"_zbin_"+to_string(this->IZ)+".fits";
    if(this->generate_fits_files){
      cout<<CYAN<<"Writing map in fits format in file "<<map_fits<<RESET<<endl;
      write_Healpix_map_to_fits(map_fits, map, (PDT)1);
    }
    if(c=='d')ngal=this->  prop.size()/this->n_columns_gal;
    if(c=='r')nran=this->  prop_r.size()/this->n_columns_ran;
    cout<<CYAN<<"Done.   "<<ngal<<"  "<<this->n_columns_gal<<RESET<<endl;
  }
  else{
    if(file_type=="fits"){
      cout<<CYAN<<"Reading map from "<<endl;
      cout<<input_file<<endl;
      cout<<RESET<<endl;
  //    read_Healpix_map_from_fits(input_file,map, 1);
      this->nside=map.Nside();
      set_healpix_pars();
    }
  }

  if(c=='d'){
    rms_ngal=map.rms();
    this->mean_number_galaxies_pix=(double)ngal_used/((double)n_observed_pixels);
  }
  if(c=='r'){
    rms_nran=map.rms();
    this->mean_number_randoms_pix=(double)nran_used/((double)n_observed_pixels);
  }

  
  if(file_type=="ascii"){
    Healpix_Map<double>map_delta(log2(this->nside), RING);
    if(this->n_z_bins!=-1)Cat2Map_delta(c,map_delta);
    
    /*
      string map_fits_delta=fits_map+"_delta_zbin_"+to_string(this->IZ)+".fits";
      cout<<CYAN<<"Writing delta_map in fits format in file "<<map_fits_delta<<RESET<<endl;
      write_Healpix_map_to_fits(map_fits_delta, map_delta, (PDT)1);
      cout<<CYAN<<"Done."<<RESET<<endl;
    */
  
 }
 

}


// *******************************************************************************************************************************************************
void Cl_FUNCTIONS::Cat2Map(char c, Healpix_Map<double>&map){
  
  int ira=(c=='r'? this->i_alpha_ran:this->i_alpha);
  int idec=(c=='r'? this->i_delta_ran:this->i_delta);
  int iz=(c=='r'? this->i_z_ran:this->i_z);
  int iw=(c=='r'? this->i_w_ran:this->i_w);
  int iM=(c=='r'? this->i_M_ran:this->i_M);
  int ncols = (c=='r'? this->n_columns_ran:this->n_columns_gal);
  if(c=='r')nran_used=0;
  if(c=='d')ngal_used=0;
  vector<double> prop_a;
  prop_a = (c=='d'? this->  prop:this->  prop_r) ;
  vector<double>w_property(prop_a.size()/ncols,0);
  
  int nau=0;
  
  // Here we first allocate the normalized property, in case we want it
  if(this->compute_property_weighted_cl){
    double mean_p=0;
    for(int i=0;i<prop_a.size();++i){
      long pixel;
      point.phi=prop_a[ira+ncols*i]*fac;
      point.theta=0.5*M_PI-fac*prop_a[idec+ncols*i];
      double weight1=(this->use_weight? prop_a[iw+ncols*i]:1.0);
      double weight2=(this->compute_property_weighted_cl? w_property[i]:1.0);
      pixel=map.ang2pix(point);
      if(this->pixmask[pixel]==1){
	if(prop_a[iz+ncols*i]<this->z_max[this->IZ] && prop_a[iz+ncols*i]>=this->z_min[this->IZ]){     // Now do the map in the redshift bin:
	  ++nau;
	  if(this->property_i_M_type=="original"){
	    w_property[i]=prop_a[iM+ncols*i];
	    mean_p+=prop_a[iM+ncols*i];
	  }
	  else if (this->property_i_M_type=="log"){
	    w_property[i]=log10(prop_a[iM+ncols*i]);
	    mean_p+=w_property[i];
	  }
	  else if (this->property_i_M_type=="minus_power_ten"){
	    w_property[i]=pow(10,-2.*prop_a[iM+ncols*i]/5.);
	    mean_p+=w_property[i];
	  }
	}
      }
    }
    mean_p/=(double)nau;
    for(int i=0;i<w_property.size();++i)w_property[i]/=mean_p;
    cout<<CYAN<<"Mean value of property "<<mean_p<<RESET<<endl;  
  }
  else {
    for(int i=0;i<w_property.size();++i)w_property[i]=1.0;
  }
  
  
 
  

  vector<double>map_aux(map.Npix(),0);
  nau=0;
  map.fill(0.0);
  for(int i=0;i<prop_a.size()/ncols;i++){
    // These lines
    // allows us to fill the vector map_aux
    // which will tell us the number of empty pixels
    // for all redshifts, used to compute 
    // the total surveyed area.
    long pixel;
    point.phi=prop_a[ira+ncols*i]*fac;
    point.theta=0.5*M_PI-prop_a[idec+ncols*i]*fac;

    double weight1=(this->use_weight? prop_a[iw+ncols*i]:1.0);
    double weight2=(this->compute_property_weighted_cl? w_property[i]:1.0);

    pixel=map.ang2pix(point);
    map_aux[pixel]++;

    if(pixmask[pixel]==1){ 
      //Just to be sure, we go throgh pixels allowed by the mask
      //If the catalogue is already masked with the same mask, fine
      //but there might be cases in which we use a non masked catalogue.
      if(prop_a[iz+ncols*i]< z_max[this->IZ] && prop_a[iz+ncols*i]>=z_min[this->IZ]){     // Now do the map in the redshift bin:
	++nau;
	map[pixel]+=weight1*weight2;
      }
    }
  }
  prop_a.clear();
  if(c=='d')this->ngal_used=nau;
  if(c=='r')this->nran_used=nau;
  

  // Build here theta_new and phi_new
  // from the map of the randoms
  if(this->use_random_cat){
    if(c=='r'){
      // Create an auxiliary Healpix map to
      // copy the mask and find the angles

      int nr=4*nside-1; 
      
      theta_new.resize(nr);
      fill(theta_new.begin(), theta_new.end(), 0);
      
      phi_new.resize(n_pixels);
      fill(phi_new.begin(), phi_new.end(), 0);
      
      Healpix_Map<double>map_aux(log2(nside), RING); 

      for(int i=0;i<n_pixels;i++)map_aux[i]=map[i];
      pointing point_rev;
      int ip=0;
      int no_pix=0;
      for(int ir=0;ir<nr;ir++){
	int Npixels_ring=npix_ring(nside, ir);
	for(int ipr=0;ipr<Npixels_ring;ipr++){
          point_rev=map_aux.pix2ang(ip);
	  theta_new[ir]=point_rev.theta;
	  phi_new[ip]=point_rev.phi;
          if(map_aux[ip]==0)++no_pix;
	  ++ip;
	}
      }
      n_observed_pixels=map.Npix()-no_pix;
      n_total_pixels=map.Npix();
    }
  }
  cout<<"Cat2Map done"<<endl;

}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void Cl_FUNCTIONS::Cat2Map_noz(char c, Healpix_Map<double>&map){
  
  int ira=(c=='r'? this->i_alpha_ran:this->i_alpha);
  int idec=(c=='r'? this->i_delta_ran:this->i_delta);
  int ncols = (c=='r'? this->n_columns_ran:this->n_columns_gal);
  if(c=='r')nran_used=0;
  if(c=='d')ngal_used=0;
  vector<double > prop_a;
  prop_a = (c=='d'? this->  prop:this->  prop_r) ;
  
  int nau=0;
  

  vector<double>map_aux(map.Npix(),0);
  nau=0;
  map.fill(0.0);
  cout<<prop_a.size()/ncols<<endl;
  for(int i=0;i<prop_a.size()/ncols;i++){
    // These lines
    // allows us to fill the vector map_aux
    // which will tell us the number of empty pixels
    // for all redshifts, used to compute 
    // the total surveyed area.
    long pixel;
    point.phi=prop_a[ira+ncols*i]*fac;
    point.theta=0.5*M_PI-prop_a[idec+ncols*i]*fac;

    pixel=map.ang2pix(point);
    map_aux[pixel]++;

    if(pixmask[pixel]==1){ 
      //Just to be sure, we go throgh pixels allowed by the mask
      //If the catalogue is already masked with the same mask, fine
      //but there might be cases in which we use a non masked catalogue.
      ++nau;
      map[pixel]++;
    } 
  }
  prop_a.clear();
  if(c=='d')this->ngal_used=nau;
  if(c=='r')this->nran_used=nau;
  
  
  // Build here theta_new and phi_new
  // from the map of the randoms
  if(this->use_random_cat){
    if(c=='r'){
      // Create an auxiliary Healpix map to
      // copy the mask and find the angles

      int nr=4*nside-1; 
      
      theta_new.resize(nr);
      fill(theta_new.begin(), theta_new.end(), 0);
      
      phi_new.resize(n_pixels);
      fill(phi_new.begin(), phi_new.end(), 0);
      
      Healpix_Map<double>map_aux(log2(nside), RING); 

      for(int i=0;i<n_pixels;i++)map_aux[i]=map[i];
      pointing point_rev;
      int ip=0;
      int no_pix=0;
      for(int ir=0;ir<nr;ir++){
	int Npixels_ring=npix_ring(nside, ir);
	for(int ipr=0;ipr<Npixels_ring;ipr++){
          point_rev=map_aux.pix2ang(ip);
	  theta_new[ir]=point_rev.theta;
	  phi_new[ip]=point_rev.phi;
          if(map_aux[ip]==0)++no_pix;
	  ++ip;
	}
      }
      n_observed_pixels=map.Npix()-no_pix;
      n_total_pixels=map.Npix();
    }
  }
  cout<<"Cat2Map_noz done"<<endl;
  
}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// Compute the overdensity in pixels to write it in a map
void Cl_FUNCTIONS::Cat2Map_delta(char c, Healpix_Map<double>&map){
  int ira=(c=='r'? this->i_alpha_ran:this->i_alpha);
  int idec=(c=='r'? this->i_delta_ran:this->i_delta);
  int iz=(c=='r'? this->i_z_ran:this->i_z);
  int ncols=(c=='r'? this->n_columns_ran:this->n_columns_gal);
  vector<double>prop_a=(c=='r'? this->  prop_r:this->  prop);
  map.fill(0.0);
  for(int i=0;i<prop_a.size();i++){
    long pixel;
    point.phi=prop_a[ira+ncols*i]*fac;
    point.theta=0.5*M_PI-prop_a[idec+ncols*i]*fac;
    pixel=map.ang2pix(point);
    if(pixmask[pixel]==1){
      if(prop_a[iz+ncols*i]< z_max[this->IZ] && prop_a[iz+ncols*i]>=z_min[this->IZ]){     // Now do the map in the redshift bin:
	map[pixel]++;
      }
    }
  }
  prop_a.clear();
  
  string fileo= this->name_output_dir+"2MPZ_delta_pix_dist_"+this->hemis+"_"+this->coord+"_zbin_"+to_string(this->IZ)+".dat";
  ofstream sout; sout.open(fileo.c_str());
  double ngala=this->mean_number_galaxies_pix;
  for(int i=0;i<map.Npix();++i)map[i]=this->mean_number_galaxies_pix==0.0 ? 0.0: ((double)map[i]-ngala)/ngala ;
  for(int i=0;i<map.Npix();++i)if(map[i]!=-1)sout<<log(1+map[i])<<endl;
  int nnpix=0;
  for(int i=0;i<map.Npix();++i)if(map[i]!=-1)nnpix++;
  cout<<CYAN<<"File "<<fileo<<" written with log (1+delta)"<<RESET<<endl;

  double mean_lg=0;
  double var_lg=0;
  for(int i=0;i<map.Npix();++i)if(map[i]!=-1)mean_lg+=log(1+map[i])/((double)nnpix);
  for(int i=0;i<map.Npix();++i)if(map[i]!=-1)var_lg+=pow(log(1+map[i]) - mean_lg,2)/((double)nnpix);
  var_lg=sqrt(var_lg);

  cout<<CYAN<<"Mean log normal = "<<mean_lg<<"  Variance log normal = "<<var_lg<<RESET<<endl;

  cout<<"Cat2Map_delta done"<<endl;
  sout.close();
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void Cl_FUNCTIONS::get_pars_cl_estimator(){
  cout<<RED;
  this->sky_fraction=(double)n_observed_pixels/(double)(n_total_pixels);  
  this->area_pixel=4.*M_PI/(double)n_total_pixels;
  this->area_survey=((double)n_observed_pixels)*area_pixel;
  this->shot_noise=area_survey/((double)ngal_used);
  this->alpha=((double)ngal_used)/((double)nran_used);
}
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void Cl_FUNCTIONS::set_mean_ngal_pix(int nb){
  this->Mean_ngal_pix.resize(nb+1,0);
  this->Mean_ngal.resize(nb+1,0);
  this->Shot_Noise.resize(nb+1,0);
}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************


void Cl_FUNCTIONS::get_pars_mask(){
  cout<<RED;
  this->area_pixel=4.*M_PI/(double)n_total_pixels;
  this->area_survey=(double)n_observed_pixels*area_pixel;
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************


void Cl_FUNCTIONS::write_pars_cl_estimator(){
  cout<<CYAN<<"*******************************************************"<<endl;  
  cout<<"Information computed from the galaxy catalog: "<<endl;
  if(this->n_z_bins!=-1){
    cout<<"Redshift bin "<<IZ<<endl;
    cout<<"z min = "<<this->z_min[IZ]<<endl;
    cout<<"z max = "<<this->z_max[IZ]<<endl;
  }
  cout<<"Number of galaxies in catalogue= "<<this->ngal<<endl;
  cout<<"Number of galaxies used = "<<ngal_used<<endl;
  if(use_random_cat)cout<<"Number of randoms in catalogue = "<<nran<<endl;
  if(use_random_cat)cout<<"Number of randoms used = "<<nran_used<<endl;
  if(use_random_cat)cout<<"alpha ="<<alpha<<endl;
  
  cout<<"Number of pixels in the mask = "<<n_total_pixels<<endl;
  cout<<"Number of pixels in the suveyed area = "<<n_observed_pixels<<endl;
  cout<<"Skyfraction = "<<sky_fraction<<endl;  
  cout<<"Area survey = "<<area_survey<<endl;  
  cout<<"Area pixel = "<<area_pixel<<endl;  
  cout<<"Shot-noise = "<<shot_noise<<endl;
  cout<<"Mean number of galaxies in pixels = "<<mean_number_galaxies_pix<<endl;
  cout<<"Mean number of galaxies= "<<ngal_used/area_survey<<endl;
  cout<<"RMS  = "<<rms_ngal<<endl;
  if(use_random_cat)cout<<"Mean number of randoms  in pixels = "<<mean_number_randoms_pix<<endl;
  if(use_random_cat) cout<<"RMS  = "<<rms_nran<<endl;
  std::cout<<"*******************************************************"<<endl;  
  cout<<RESET;

}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void Cl_FUNCTIONS::Map2Alm(Healpix_Map<double>map, Alm<xcomplex <double> >&alm,  Alm<xcomplex <double> >&Ilm, arr<arr<double> >&Jlm){

  time_t start;  
  time (&start);
  int iring=0;

  for(int i=0;i<Lmax+1;i++)for(int im=0;im<Lmax+1;im++)alm(i,im).real(0);
  for(int i=0;i<Lmax+1;i++)for(int im=0;im<Lmax+1;im++)alm(i,im).imag(0);

  for(int ir=0;ir<Nrings;ir++){
    comp_time(start,Nrings,ir);
    double x=cos(theta_new[ir]);
    int Npixels_ring=npix_ring(nside, ir);
    for(int m=0;m<=Lmax;m++){
      double *Pl= new double[Lmax-m+1];
      gsl_sf_legendre_sphPlm_array(Lmax,m,x,Pl);
      for(int ipr=0;ipr<Npixels_ring;ipr++){
	int ip=ipr+iring;      
      	if(pixmask[ip]==1){	
   	  double cc=cos(double(m)*phi_new[ip]);
  	  double ss=sin(double(m)*phi_new[ip]);
   	  int il=0;
   	  for(int l=m;l<=Lmax;l++){	
  	    double Plm=Pl[il];
   	    double Yreal=Plm*cc;
   	    double Yimag=Plm*ss;
	    /*
	    alm(l,m).re+= area_pixel*Yreal*map[ip];
   	    alm(l,m).im+=-area_pixel*Yimag*map[ip];
	    Ilm(l,m).re+= area_pixel*Yreal;
	    Ilm(l,m).im+=-area_pixel*Yimag;*/
	    Jlm[l][m]  += area_pixel*pow(Plm,2);
	    ++il;
   	  }
   	}
      }
      delete[] Pl;
    }
    iring+=Npixels_ring;
  }
}



// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************


void Cl_FUNCTIONS::Map2Alm(Healpix_Map<double>map, Alm<xcomplex <double> >&alm){

  // This function does the job of the map2alm pf Healpix in longer time

  time_t start;
  time (&start);
  int iring=0;

  for(int i=0;i<Lmax+1;i++)for(int im=0;im<Lmax+1;im++)alm(i,im).real(0);
  for(int i=0;i<Lmax+1;i++)for(int im=0;im<Lmax+1;im++)alm(i,im).imag(0);

  for(int ir=0;ir<Nrings;ir++){
    double x=cos(theta_new[ir]);
    int Npixels_ring=npix_ring(nside,ir);

    double a1=0,a2=0;
    for(int m=0;m<=Lmax;m++){
      double *Pl= new double[Lmax-m+1];
      gsl_sf_legendre_sphPlm_array(Lmax,m,x,Pl);
      for(int ipr=0;ipr<Npixels_ring;ipr++){
	int ip=ipr+iring;      
	if(pixmask[ip]==1){	
	  double cc=cos(double(m)*phi_new[ip]);
  	  double ss=sin(double(m)*phi_new[ip]);
   	  int il=0;
   	  for(int l=m;l<=Lmax;l++){	
  	    double Plm=Pl[il];
   	    double Yreal=Plm*cc;
   	    double Yimag=Plm*ss;
	    a1+=area_pixel*Yreal*map[ip];
	    alm(l,m).real(a1); 
	    a2+=-area_pixel*Yimag*map[ip];
   	    alm(l,m).imag(a2);
   	    ++il;
   	  }
	}
      }
      delete[] Pl;
    }
    iring+=Npixels_ring;
  }



}



// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void Cl_FUNCTIONS::set_ilm_jlm(){
  this->Jlm.resize(this->Lmax+1);
  for(int i=0;i<Jlm.size();++i)this->Jlm[i].resize(this->Lmax+1,1.0);
}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void Cl_FUNCTIONS::Map2Ilm_Jlm(Alm<xcomplex <double> >&Ilm){

  time_t start;  
  time (&start);

  for(int i=0;i<Lmax+1;i++)for(int im=0;im<Lmax+1;im++)Ilm(i,im).real(0);
  for(int i=0;i<Lmax+1;i++)for(int im=0;im<Lmax+1;im++)Ilm(i,im).imag(0);


  int iring=0;
  // ******************************************************************** 
  // In this part we get the Ilm using HEALPIX functions 
  
  cout<<CYAN<<"Computing Ilm...";
  Healpix_Map<double>pmask(log2(nside), RING);
  arr<double>weight_data(2*pmask.Nside(),1.0);
  for(int i=0;i<n_pixels;i++)pmask[i]=this->pixmask[i];
  //map2alm(pmask,Ilm,weight_data,false);
  map2alm_iter(pmask,Ilm,NUM_ITER_MAP2ALM, weight_data);


  cout<<"Done."<<RESET<<endl;
  
  // ********************************************************************

  if(this->type_of_P_estimator=="D"){
      cout<<CYAN<<"Computing Jlm";
      if(this->compute_jlm){
       for(int l=this->Lmin; l<this->Lmax;++l)for(int m=this->Lmin; m<this->Lmax;++m)this->Jlm[l][m]=0;
       cout<<this->compute_jlm<<"  "<<Nrings<<endl;

       for(int ir=0;ir<Nrings;ir++){
         comp_time(start,Nrings,ir);
         double x=cos(theta_new[ir]);
         int Npixels_ring=npix_ring(nside, ir);
         for(int m=0;m<=Lmax;m++){
           double *Pl= new double[this->Lmax-m+1];
           gsl_sf_legendre_sphPlm_array(this->Lmax,m,x,Pl);
           for(int ipr=0;ipr<Npixels_ring;ipr++){
             int ip=ipr+iring;
             if(pixmask[ip]==1){
              int il=0;
              for(int l=m;l<=this->Lmax;l++){
                double Plm=Pl[il];
                this->Jlm[l][m]  += this->area_pixel*pow(Plm,2);
                ++il;
              }
            }
          }
          delete[] Pl;
        }
        iring+=Npixels_ring;
      }

      ofstream jlm_file;
      jlm_file.open(this->output_file_Jlm.c_str());
      for(int l=this->Lmin;l<=this->Lmax;++l){
        for(int m=0;m<=this->Lmax;++m){
          jlm_file<<l<<"\t"<<m<<"\t"<<this->Jlm[l][m]<<endl;
        }
      }
      jlm_file.close();
     }
  
    else if(!this->compute_jlm){
      vector<vector<double> > propJ;
      Fmd.read_file(this->output_file_Jlm,propJ);
      int ik=-1;
      for(int l=this->Lmin;l<=this->Lmax;++l){
        for(int m=0;m<=this->Lmax;++m){
          ik++;this->Jlm[l][m]= propJ[ik][2];
        }
      }
      propJ.clear();
    }
    cout<<"Done"<<RESET<<endl;
    }

}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void Cl_FUNCTIONS::Map2Jlm(Healpix_Map<double>&map_ran){

  time_t start;  
  time (&start);
  for(int i=0;i<Lmax+1;i++)for(int im=0;im<Lmax+1;im++)this->Jlm[i][im]=0;
  int iring=0;
  for(int ir=0;ir<Nrings;ir++){
    comp_time(start,Nrings,ir);
    double x=cos(theta_new[ir]);
    int Npixels_ring=npix_ring(nside, ir);
    for(int m=0;m<=Lmax;m++){
      double *Pl= new double[Lmax-m+1];
      gsl_sf_legendre_sphPlm_array(Lmax,m,x,Pl);
      for(int ipr=0;ipr<Npixels_ring;ipr++){
	int ip=ipr+iring;      
	double cc=cos(double(m)*phi_new[ip]);
	double ss=sin(double(m)*phi_new[ip]);
	int il=0;
	for(int l=m;l<=Lmax;l++){	
	  double Plm=Pl[il];
	  double Yreal=Plm*cc;
	  double Yimag=Plm*ss;
	  this->Jlm[l][m]+= area_pixel*pow(Plm,2)*map_ran[ip];
	  ++il;
	}
      }
      delete[] Pl;
    }
    iring+=Npixels_ring;
  }
}



// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************



// This is the one

void  Cl_FUNCTIONS::Alm2Cl(int iz, int jz,  Alm<xcomplex <double> >&Ilm, vector<double> &Cl){
  double sn=(iz==jz? Shot_Noise[iz]:0);
  double norm_fsky=1./sn; // The mean surface number density
  sn = (this->shot_noise_correction=="yes"? sn:0);
  double FSky = this->type_of_P_estimator=="D" ? 1.0 : this->sky_fraction;

  double Mean_ngal_pix_i=Mean_ngal_pix[iz];
  double Mean_ngal_pix_j=Mean_ngal_pix[jz];


  for(int l=this->Lmin;l<=this->Lmax;l++)this->Wl[l]=0;
  for(int l=this->Lmin;l<=this->Lmax;l++)Cl[l]=0;
  
  if(this->sky == "masked_sky"){
    if(this->sampling=="H"){
      
      for(int l=this->Lmin;l<=this->Lmax;l++){
	for(int m=0;m<=l;m++){

	  double jlm=this->type_of_P_estimator=="D" ? this->Jlm[l][m] : 1.0 ;
	  double Br_z1=(this->Blm[l][m][iz].real()-Mean_ngal_pix_i*Ilm(l,m).real())/Mean_ngal_pix_i;
	  double Bi_z1=(this->Blm[l][m][iz].imag()-Mean_ngal_pix_i*Ilm(l,m).imag())/Mean_ngal_pix_i;
	  double Br_z2=(this->Blm[l][m][jz].real()-Mean_ngal_pix_j*Ilm(l,m).real())/Mean_ngal_pix_j;
	  double Bi_z2=(this->Blm[l][m][jz].imag()-Mean_ngal_pix_j*Ilm(l,m).imag())/Mean_ngal_pix_j;
	  Cl[l]+=(Br_z1*Br_z2+Bi_z1*Bi_z2)/jlm;
	  this->Wl[l]+=norm(Ilm(l,m));
	  if(m==0){Cl[l]/=2.; this->Wl[l]/=2.;}
	}
	
  //      double cor=pow(this-> pixel_window[l],-1);
        Cl[l]=Cl[l]/(FSky*(l+0.5))-sn;
        this->Wl[l]=Wl[l]/(l+0.5);
	this->lvec[l]=(double)l;
      }
    }    
    else{ //if direct sum. Here we are not sure about the implementation of the Ilm, which comes from a Masked. Should we use randoms?
      for(int l=this->Lmin;l<=this->Lmax;l++){
	for(int m=0;m<=l;m++){
	  double jlm=this->type_of_P_estimator=="D" ? this->Jlm[l][m] : 1.0 ;
	  double Br_z1=(this->Blm[l][m][iz].real()-Mean_ngal[iz]*Ilm(l,m).real())/Mean_ngal[iz];
	  double Bi_z1=(this->Blm[l][m][iz].imag()-Mean_ngal[iz]*Ilm(l,m).imag())/Mean_ngal[iz];
	  double Br_z2=(this->Blm[l][m][jz].real()-Mean_ngal[jz]*Ilm(l,m).real())/Mean_ngal[jz];
	  double Bi_z2=(this->Blm[l][m][jz].imag()-Mean_ngal[jz]*Ilm(l,m).imag())/Mean_ngal[jz];
	  Cl[l]+=(Br_z1*Br_z2+Bi_z1*Bi_z2)/jlm;
          this->Wl[l]+=norm(Ilm(l,m));
          if(m==0){Cl[l]/=2.;this->Wl[l]/=2.;}
	}
	Cl[l]=Cl[l]/(FSky*(l+0.5))-sn;
	this->Wl[l]/=(l+0.5);
	this->lvec[l]=(double)l;
      }
    }
  }

  
  else {  // If full sky. In this case, the estimator D and K are equivalent, Jlm=1, Ilm=0
    
    if(this->sampling=="H"){
      
      for(int l=this->Lmin;l<=this->Lmax;l++){
	for(int m=0;m<=l;m++){
	  double Br_z1=this->Blm[l][m][iz].real()/Mean_ngal_pix[iz];
	  double Bi_z1=this->Blm[l][m][iz].imag()/Mean_ngal_pix[iz];
	  double Br_z2=this->Blm[l][m][jz].real()/Mean_ngal_pix[jz];
	  double Bi_z2=this->Blm[l][m][jz].imag()/Mean_ngal_pix[jz];
	  
	  Cl[l]+=(Br_z1*Br_z2+Bi_z1*Bi_z2);
	  this->Wl[l]+=norm(Ilm(l,m));
	  if(m==0)Cl[l]/=2.;
	  if(m==0) this->Wl[l]/=2.;
	}
        Cl[l]=(Cl[l]/(l+0.5))-sn;  // FOR D
        this->Wl[l]=Wl[l]/(l+0.5);
	this->lvec[l]=(double)l;
      }
    }
    
    else{ //if Direct summation && full sky
      for(int l=this->Lmin;l<=this->Lmax;l++){
	for(int m=0;m<=l;m++){
	  double Br_z1=this->Blm[l][m][iz].real();
	  double Bi_z1=this->Blm[l][m][iz].imag();
	  double Br_z2=this->Blm[l][m][jz].real();
	  double Bi_z2=this->Blm[l][m][jz].imag();
	  
	  Cl[l]+=(Br_z1*Br_z2+Bi_z1*Bi_z2);
	  this->Wl[l]+=norm(Ilm(l,m));
	  if(m==0)Cl[l]/=2.;
	  if(m==0) this->Wl[l]/=2.;
	}
	Cl[l]=Cl[l]/(l+0.5)/pow(norm_fsky,2)-sn;  // FOR D
        this->Wl[l]/=(l+0.5);  // Applied only to the window function correctly?
	this->lvec[l]=(double)l;
      }
    }
  }
}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void Cl_FUNCTIONS::Cl_bins(vector<double>Cl,  vector<double>&Clbin){

  cout<<BLUE<<endl;
  if(this->bin_type=="linear"){ 
    for(int i=0;i<(signed)lbin.size();i++)this->nmodes[i]=0;
    for(int i=0;i<(signed)lbin.size();i++)Clbin[i]=0;
    
    double deltal=((double)(this->Lmax-this->Lmin))/((double)this->lbin.size()); 
    for(int l=this->Lmin;l<=this->Lmax;l++){
      int il=(int)floor((double)(l-this->Lmin)/((double)deltal));
      if(il==lbin.size())il--;
      Clbin[il]+=(l+0.5)*Cl[l];
      this->nmodes[il]+=(l+0.5);
    }
    for(int i=0;i<lbin.size();++i)Clbin[i]/=this->nmodes[i];
  }
  else{
    if(this->bin_type=="log"){
      double deltal=log10(Lmax/1.0)/((double)this->lbin.size());
      for(int i=0;i<lbin.size();++i)this->nmodes[i]=0;
      for(int l=0;l<=Lmax;l++){
        int il=(int)floor((log10(l)-log10(1.0))/deltal);
        if(il==lbin.size())il--;
        Clbin[il]+=(l+0.5)*Cl[l];
	this->nmodes[il]+=(l+0.5);
      }
      for(int i=0;i<lbin.size();++i)Clbin[i]/=this->nmodes[i];
    }
  }
  cout<<RESET;
}

// *******************************************************************************************************************************************************
void Cl_FUNCTIONS::eCl_bins(vector<double>Cl,  vector<double>&Clbin){

  cout<<BLUE<<endl;
  if(this->bin_type=="linear"){
    for(int i=0;i<(signed)lbin.size();i++)this->nmodes[i]=0;
    for(int i=0;i<(signed)lbin.size();i++)Clbin[i]=0;

    double deltal=((double)(this->Lmax-this->Lmin))/((double)this->lbin.size());
    for(int l=this->Lmin;l<=this->Lmax;l++){
      int il=(int)floor((double)(l-this->Lmin)/((double)deltal));
      if(il==lbin.size())il--;
      Clbin[il]+=pow((l+0.5)*Cl[l],2);
      this->nmodes[il]+=(l+0.5);
    }
    for(int i=0;i<lbin.size();i++)Clbin[i]=sqrt(Clbin[i])/this->nmodes[i];
  }
  else{
    if(this->bin_type=="log"){
      double deltal=log10(Lmax/1.0)/((double)this->lbin.size());
      for(int i=0;i<lbin.size();++i)this->nmodes[i]=0;
      for(int l=0;l<=Lmax;l++){
        int il=(int)floor((log10(l)-log10(1.0))/deltal);
        if(il==lbin.size())il--;
        Clbin[il]+=pow(Cl[l],2);
    this->nmodes[il]++;
      }
      for(int i=0;i<lbin.size();i++)Clbin[i]=sqrt(Clbin[i])/this->nmodes[i];
    }
  }
  cout<<RESET;
}
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************



// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void Cl_FUNCTIONS::W3J(){
  cout<<RED<<"Computing Wigner Symbols"<<endl;
  Wigner3J.resize(Lmax+1);
  for(int i=0;i<Lmax+1;i++){
    Wigner3J[i].resize(Lmax+1);
    for(int j=0;j<Lmax+1;j++){
      Wigner3J[i][j].resize(Lmax+1);
    }
  }

  time_t start;  
  time (&start);
  for(int i=0;i<Lmax+1;i++){
    for(int j=0;j<Lmax+1;j++){
      for(int k=0;k<Lmax+1;k++){
	int J=i+j+k;
	if(abs(i-j) <= k &&  k<=i+j){
	  if(J%2==0){
	    Wigner3J[i][j][k]= sqrt(factorial(J-2*i)*factorial(J-2*j)*factorial(J-2*k)/factorial(J+1))*(factorial(J/2))/(factorial(J/2-i)*factorial(J/2-j)*factorial(J/2-k));
	  }
	  else{
	    Wigner3J[i][j][k]=0;
	  }
	}
      }
    }
  }
  
  cout<<"DONE"<<RESET<<endl;								       
  
}
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void Cl_FUNCTIONS::Mixing_matrix(){
  cout<<"Computing the mixing matrix ..."<<endl;
  for(int l1=Lmin;l1<Lmax;l1++){
    for(int l2=Lmin;l2<Lmax+1;l2++){
      double Rbis=0;
      for(int l3=Lmin;l3<Lmax+1;l3++){
       	if(abs(l1-l2) <= l3 &&  l3<=l1+l2){
	  double WC = Wigner3J[l1][l2][l3];
	  Rbis+=(2.*l3+1)*this->Wl[l3]*(WC*WC);
	}
      }
      this->R[l1][l2]=((2.*l2+1)/(4.*M_PI))*Rbis;
    }
  }
  cout<<"Done"<<RESET<<endl;
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************


void Cl_FUNCTIONS::Rll_bins(){


  // Only available for the Rll written in term of the Wigner symbols

  double deltal;
  
  if(this->bin_type=="linear"){ 
    
    for(int i=0;i<(signed)lbin.size();i++)this->nmodes[i]=0;
    
    deltal=(this->Lmax-this->Lmin)/((double)lbin.size()); //here and just here I define N_L_bins as a double instead of an integer
    cout<<RESET;

    for(int l1=this->Lmin;l1<this->Lmax;++l1){
      int lbina=floor((l1-this->Lmin)/deltal);

      if(lbina==lbin.size())lbina--;
      this->nmodes[lbina]+=(2.*l1+1);
      for(int l2=this->Lmin;l2<this->Lmax+1;++l2){
    	double Rbis=0;
    	for(int l3=this->Lmin;l3<this->Lmax+1;l3++){
    	  if(abs(l1-l2) <= l3 &&  l3<=l1+l2){
    	    double WC = this->Wigner3J[l1][l2][l3];
  	    Rbis+=(2*l1+1)*(2.*l3+1)*this->Wl[l3]*(WC*WC);
  	  }
  	}
        this->Rll_bin[lbina][l2]+=((2.*l2+1)/(4.*M_PI))*Rbis;
      }
    }
    for(int i=0;i<lbin.size();i++)for(int l2=this->Lmin;l2<this->Lmax+1;l2++)Rll_bin[i][l2]/=this->nmodes[i];   }
  else{
    if(this->bin_type=="log"){
      
      deltal=log10(Lmax/1.0)/((double)lbin.size());
      cout<<RESET;
      for(int i=0;i<lbin.size();i++)this->nmodes[i]=0;
      
      for(int l1=this->Lmin;l1<this->Lmax;l1++){
  	int lbina=floor((log(l1)-log(1.0))/deltal);
	this->nmodes[lbina]+=(2.*l1+1);
  	for(int l2=this->Lmin;l2<this->Lmax+1;l2++){
  	  double Rbis=0;
  	  for(int l3=this->Lmin;l3<this->Lmax+1;l3++){
  	    if(abs(l1-l2) <= l3 &&  l3<=l1+l2){
  	      double WC = this->Wigner3J[l1][l2][l3];
              Rbis+=(2.*l1+1)*(2.*l3+1)*this->Wl[l3]*(WC*WC);
  	    }
  	  }
          this->Rll_bin[lbina][l2]+=((2.*l2+1)/(4.*M_PI))*Rbis;
  	}
      }
      
      for(int i=0;i<lbin.size();i++)for(int l2=this->Lmin;l2<this->Lmax+1;l2++)Rll_bin[i][l2]/=this->nmodes[i]; 
    }
  }
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void Cl_FUNCTIONS::get_mixing_matrix(){

  // compute the wigner symbols
  this->W3J();

  this->Mixing_matrix();
  double factor_fs = (this->mixing_matrix_exact? 1.0:this->sky_fraction);
  for(int l1=this->Lmin;l1<=this->Lmax;l1++)for(int l2=this->Lmin;l2<=this->Lmax;l2++)this->R[l1][l2]/=factor_fs;
  for(int l2=this->Lmin;l2<=this->Lmax;l2++){
    vector<double> slice (this->Lmax+1,0);
    string rfile=this->output_file_window+"_MixingMatrix_l_"+to_string(l2)+".dat";
    for(int l1=this->Lmin;l1<this->Lmax;l1++)slice[l1]=this->R[l1][l2];
    this->Fmi.write_to_file(rfile, this->lvec, slice);
  }
  string rfile;
  if(this->mixing_matrix_exact)rfile=this->output_file_window+"_MixingMatrix_exact_nside_"+to_string(this->nside)+".dat";
  else rfile=this->output_file_window+"_MixingMatrix_nside_"+to_string(this->nside)+".dat";
  this->Fmi.write_to_file(rfile, this->lvec, this->lvec, this->R);

  this->Rll_bins();
  for(int l1=0;l1<this->N_L_bins;l1++)for(int l2=this->Lmin;l2<=this->Lmax;l2++)this->Rll_bin[l1][l2]/=factor_fs;
  rfile=this->output_file_window+"_MixingMatrix_lbins_"+to_string(this->N_L_bins)+"_nside_"+to_string(this->nside)+".dat";
  this->Fmi.write_to_file(rfile, this->lbin, this->lvec, this->Rll_bin);


  for(int l1=0;l1<this->N_L_bins;l1++){
    vector<double> slice3 (this->Lmax+1,0);
    for(int l2=this->Lmin;l2<=this->Lmax;l2++)slice3[l2]=this->Rll_bin[l1][l2];
    rfile=this->output_file_window+"_MixingMatrix_lbin_"+to_string(l1)+"_nside_"+to_string(this->nside)+".dat";
    this->Fmi.write_to_file(rfile, this->lvec,slice3);
  }

}





// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// ######################################################################
// ######################################################################
// ######################################################################
// ######################################################################
double Cl_FUNCTIONS::get_min(char c, int ii){

  vector<double> pro=(c=='d'? this->  prop:this->  prop_r);
  int ncols=(c=='d'? this->n_columns_gal:this->n_columns_ran);
  double ZX=10000;
  double zmx;
  for(int i=0;i<pro.size()/ncols;i++){
    zmx=min(ZX,pro[ii+ncols*i]);
    ZX=zmx;
  }
  return zmx;
}
// ######################################################################
// ######################################################################
double Cl_FUNCTIONS::get_max(char c, int ii ){
  vector<double> pro=(c=='d'? this->  prop  :this->  prop_r);
  int ncols=(c=='d'? this->n_columns_gal:this->n_columns_ran);
  double ZX=-100;
  double zmx;
  for(int i=0;i<pro.size()/ncols;i++){
    zmx=max(ZX,pro[ii+ncols*i]);
    ZX=zmx;
  }
  return zmx;
}

// ######################################################################
// ######################################################################
void Cl_FUNCTIONS::read_input_cats(string cat, string file){
   if(cat=="g")this->n_columns_gal=Fmd.read_file2(file, this->  prop);
   else if(cat=="r")this->n_columns_ran=Fmd.read_file2(file, this->  prop_r);
}



// ######################################################################
// ######################################################################

void Cl_FUNCTIONS::get_zbins_same_ngal(int nbins, int iiz, double zmn, double zmx, vector<vector<double> >&zzv){

  // ************************************
  // Construct the dNdz with many bins
  int na=20000; //This value is critical to avoid seg foults. The higher, the best.
  double delta=(zmx-zmn)/((double)na);
  vector<int>dn(na,0);
  for(int i=0;i<this->  prop.size()/this->n_columns_gal;++i){
    if(this->  prop[iiz+this->n_columns_gal*i]<zmx && this->  prop[iiz+this->n_columns_gal*i]>=zmn){
      int iza=floor((this->  prop[iiz+this->n_columns_gal*i]-zmn)/delta);
      dn[iza]++;
    }
  }

  int Nca=(int)(floor)(  prop.size()/this->n_columns_gal/((double)nbins)); //Desired number of galaxies per redshift bin:
  cout<<"Number of galaxies per redshift bin = "<<Nca<<endl;
  vector<double>zan(na,0);
  for(int i=0;i<dn.size();++i)zan[i]=zmn+(i+0.5)*delta;

  // Set the full z-interval
  zzv[0][0]=zmn;
  zzv[0][1]=zmx;
  if(nbins==1){
    zzv[1][0]=zmn;
    zzv[1][1]=zmx;
  }

  if(nbins>1){
    for(int ib=1;ib<=nbins;++ib){
      int caa=0;
      vector<double>zaux;
      for(int i=0;i<dn.size();++i){
        caa+=dn[i];  //Cumulative number of galaxies
        if((caa>=Nca*(ib-1)) &&  (caa<Nca*ib)) zaux.push_back(zan[i]);
      }
      zzv[ib][0]=zaux[0]-0.5*delta;  //Allocate the zmin of the ib zbin
      zzv[ib][1]=zaux[zaux.size()-1]+0.5*delta; //Allocate the zmax of the ib zbin
      zaux.clear();
    }
  }
 dn.clear();
 return ;
}


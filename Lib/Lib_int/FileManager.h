#ifndef __FILE_MANAGER__
#define __FILE_MANAGER__

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
# include <omp.h>
using namespace std;

template<class T> class FILE_MANAGER{
 protected:
  string line;
  ifstream dina;
  ofstream out; 

public:
  FILE_MANAGER(){}
  ~FILE_MANAGER(){}

  void read_file(string, vector<vector<T> > &);

  int read_file2(string,vector<T> &);
  
  void write_to_file(string, vector<T>& , vector<double> &);

  void write_to_file(string, vector<T>& , vector<double> &, vector<double> &);

  void write_to_file(string, vector<T>& , vector<double> &, vector<double> &,vector<double>&);

  void write_to_file(string, vector<T>&,  vector<T> &, vector<T> &, vector<vector<double> >&);

  void write_to_file(string, vector<T>&,  vector<T> &, vector<vector<double> >&); 

  void write_to_file(string, vector<double>,  vector<T>, vector<vector<double> >);

  void write_to_file_n(string fname, vector<T> &, vector<T> &, vector<vector<double> >&);

  void write_to_file(string, vector<T>& , vector<T> &, vector<T> &, vector<T> &, vector<T> &,vector<double> &);

  void write_to_file(string, vector<T>& , vector<T> &, vector<T> &, vector<T> &, vector<double> &);

  void write_to_file_nn(string, vector<double>& , vector<T> &, vector<double> &);

};



// *********************************************************************************************
// **********************************************************************************************

template<class T> void FILE_MANAGER<T>::read_file(string fname, vector< vector<T> > &prop_ob){
  time_t start;
  time(&start);
  
  dina.open(fname.c_str(), ios::in);
  assert(dina);
  T aa;
  cout<<CYAN<<"*********************************************************"<<RESET<<endl;
  cout<<CYAN<<"*********************************************************"<<RESET<<endl;
  cout<<BLUE<<"Reading input file"<<RESET<<endl;
  cout<<GREEN<<fname<< " . . . "<<RESET<<endl;
  while(getline (dina, line)){
    stringstream ss(line);
    vector <T> nc;
    while(ss>>aa)nc.push_back (aa); //->vector for each galaxy with properties
    prop_ob.push_back (nc);         //->fill array with all galaxies
  }
  cout<<BLUE<<"with "<<prop_ob.size()<<" lines and "<<prop_ob[0].size()<<" columns "<<RESET<<endl;
  dina.clear();
  dina.close();

  cout<<RED;
  time_t end; time (&end);
  double diff = difftime(end,start);
  if (diff<60) cout <<"Lapse: "<<diff<<" secs"<<endl;
  else if (diff<3600) cout <<"Lapse: "<<diff/60<<" minutes"<<endl;
  else cout<<"Lapse: "<<diff/3600<<" hours"<<endl;
  cout<<RESET;

}

template<class T> int FILE_MANAGER<T>::read_file2(string fname, vector<T>&prop_ob){

  time_t start;
  time(&start);

  int items_per_line=0;
  string line;
  double aa;
  ifstream dina;
  dina.open(fname.c_str(), ios::in);
  assert(dina);

  
  if(!dina.eof()) {
    getline (dina, line);
    stringstream ss(line);
    while(ss>>aa)items_per_line++;
  }
  dina.seekg(0, dina.beg);
  

  cout<<CYAN;
  cout<<"Reading input ile  "<<fname<<endl;
  while(getline (dina, line)){
    stringstream ss(line);
    while(ss>>aa)prop_ob.push_back (aa); //->vector for each galaxy with properties
  }
  
  cout<<"Number of lines = "<<prop_ob.size()/items_per_line<<endl;
  cout<<"Number of columns = "<<items_per_line<<endl;
  dina.clear();
  dina.close();

  cout<<RED;
  time_t end; time (&end);
  double diff = difftime(end,start);
  if (diff<60) cout <<"Lapse: "<<diff<<" secs"<<endl;
  else if (diff<3600) cout <<"Lapse: "<<diff/60<<" minutes"<<endl;
  else cout<<"Lapse: "<<diff/3600<<" hours"<<endl;
  cout<<RESET;

  return items_per_line;

}

// *********************************************************************************************
// **********************************************************************************************


template<class T>  void FILE_MANAGER<T>::write_to_file(string fname, vector<T>& kve, vector<double> &pk)
{
  
  ofstream out; 
  out.open(fname.c_str() , ios::out); 
  out.precision(10); 
  out.setf(ios::showpoint); 
  out.setf(ios::scientific); 
  out.width(12); 
  for(int i=0;i<kve.size();i++)out<<kve[i]<<"\t"<<pk[i]<<endl;
  out.close(); 
  cout<<BLUE<<"Output file:"<<RESET<<endl;
  cout<<GREEN<<fname<<RESET<<endl;
  cout<<BLUE<<"with "<<kve.size()<<" lines"<<RESET<<endl;
  cout<<"***********************************************************************************"<<endl;
  return; 
}

// *********************************************************************************************
// **********************************************************************************************


template<class T>  void FILE_MANAGER<T>::write_to_file(string fname, vector<T>& kve, vector<double> &pk, vector<double> &pt)
{
  
  ofstream out; 
  out.open(fname.c_str() , ios::out); 
  out.precision(10); 
  out.setf(ios::showpoint); 
  out.setf(ios::scientific); 
  out.width(12); 
  for(int i=0;i<kve.size();i++)out<<kve[i]<<"\t"<<pk[i]<<"\t"<<pt[i]<<endl;
  out.close(); 
  cout<<CYAN;
  cout<<"Output file:"<<endl;
  cout<<fname<<endl; 
  cout<<"with "<<kve.size()<<" lines"<<endl; 
  cout<<"***********************************************************************************"<<endl;
  cout<<RESET;
  return; 
}

// *********************************************************************************************
// **********************************************************************************************

template<class T>  void FILE_MANAGER<T>::write_to_file(string fname, vector<T>& kve, vector<double> &pk, vector<double> &pt, vector<double> &pta)
{
  
  ofstream out; 
  out.open(fname.c_str() , ios::out); 
  out.precision(10); 
  out.setf(ios::showpoint); 
  out.setf(ios::scientific); 
  out.width(12); 
  for(int i=0;i<kve.size();i++)out<<kve[i]<<"\t"<<pk[i]<<"\t"<<pt[i]<<"\t"<<pta[i]<<endl;
  out.close(); 
  cout<<CYAN;
  cout<<"Output file:"<<endl;
  cout<<fname<<endl; 
  cout<<"with "<<kve.size()<<" lines"<<endl; 
  cout<<"***********************************************************************************"<<endl;
  cout<<RESET;
  return; 
}

// *********************************************************************************************
// **********************************************************************************************

template <class T> void FILE_MANAGER<T>::write_to_file(string fname, vector<T> &xx, vector<T> &yy, vector<vector<double> >&xy){
  out.open(fname.c_str() , ios::out);
  out.precision(10);
  out.setf(ios::showpoint);
  out.setf(ios::scientific);
  out.width(12);
  for(int i=0;i<xx.size();i++)for(int j=0;j<yy.size();j++)out<<xx[i]<<"\t"<<yy[j]<<"\t"<<xy[i][j]<<endl;
  out.close();
  cout<<CYAN;
  cout<<"Output file:"<<endl;
  cout<<fname<<endl;
  cout<<"with "<<xx.size()<<" lines"<<endl;
  cout<<"***********************************************************************************"<<endl;
  cout<<RESET;
}


template <class T>  void FILE_MANAGER<T>::write_to_file_n(string fname, vector<T> &xx, vector<T> &yy, vector<vector<double> >&xy){
  out.open(fname.c_str() , ios::out);
  out.precision(10);
  out.setf(ios::showpoint);
  out.setf(ios::scientific);
  out.width(12);
  for(int i=0;i<xx.size();i++)for(int j=0;j<yy.size();j++)out<<xx[i]<<"\t"<<yy[j]<<"\t"<<xy[i][j]<<endl;
  out.close();
  cout<<CYAN;
  cout<<"Output file:"<<endl;
  cout<<fname<<endl;
  cout<<"with "<<xx.size()<<" lines"<<endl;
  cout<<"***********************************************************************************"<<endl;
  cout<<RESET;
}


// **********************************************************************************************
// **********************************************************************************************


template <class T> void FILE_MANAGER<T>::write_to_file(string fname, vector<double> xx, vector<T> yy, vector<vector<double> >xy){
  out.open(fname.c_str() , ios::out);
  out.precision(10);
  out.setf(ios::showpoint);
  out.setf(ios::scientific);
  out.width(12);
  for(int i=0;i<xx.size();i++)for(int j=0;j<yy.size();j++)out<<xx[i]<<"\t"<<yy[j]<<"\t"<<xy[i][j]<<endl;
  out.close();
  cout<<CYAN;
  cout<<"Output file:"<<endl;
  cout<<fname<<endl;
  cout<<"with "<<xx.size()<<" lines"<<endl;
  cout<<"***********************************************************************************"<<endl;
  cout<<RESET;
}

// **********************************************************************************************
// **********************************************************************************************

template <class T> void FILE_MANAGER<T>::write_to_file(string fname, vector<T> &xx, vector<T> &yy, vector<T> &zz, vector<vector<double> >&xy){
  out.open(fname.c_str() , ios::out);
  out.precision(10);
  out.setf(ios::showpoint);
  out.setf(ios::scientific);
  out.width(12);
  for(int i=0;i<xx.size();i++)for(int j=0;j<yy.size();j++)out<<xx[i]<<"\t"<<yy[j]<<"\t"<<zz[j]<<"\t"<<xy[i][j]<<endl;
  out.close();
  cout<<CYAN;
  cout<<"Output file:"<<endl;
  cout<<fname<<endl;
  cout<<"with "<<xx.size()<<" lines"<<endl;
  cout<<"***********************************************************************************"<<endl;
  cout<<RESET;
}


// *********************************************************************************************
// **********************************************************************************************


template<class T>  void FILE_MANAGER<T>::write_to_file(string fname, vector<T>& kve, vector<T> &p1, vector<T> &p2, vector<T> &p3, vector<T> &p4, vector<double> &p5 )
{
  
  out.open(fname.c_str() , ios::out); 
  out.precision(10); 
  out.setf(ios::showpoint); 
  out.setf(ios::scientific); 
  out.width(12); 
  for(int i=0;i<kve.size();i++)out<<kve[i]<<"\t"<<p1[i]<<"\t"<<p2[i]<<"\t"<<p3[i]<<"\t"<<p4[i]<<"\t"<<p5[i]<<endl;
  out.close();
  cout<<CYAN;
  cout<<"Output file:"<<endl;
  cout<<fname<<endl; 
  cout<<"with "<<kve.size()<<" lines"<<endl; 
  cout<<"***********************************************************************************"<<endl;
  cout<<RESET;
  return; 
}


// *********************************************************************************************
// **********************************************************************************************


template<class T>  void FILE_MANAGER<T>::write_to_file(string fname, vector<T>& kve, vector<T> &p1, vector<T> &p2, vector<T> &p3, vector<double> &p4){
  
  out.open(fname.c_str() , ios::out); 
  out.precision(10); 
  out.setf(ios::showpoint); 
  out.setf(ios::scientific); 
  out.width(12); 
  for(int i=0;i<kve.size();i++)out<<kve[i]<<"\t"<<p1[i]<<"\t"<<p2[i]<<"\t"<<p3[i]<<"\t"<<p4[i]<<endl;
  out.close(); 
  cout<<CYAN;
  cout<<"Output file:"<<endl;
  cout<<fname<<endl; 
  cout<<"with "<<kve.size()<<" lines"<<endl; 
  cout<<"***********************************************************************************"<<endl;
  cout<<RESET;
  return; 
}

// *********************************************************************************************
// **********************************************************************************************

template<class T>  void FILE_MANAGER<T>::write_to_file_nn(string fname, vector<double>& xx, vector<T> &p1, vector<double> &p2)
{

  out.open(fname.c_str() , ios::out);
  out.precision(10);
  out.setf(ios::showpoint);
  out.setf(ios::scientific);
  out.width(12);
  for(int i=0;i<xx.size();i++)out<<xx[i]<<"\t"<<p1[i]<<"\t"<<p2[i]<<endl;
  out.close();
  cout<<"Output file:"<<endl;
  cout<<fname<<endl;
  cout<<"with "<<xx.size()<<" lines"<<endl;
  cout<<"***********************************************************************************"<<endl;
  return;
}




// **********************************************************************************************
// **********************************************************************************************
#endif

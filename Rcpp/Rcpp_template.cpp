///////////////////
// C++ libraries //
///////////////////

#include <Rcpp.h>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <ctype.h>      // For isalnum()
#include <bits/stdc++.h> 

using namespace Rcpp;
using namespace std;


///////////////////////
// Declare variables //
///////////////////////
typedef vector<string> stringList;
typedef vector<int> numList;



///////////////////////
// Declare functions //
///////////////////////

// Given a string 's' and a character 'c',
// counts the number of times 'c' appears in the input string 's'
int count_char_in_string(string s, string c) {
  
  int count = 0;
  
  // Iterate over the input string
  for (int i = 0; i < s.size(); i++) {
    
    // Increase the counter when if finds the input character 'c'
    if (s[i] == c[0]) {
      count++;
    }
  }
  
  return count;
}

// Factorial function: due to the restriction of data type, this supports as max factorial(20)
unsigned long long int factorial(int n) {
  if (n == 0 || n == 1)
    return 1;
  else
    return n * factorial(n - 1);
}

// Logarithm of any base
double log_base(double base, double x) {
  return log(x) / log(base);
}


// This is the core function to calculate the WF-score
// The input is one sequence
double WF(const string x){
  
  double ws_score = 0;
  double N        = x.size();
  int alpha_size  = 4;  // DNA alphabet
  
  // Count the number of times each nucleotide appears in the input strin
  int As = count_char_in_string(x, "A");
  int Cs = count_char_in_string(x, "C");
  int Gs = count_char_in_string(x, "G");
  int Ts = count_char_in_string(x, "T");
  
  // cout << "A:" << As << "\n";
  // cout << "C:" << Cs << "\n";
  // cout << "G:" << Gs << "\n";
  // cout << "T:" << Ts << "\n";
  
  // Wottoon-Federhen formula
  unsigned long long int Nt_prod  = factorial(As) * factorial(Cs) * factorial(Gs) * factorial(Ts);
  // cout << "Nt product:" << Nt_prod << "\n";
  ws_score        = log_base(alpha_size, factorial(N) / Nt_prod) / N;
  
  return(ws_score);
}




/////////////////////////////
// Vectorize core function //
/////////////////////////////

// [[Rcpp::export]]
// This is the 'skeleton' to vectorize the function
// The input must be one array, then the core function is applied to each element
// This is the function called from the R code, indicated by [[Rcpp::export]]
vector<double> WF_complexity_score(stringList s1)
{
  
  // Number of input elements
  int size1 = s1.size();
  
  // Iterate over each element
  vector<double> dd(size1);
  for (int i = 0; i < size1; i++) {
    
    // Calculate the WF score on each string
    double d = WF(s1[i]);
    dd[i] = d;
  }
  
  return dd;
  
}

// The idea of using Rcpp code is to create a small function that runs much faster
// in C++ than in R, it is necessary to create an entire code, just the functions
// that need to be called from R.
//
// In my experience, the best solution is to create that core function and finally 
// create a vectorized version of the core function, i.e., call the function through
// a for loop to apply the function in vectors
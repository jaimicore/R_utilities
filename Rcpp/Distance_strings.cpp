#include <Rcpp.h>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <ctype.h>      // For isalnum()
#include <bits/stdc++.h> 

using namespace Rcpp;
using namespace std;

typedef vector<string> stringList;
typedef vector<int> numList;


vector<char> Char2Vec(const string& str)
{
  vector<char> result;
  
  // For each character in the string
  for (char ch : str)
  {
    // Copy only alphabetical characters and numeric digits
    if (isalnum(ch))
    {
      result.push_back(ch);
    }
  }
  
  return result;
}



double CalcDist(const string& x, const string& y){
  
  vector<char> xx = Char2Vec(x);
  vector<char> yy = Char2Vec(y);

  
  double match    = 0;
  double distance = 0;
  for( int i=0; i < xx.size(); i++) {
    if (xx[i] == yy[i]) {
      match += 1;
    }
  }

  distance = match/xx.size();

  return(distance);
}



// [[Rcpp::export]]
vector<double> DistanceStrings(stringList s1, stringList s2)
{
  
  // Number of input elements
  int size1 = s1.size();
  int size2 = s2.size();
  
  
  // This works //
  // NumericVector string_split(stringList s1, stringList s2)
  // Rcpp::NumericVector xx(2);
  // xx[0] = size1;
  // xx[1] = size2;
  // return xx;
  
  
  // Iterate over each pair of elements
  vector<double> dd(size1);
  for (int i = 0; i < size1; i++)
  {

    // Extract the string from the list
    String string1 = s1[i];
    String string2 = s2[i];
    

    // Calculate the distance on each pair of strings
    double d = CalcDist(string1, string2);
    dd[i] = d;
  }

  return dd;
  
}
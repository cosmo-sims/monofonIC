#pragma once

#include <string>
#include <sstream>
#include <iomanip>


template <typename T>
std::string str(const T &x)
{
  std::ostringstream os;
  os << x;
  return os.str();
}
// specilization
template <>
std::string str(const float &x)
{
  std::ostringstream os;
  os << std::setprecision(8) << x;
  return os.str();
}
template <>
std::string str(const double &x)
{
  std::ostringstream os;
  os << std::setprecision(16) << x;
  return os.str();
}
template <>
std::string str(const bool &x)
{
  {
    return x ? "yes" : "no";
  }
}

template <>
std::string str(const std::string &x) { return x; }

std::string str(const char *s) { return std::string(s); }

class ClassParams{
public:

  ClassParams(){};
  ClassParams( const ClassParams& o):pars(o.pars){};

  //use this to add a CLASS variable
  template<typename T> unsigned add(const std::string& key,const T& val){
  pars.push_back(std::make_pair(key,str(val)));
  return pars.size();
  }

  //accesors
  inline unsigned size() const {return pars.size();}
  inline std::string key(const unsigned& i) const {return pars[i].first;}
  inline std::string value(const unsigned& i) const {return pars[i].second;}


private:
  std::vector<std::pair<std::string,std::string> > pars;
};

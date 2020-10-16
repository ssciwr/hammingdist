#include<Eigen/Dense>

#include<iostream>

int encoded[4] = 
  { 1<<0, 1<<1, 1<<2, 1<<3 };

struct Genome
{
  Genome() : val(0) {}
  Genome(const int& val_) : val(val_) {}
  Genome& operator=(const int& val_)
  {
    val = val_;
    return *this;
  }

  inline Genome operator* (const Genome& other) const
  {
    // This seems to compile to pretty efficient assembly
    // https://godbolt.org/z/nPMxeG
    return ((val & other.val) == 0) ? 1 : 0;
  }

  inline Genome operator+ (const Genome& other) const
  {
    return Genome(this->val + other.val);
  }

  Genome& operator+= (const Genome& other)
  {
    val += other.val;
    return *this;
  }

  int val;
};

template<typename Stream>
Stream& operator<<(Stream& stream, const Genome& g)
{
  stream << g.val;
  return stream;
}

namespace Eigen {

  template<>
  class NumTraits<Genome>
  {
    public:
    using Real = Genome;
    using Literal = Genome;
    
    enum {
      IsComplex = 0,
      IsInteger = 0,
      IsSigned = 0,
      RequireInitialization = 1,
      ReadCost = 1,
      AddCost = 3,
      MulCost = 3
    };

    static constexpr int digits10()
    {
      return 0;
    }
  };
}


int main(int argc, char** argv)
{
  Eigen::Matrix<Genome, 2, 2> data;
  
  data(0, 0) = encoded[0];
  data(0, 1) = encoded[1];
  data(1, 0) = encoded[0];
  data(1, 1) = encoded[0];

  Eigen::Matrix<Genome, 2, 2> prod;
  std::cout << data * data.transpose() << std::endl;
  
  return 0;
}

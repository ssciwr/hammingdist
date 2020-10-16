#include<blaze/Math.h>
#include<iostream>

char encoded[4] = 
  { 1<<0, 1<<1, 1<<2, 1<<3 };

struct Genome
{
  Genome() : val(0) {}
  Genome(const char& val_) : val(val_) {}
  Genome& operator=(const char& val_)
  {
    val = val_;
    return *this;
  }

  inline int operator* (const Genome& other) const
  {
    // This seems to compile to pretty efficient assembly
    // https://godbolt.org/z/nPMxeG
    return ((val & other.val) == 0) ? 1 : 0;
  }

  char val;
};

template<typename Stream>
Stream& operator<<(Stream& stream, const Genome& g)
{
  stream << g.val;
  return stream;
}

int main(int argc, char** argv)
{
  blaze::StaticMatrix<Genome, 2, 2> data;
  
  data(0, 0) = Genome(encoded[0]);
  data(0, 1) = Genome(encoded[1]);
  data(1, 0) = Genome(encoded[0]);
  data(1, 1) = Genome(encoded[0]);

  auto prod = data * blaze::trans(data);

  std::cout << prod << std::endl;
  
  return 0;
}

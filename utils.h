#include <cassert>
#include <iostream>
#include <vector>

int checkResults(const std::vector<int>&    r,
                 const std::vector<int>&    c,
                 const std::vector<double>& v,
                 const std::vector<int>&    p,
                 const std::vector<int>&    r_true,
                 const std::vector<int>&    c_true,
                 const std::vector<double>& v_true,
                 const std::vector<int>&    p_true)
{
  int failures = 0;

  // Check permutation
  assert(p.size() == p_true.size());
  {
    for (size_t i = 0; i < p.size(); ++i)
    {
      if (p[i] != p_true[i])
      {
        std::cout << "ERROR: Permutation mismatch at "
                  << i << ": " << p[i] << " vs " << p_true[i] << "\n";
        failures++;
      }
    }
  }

  // Check rows
  assert(r.size() == r_true.size());
  {
    for (size_t i = 0; i < r.size(); ++i)
    {
      if (r[i] != r_true[i])
      {
        std::cout << "ERROR: Row mismatch at "
                  << i << ": " << r[i] << " vs " << r_true[i] << "\n";
        failures++;
      }
    }
  }

  // Check columns
  assert(c.size() == c_true.size());
  {
    for (size_t i = 0; i < c.size(); ++i)
    {
      if (c[i] != c_true[i])
      {
        std::cout << "ERROR: Column mismatch at "
                  << i << ": " << c[i] << " vs " << c_true[i] << "\n";
        failures++;
      }
    }
  }

  // Check values
  assert(v.size() == v_true.size());
  {
    for (size_t i = 0; i < v.size(); ++i)
    {
      if (std::abs(v[i] - v_true[i]) > 1e-10)
      {
        std::cout << "ERROR: Value mismatch at "
                  << i << ": " << v[i] << " vs " << v_true[i] << "\n";
        failures++;
      }
    }
  }

  return failures;
}

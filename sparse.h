#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

class CooMatrix
{
public:
  CooMatrix(int m, int n)
  {
    m_ = m;
    n_ = n;
  }

  ~CooMatrix() = default;

  void addEntries(const std::vector<int>&    r,
                  const std::vector<int>&    c,
                  const std::vector<double>& v)
  {
    int nnz = r.size();
    assert(nnz == static_cast<int>(c.size())
           && nnz == static_cast<int>(v.size()));

    for (int i = 0; i < nnz; ++i)
    {
      nnz_dup_++;

      std::pair<int, int> key = {r[i], c[i]};

      auto it = entries_.find(key);
      if (it == entries_.end())
      {
        // New entry, add to arrays
        int idx       = nnz_;
        entries_[key] = idx;

        r_.push_back(r[i]);
        c_.push_back(c[i]);
        v_.push_back(v[i]);

        perm_.push_back(idx);
        nnz_++;
      }
      else
      {
        // Duplicate entry, sum the value
        int idx = it->second;
        v_[idx] += v[i];
        perm_.push_back(idx);
      }
    }
  }

  void update(const std::vector<double>& vnew)
  {
    std::fill(v_.begin(), v_.end(), 0.0);

    for (int i = 0; i < perm_.size(); ++i)
    {
      int idx = perm_[i];
      v_[idx] += vnew[i];
    }
  }

  const std::vector<int>& getRows() const
  {
    return r_;
  }

  const std::vector<int>& getColumns() const
  {
    return c_;
  }

  const std::vector<double>& getValues() const
  {
    return v_;
  }

  int getNnz() const
  {
    return nnz_;
  }

  int getNnzDuplicates() const
  {
    return nnz_dup_;
  }

  const std::vector<int>& getPermutation() const
  {
    return perm_;
  }

  void printInfo() const
  {
    std::cout << "\nCOO entries:\n";
    std::cout << "  row  col  value\n";
    for (int i = 0; i < nnz_; ++i)
    {
      std::cout << "  " << std::setw(3) << r_[i] << "  "
                << std::setw(3) << c_[i] << "  "
                << std::fixed << std::setprecision(1) << std::setw(5) << v_[i] << "\n";
    }
    std::cout << "\n";
  }

private:
  int m_;
  int n_;
  int nnz_dup_{0};
  int nnz_{0};

  std::vector<int>    r_;
  std::vector<int>    c_;
  std::vector<double> v_;

  std::vector<int> perm_;

  std::map<std::pair<int, int>, int> entries_;
};

class CsrMatrix
{
public:
  CsrMatrix(int m, int n, int nnz)
  {
    m_   = m;
    n_   = n;
    nnz_ = nnz;

    r_ptr_.resize(m_ + 1, 0);
    map_.resize(nnz_, 0);
    c_.resize(nnz_, 0);
    v_.resize(nnz_, 0);
  }

  ~CsrMatrix() = default;

  const std::vector<int>& get_row_ptr() const
  {
    return r_ptr_;
  }

  const std::vector<int>& get_column_ind() const
  {
    return c_;
  }

  const std::vector<double>& get_value() const
  {
    return v_;
  }

  const std::vector<int>& get_map2csr() const
  {
    return map_;
  }

  int get_nnz() const
  {
    return nnz_;
  }

  void compress(const std::vector<int>&    r,
                const std::vector<int>&    c,
                const std::vector<double>& v)
  {
    assert(static_cast<std::size_t>(nnz_) == r.size() && r.size() == c.size() && c.size() == v.size());

    std::vector<int> r_buffer(m_, 0);
    for (int i = 0; i < nnz_; ++i)
    {
      assert(r[i] >= 0 && r[i] < m_);
      assert(c[i] >= 0 && c[i] < n_);
      r_buffer[r[i]]++;
    }

    r_ptr_[0] = 0;
    for (int i = 0; i < m_; ++i)
    {
      r_ptr_[i + 1] = r_ptr_[i] + r_buffer[i];
    }

    assert(r_ptr_[m_] == nnz_);

    std::fill(r_buffer.begin(), r_buffer.end(), 0);

    for (int i = 0; i < nnz_; ++i)
    {
      map_[i]     = r_ptr_[r[i]] + r_buffer[r[i]];
      c_[map_[i]] = c[i];
      v_[map_[i]] = v[i];
      r_buffer[r[i]]++;
    }

    for (int i = 0; i < m_; ++i)
    {
      int r_start = r_ptr_[i];
      int r_end   = r_ptr_[i + 1];
      int r_size  = r_end - r_start;

      std::vector<std::pair<int, double>> r_data(r_size);
      for (int j = 0; j < r_size; ++j)
      {
        r_data[j] = {c_[r_start + j], v_[r_start + j]};
      }

      std::sort(r_data.begin(), r_data.end());

      for (int j = 0; j < r_size; ++j)
      {
        c_[r_start + j] = r_data[j].first;
        v_[r_start + j] = r_data[j].second;
      }
    }
  }

  // Update values in the original COO input order.
  void update(const std::vector<double>& vnew)
  {
    assert(vnew.size() == static_cast<std::size_t>(nnz_));
    for (int i = 0; i < nnz_; ++i)
    {
      v_[map_[i]] = vnew[i];
    }
  }

private:
  int m_;
  int n_;
  int nnz_;

  std::vector<int>    map_;
  std::vector<int>    r_ptr_;
  std::vector<int>    c_;
  std::vector<double> v_;
};

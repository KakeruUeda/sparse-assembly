#pragma once

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

/**
 * @brief CSR sparse matrix.
 *
 * Stores the matrix in CSR format and supports COO -> CSR conversion,
 * merging duplicate entries, and colum sorting.
 */
class CsrMatrix
{
public:
  /**
   * @brief Constructor
   *
   */
  CsrMatrix(int m, int n, int nnz)
  {
    m_      = m;
    n_      = n;
    nnzold_ = nnz;
    nnznew_ = nnz;

    rptr_.resize(m_ + 1, 0);
    map2csr_.resize(nnzold_, 0);
    cind_.resize(nnzold_, 0);
    vals_.resize(nnzold_, 0);
  }

  ~CsrMatrix() = default;

  const std::vector<int>& get_row_ptr() const
  {
    return rptr_;
  }

  const std::vector<int>& get_column_ind() const
  {
    return cind_;
  }

  const std::vector<double>& get_value() const
  {
    return vals_;
  }

  const std::vector<int>& get_map2csr() const
  {
    return map2csr_;
  }

  /**
   * @brief Convert COO to CSR.
   *
   * Builds CSR arrays from COO input. Duplicates are kept.
   *
   * @ref Timony A. Davis, Direct Methods for Sparse Linear Systems, SIAM.
   */
  void compress(const std::vector<int>&    r,
                const std::vector<int>&    c,
                const std::vector<double>& v)
  {
    assert(static_cast<std::size_t>(nnzold_) == r.size()
           && r.size() == c.size()
           && c.size() == v.size());

    // Count nnz per row
    for (int k = 0; k < nnzold_; ++k)
    {
      rptr_[r[k] + 1]++;
    }

    // Cumulative sum to get row pointers
    for (int i = 0; i < m_; ++i)
    {
      rptr_[i + 1] += rptr_[i];
    }
    assert(rptr_[m_] == nnzold_);

    // Next write position in each row
    std::vector<int> rpos(rptr_.begin(), rptr_.end() - 1);
    for (int n = 0; n < nnzold_; ++n)
    {
      int i = r[n];
      int p = rpos[i]++;

      map2csr_[n] = p;
      cind_[p]    = c[n];
      vals_[p]    = v[n];
    }
  }

  /**
   * @brief Remove duplicate entries.
   *
   * Combines duplicate entries by summing values.
   *
   * @ref Timony A. Davis, Direct Methods for Sparse Linear Systems, SIAM.
   *
   */
  void deduplicate()
  {
    // Mapping from old position to new position
    std::vector<int> map2new(nnzold_);

    // Compute deduplicated nnz_
    nnznew_ = 0;

    // Last written position per column
    std::vector<int> pos(n_, -1);
    for (int i = 0; i < m_; ++i)
    {
      int q = nnznew_;
      for (int p = rptr_[i]; p < rptr_[i + 1]; ++p)
      {
        int j = cind_[p];
        if (pos[j] >= q)
        {
          // Duplicate: map to existing position
          map2new[p] = pos[j];
          vals_[pos[j]] += vals_[p];
        }
        else
        {
          // Not duplicate: create new position
          map2new[p]     = nnznew_;
          pos[j]         = nnznew_;
          cind_[nnznew_] = j;
          vals_[nnznew_] = vals_[p];
          nnznew_++;
        }
      }
      rptr_[i] = q;
    }
    rptr_[m_] = nnznew_;

    // Update map2csr_ to point to new positions
    for (int n = 0; n < nnzold_; ++n)
    {
      map2csr_[n] = map2new[map2csr_[n]];
    }
  }

  /**
   * @brief Sort column indices within each row.
   *
   * Uses double transpose algorithm. O(m + n + nnz) complexity.
   *
   * @ref Timony A. Davis, Direct Methods for Sparse Linear Systems, SIAM.
   *
   */
  void sort()
  {
    // Build mapping from old position to new position
    std::vector<int> map2new(nnznew_);

    // First transpose
    // (m x n) -> (n x m)

    std::vector<int>    cptr(n_ + 1);
    std::vector<int>    rind(nnznew_);
    std::vector<double> vals(nnznew_);
    std::vector<int>    tmp(nnznew_); // intermidiate mapper

    // Count entries per column
    for (int p = 0; p < nnznew_; ++p)
    {
      cptr[cind_[p] + 1]++;
    }

    // Cumulative sum to get column pointers
    for (int j = 0; j < n_; ++j)
    {
      cptr[j + 1] += cptr[j];
    }
    assert(cptr[n_] == nnznew_);

    // Next write position in each column
    std::vector<int> cpos(cptr.begin(), cptr.end() - 1);

    // Transpose entries
    for (int i = 0; i < m_; ++i)
    {
      for (int p = rptr_[i]; p < rptr_[i + 1]; ++p)
      {
        int j   = cind_[p];
        int q   = cpos[j]++;
        rind[q] = i;
        vals[q] = vals_[p];
        tmp[q]  = p;
      }
    }

    // Second transpose
    // (n x m) -> (m x n)

    // Reset row pointers
    std::fill(rptr_.begin(), rptr_.end(), 0);

    // Count entries per row
    for (int p = 0; p < nnznew_; ++p)
    {
      rptr_[rind[p] + 1]++;
    }

    // Cumulative sum to get row pointers
    for (int i = 0; i < m_; ++i)
    {
      rptr_[i + 1] += rptr_[i];
    }
    assert(rptr_[m_] == nnznew_);

    // Next position in each row
    std::vector<int> rpos(rptr_.begin(), rptr_.end() - 1);

    // Place entries back (columns are in order 0..n-1)
    for (int j = 0; j < n_; ++j)
    {
      for (int p = cptr[j]; p < cptr[j + 1]; ++p)
      {
        int i           = rind[p];
        int q           = rpos[i]++;
        cind_[q]        = j;
        vals_[q]        = vals[p];
        map2new[tmp[p]] = q;
      }
    }

    // Update map2csr_ to point to new positions
    for (int n = 0; n < nnzold_; ++n)
    {
      map2csr_[n] = map2new[map2csr_[n]];
    }
  }

  /**
   * @brief Update matrix values.
   *
   * Uses mapping to update values without changing the sparsity pattern.
   */
  void update(const std::vector<double>& v)
  {
    assert(map2csr_.size() == v.size());

    // Reset values
    std::fill(vals_.begin(), vals_.end(), 0.0);

    // Update values with mapping
    for (int n = 0; n < nnzold_; ++n)
    {
      vals_[map2csr_[n]] += v[n];
    }
  }

  /**
   * @brief Print matrix entries.
   *
   */
  void printEntries() const
  {
    // Print row pointers
    std::cout << "Row pointers:   ";
    for (int i = 0; i <= m_; ++i)
    {
      std::cout << std::setw(6) << rptr_[i];
    }
    std::cout << "\n";

    // Print column indices
    std::cout << "Column indices: ";
    for (int i = 0; i < m_; ++i)
    {
      for (int p = rptr_[i]; p < rptr_[i + 1]; ++p)
      {
        std::cout << std::setw(6) << cind_[p];
      }
    }
    std::cout << "\n";

    // Print values
    std::cout << "Values:         ";
    for (int i = 0; i < m_; ++i)
    {
      for (int p = rptr_[i]; p < rptr_[i + 1]; ++p)
      {
        std::cout << std::setw(6) << std::fixed << std::setprecision(1) << vals_[p];
      }
    }
    std::cout << "\n";
  }

private:
  int m_;      // Number of rows
  int n_;      // Number of column
  int nnzold_; // Nonzeros before deduplication
  int nnznew_; // Nonzeros after deduplication

  std::vector<int>    map2csr_; // COO -> CSR mapper
  std::vector<int>    rptr_;    // CSR row pointer
  std::vector<int>    cind_;    // CSR column indices
  std::vector<double> vals_;    // CSR values
};

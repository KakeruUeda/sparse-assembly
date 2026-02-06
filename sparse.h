#pragma once

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <iostream>

/**
 * @brief COO sparse matrix.
 *
 * Stores the matrix in COO format and supports sorting and deduplication.
 * After processing, can provide row pointers for CSR conversion.
 */
class CooMatrix
{
public:
  /**
   * @brief Constructor
   *
   * @param m Number of rows
   * @param n Number of columns
   * @param nnz Number of nonzeros
   */
  CooMatrix(int m, int n, int nnz)
  {
    m_      = m;
    n_      = n;
    nnzdup_ = nnz;
    nnz_    = nnz;

    rind_ = new int[nnzdup_];
    cind_ = new int[nnzdup_];
    vals_ = new double[nnzdup_];
    ms_   = new int[nnzdup_];
  }

  ~CooMatrix()
  {
    delete[] rind_;
    delete[] cind_;
    delete[] vals_;
    delete[] rptr_;
    delete[] ms_;
    delete[] msnd_;
  }

  const int* getRowInd() const
  {
    return rind_;
  }

  const int* getColInd() const
  {
    return cind_;
  }

  const double* getValues() const
  {
    return vals_;
  }

  const int* getRowPtr() const
  {
    return rptr_;
  }

  const int* getMs() const
  {
    return ms_;
  }

  const int* getMsnd() const
  {
    return msnd_;
  }

  int getNnz() const
  {
    return nnz_;
  }

  int getNnzDup() const
  {
    return nnzdup_;
  }

  /**
   * @brief Add entries to COO matrix.
   *
   * @param r Row indices
   * @param c Column indices
   * @param v Values
   */
  void addEntries(const int* r, const int* c, const double* v)
  {
    for (int k = 0; k < nnzdup_; ++k)
    {
      rind_[k] = r[k];
      cind_[k] = c[k];
      vals_[k] = v[k];
    }
  }

  /**
   * @brief Sort entries by (row, column).
   *
   */
  void sort()
  {
    int* idx = new int[nnz_];
    for (int k = 0; k < nnz_; ++k)
    {
      idx[k] = k;
    }

    // Sort indices by (row, col)
    std::sort(idx, idx + nnz_, [this](int a, int b)
              { if (rind_[a] != rind_[b])
                { return rind_[a] < rind_[b]; }
                  return cind_[a] < cind_[b]; });

    // Apply permutation to arrays
    int*    rtmp = new int[nnz_];
    int*    ctmp = new int[nnz_];
    double* vtmp = new double[nnz_];

    for (int k = 0; k < nnz_; ++k)
    {
      rtmp[k] = rind_[idx[k]];
      ctmp[k] = cind_[idx[k]];
      vtmp[k] = vals_[idx[k]];
    }

    // Copy back
    for (int k = 0; k < nnz_; ++k)
    {
      rind_[k] = rtmp[k];
      cind_[k] = ctmp[k];
      vals_[k] = vtmp[k];
      ms_[k]   = idx[k];
    }

    delete[] idx;
    delete[] rtmp;
    delete[] ctmp;
    delete[] vtmp;
  }

  /**
   * @brief Remove duplicate entries by summing values.
   *
   */
  void deduplicate()
  {
    // Assumes entries are sorted.
    rptr_ = new int[m_ + 1];
    msnd_ = new int[nnzdup_];

    if (nnzdup_ == 0)
    {
      nnz_ = 0;
      return;
    }

    msnd_[0] = 0;
    rptr_[rind_[0] + 1]++;

    int w = 0;
    for (int r = 1; r < nnzdup_; ++r)
    {
      if (rind_[r] == rind_[w] && cind_[r] == cind_[w])
      {
        // Duplicate, sum values
        vals_[w] += vals_[r];
        msnd_[r] = w;
      }
      else
      {
        // New entry, advance write position
        w++;
        rind_[w] = rind_[r];
        cind_[w] = cind_[r];
        vals_[w] = vals_[r];
        msnd_[r] = w;

        rptr_[rind_[w] + 1]++;
      }
    }
    nnz_ = w + 1;

    // Cumulative sum
    for (int i = 0; i < m_; ++i)
    {
      rptr_[i + 1] += rptr_[i];
    }
    assert(rptr_[m_] == nnz_);
  }

  /**
   * @brief Print COO entries.
   */
  void printEntries() const
  {
    std::cout << "COO entries (nnz = " << nnz_ << "):\n";
    for (int k = 0; k < nnz_; ++k)
    {
      std::cout << "  (" << rind_[k] << ", " << cind_[k] << ") = "
                << std::fixed << std::setprecision(1) << vals_[k] << "\n";
    }
  }

private:
  int m_;      // Number of rows
  int n_;      // Number of columns
  int nnzdup_; // nnz (with duplicates)
  int nnz_;    // nnz (after deduplication)

  int*    rind_{nullptr}; // Row indices
  int*    cind_{nullptr}; // Column indices
  double* vals_{nullptr}; // Values
  int*    rptr_{nullptr}; // Row pointer

  int* ms_{nullptr};   // Mapping from sorted to original
  int* msnd_{nullptr}; // Mapping from deduplicated to sorted
};

/**
 * @brief CSR sparse matrix.
 *
 * Stores the matrix in CSR format from sorted/deduplicated COO entries.
 */
class CsrMatrix
{
public:
  /**
   * @brief Constructor
   *
   * @param m Number of rows
   * @param n Number of columns
   * @param nnz Number of nonzeros
   */
  CsrMatrix(int m, int n, int nnz)
  {
    m_   = m;
    n_   = n;
    nnz_ = nnz;

    rptr_ = new int[m_ + 1]{};
    cind_ = new int[nnz_]{};
    vals_ = new double[nnz_]{};
  }

  ~CsrMatrix()
  {
    delete[] rptr_;
    delete[] cind_;
    delete[] vals_;
  }

  const int* getRowPtr() const
  {
    return rptr_;
  }

  const int* getColInd() const
  {
    return cind_;
  }

  const double* getValues() const
  {
    return vals_;
  }

  int getNnz() const
  {
    return nnz_;
  }

  /**
   * @brief Build CSR from COO entries.
   *
   * @param rptr Row indices from COO
   * @param cind Column indices from COO
   * @param vals Values from COO
   */
  void addEntries(const int*    rptr,
                  const int*    cind,
                  const double* vals)
  {
    for (int i = 0; i < m_ + 1; ++i)
    {
      rptr_[i] = rptr[i];
    }
    for (int k = 0; k < nnz_; ++k)
    {
      cind_[k] = cind[k];
      vals_[k] = vals[k];
    }
  }

  /**
   * @brief Update values only (sparsity pattern unchanged).
   *
   * @param newvals New values (original COO order)
   * @param ms      Mapping: ms[sorted] = original
   * @param msnd    Mapping: msnd[sorted] = deduplicated
   * @param nnzdup  Original number of nonzeros
   */
  void updateValues(const double* newvals, const int* ms, const int* msnd, int nnzdup)
  {
    // Reset values
    for (int k = 0; k < nnz_; ++k)
    {
      vals_[k] = 0.0;
    }
    // Update values: iterate through sorted positions
    for (int k = 0; k < nnzdup; ++k)
    {
      vals_[msnd[k]] += newvals[ms[k]];
    }
  }

  /**
   * @brief Print CSR entries.
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
    for (int k = 0; k < nnz_; ++k)
    {
      std::cout << std::setw(6) << cind_[k];
    }
    std::cout << "\n";

    // Print values
    std::cout << "Values:         ";
    for (int k = 0; k < nnz_; ++k)
    {
      std::cout << std::setw(6) << std::fixed << std::setprecision(1) << vals_[k];
    }
    std::cout << "\n";
  }

private:
  int m_;   // Number of rows
  int n_;   // Number of columns
  int nnz_; // Number of nonzeros

  int*    rptr_{nullptr}; // Row pointers
  int*    cind_{nullptr}; // Column indices
  double* vals_{nullptr}; // Values
};

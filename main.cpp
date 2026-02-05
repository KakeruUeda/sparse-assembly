/**
 * @brief Demonstrates COO assembly from multiple components and convert to CSR format.
 * 
 */

#include "sparse.h"
#include "utils.h"

int main()
{
  // Local COO matrices

  // Component 0
  int    rind_c0[] = {0, 0, 1, 2};
  int    cind_c0[] = {0, 2, 1, 2};
  double vals_c0[] = {2.4, 0.6, 3.1, 5.0};
  int    nnz_c0    = 4;

  // Component 1
  int    rind_c1[] = {0, 1, 1};
  int    cind_c1[] = {1, 0, 1};
  double vals_c1[] = {0.5, 4.5, 1.1};
  int    nnz_c1    = 3;

  // Local index to global index
  int map_c0[] = {0, 1, 2}; // n0, n1, int
  int map_c1[] = {1, 0};    // n1, n0

  // Total number of nonzeros (duplicates included)
  int nnzdup = nnz_c0 + nnz_c1;

  // Create global COO
  int*    rind = new int[nnzdup];
  int*    cind = new int[nnzdup];
  double* vals = new double[nnzdup];

  // Insert components
  int counter = 0;
  for (int i = 0; i < nnz_c0; ++i)
  {
    rind[counter] = map_c0[rind_c0[i]];
    cind[counter] = map_c0[cind_c0[i]];
    vals[counter] = vals_c0[i];
    counter++;
  }
  for (int i = 0; i < nnz_c1; ++i)
  {
    rind[counter] = map_c1[rind_c1[i]];
    cind[counter] = map_c1[cind_c1[i]];
    vals[counter] = vals_c1[i];
    counter++;
  }

  int m = 3;
  int n = 3;

  int results = 0;

  // COO: add, sort, deduplicate
  CooMatrix coo(m, n, nnzdup);
  coo.addEntries(rind, cind, vals);

  // Sort
  coo.sort();
  results += testSort(coo);

  // Deduplicate
  coo.deduplicate();
  results += testDeduplicate(coo);

  int nnz = coo.getNnz();

  // Create CSR from sorted/deduplicated COO
  CsrMatrix csr(m, n, nnz);

  const int*    r = coo.getRowInd();
  const int*    c = coo.getColInd();
  const double* v = coo.getValues();

  csr.buildRowPtr(r, c, v);
  results += testCsr(csr);

  std::cout << "\nCSR entries:\n";
  csr.printEntries();
  std::cout << "\n";

  // Update values
  vals_c0[0] = 4.8;  // (0,0)
  vals_c0[1] = 1.2;  // (0,2)
  vals_c0[2] = 6.2;  // (1,1)
  vals_c0[3] = 10.0; // (2,2)
  vals_c1[0] = 1.0;  // (0,1)
  vals_c1[1] = 9.0;  // (1,0)
  vals_c1[2] = 2.2;  // (1,1)

  counter = 0;
  for (int i = 0; i < nnz_c0; ++i)
  {
    vals[counter] = vals_c0[i];
    counter++;
  }
  for (int i = 0; i < nnz_c1; ++i)
  {
    vals[counter] = vals_c1[i];
    counter++;
  }

  const int* ms   = coo.getMs();
  const int* msnd = coo.getMsnd();

  csr.updateValues(vals, ms, msnd, nnzdup);
  results += testUpdate(csr);

  std::cout << "\nCSR entries (after values updated):\n";
  csr.printEntries();

  delete[] rind;
  delete[] cind;
  delete[] vals;

  if (results == 0)
  {
    std::cout << "\nAll tests PASSED\n";
    return 0;
  }
  else
  {
    std::cout << "\nSome tests FAILED\n";
    return 1;
  }
}

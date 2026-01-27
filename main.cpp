/**
 * @brief Demonstrates COO assembly from multiple components
 *        and deduplication of matrix entries.
 * @todo  Convert COO to CSR, sort column indices,
 *        and preserve the COO-to-CSR mapping.
 *
 */

#include "sparse.h"
#include "utils.h"

/* ---------------------------
  Model:
    component0 -- node0 -- component1 -- node1

  n0: node0
  n1: node1
  int: internal

  r: row index
  c: column index
  p: permutation index
  v: value

  Local matrix (c0):
          n0   int
    n0  | 1.0  0.0 |
    int | 0.0  2.0 |

  Local matrix (c1):
          n0   n1   int
    n0  | 3.0  0.0  5.0 |
    n1  | 0.0  3.0  0.0 |

  Global matrix:
          n0   n1   int
    n0  | 4.0  0.0  5.0 |
    n1  | 0.0  3.0  0.0 |
    int | 0.0  0.0  2.0 |

  Raw entries before deduplication:
    r c p v
    0 0 0 1.0   // c0  (duplicate (0,0))
    2 2 1 2.0   // c0
    0 0 2 3.0   // c1  (duplicate (0,0))
    0 2 3 5.0   // c1
    1 1 4 3.0   // c1

  After deduplication:
    r c p v
    0 0 0 4.0   // 1.0 + 3.0
    2 2 1 2.0
    0 2 3 5.0
    1 1 4 3.0
--------------------------- */

int main()
{
  // Local COO (using global indices)

  // component0
  std::vector<int>    r_c0 = {0, 2};
  std::vector<int>    c_c0 = {0, 1};
  std::vector<double> v_c0 = {1.0, 2.0};

  // component1
  std::vector<int>    r_c1 = {0, 0, 1};
  std::vector<int>    c_c1 = {0, 2, 1};
  std::vector<double> v_c1 = {3.0, 5.0, 3.0};

  int m = 3;
  int n = 3;

  // Create COO matrix
  CooMatrix coo(m, n);
  coo.addEntries(r_c0, c_c0, v_c0);
  coo.addEntries(r_c1, c_c1, v_c1);

  // Print COO matrix entries
  coo.printInfo();

  // Expected results
  std::vector<int>    p_true  = {0, 1, 0, 2, 3};      // size 5 (= nnz(c0) + nnz(c1))
  std::vector<int>    r_true  = {0, 2, 0, 1};         // size 4 (= nnz(deduplicated))
  std::vector<int>    c_true  = {0, 1, 2, 1};         // size 4 (= nnz(deduplicated))
  std::vector<double> v1_true = {4.0, 2.0, 5.0, 3.0}; // size 4 (= nnz(deduplicated))

  const auto& r  = coo.getRows();
  const auto& c  = coo.getColumns();
  const auto& v1 = coo.getValues();
  const auto& p  = coo.getPermutation();

  // Validate results against expected results
  int failures = checkResults(r, c, v1, p, r_true, c_true, v1_true, p_true);

  // Update values
  std::vector<double> vnew = {
      7.0,  // c0 (0,0) updated: 1.0 -> 7.0
      2.0,  // c0 (2,1)
      2.0,  // c1 (0,0) updated: 3.0 -> 2.0
      10.0, // c1 (0,2) updated: 5.0 -> 10.0
      3.0   // c1 (1,1)
  };
  coo.update(vnew);

  std::vector<double> v2_true = {9.0, 2.0, 10.0, 3.0};
  const auto&         v2      = coo.getValues();

  failures += checkResults(r, c, v2, p, r_true, c_true, v2_true, p_true);

  if (failures == 0)
  {
    std::cout << "Test success" << "\n";
  }
  else
  {
    std::cout << "Test failed" << "\n";
  }

  return failures != 0;
}

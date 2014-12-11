/// Causal Dynamical Triangulations in C++ using CGAL
///
/// Copyright (c) 2014 Adam Getchell
///
/// Tests for the S3 bulk action

#include <vector>

#include "gmock/gmock.h"
#include "S3Triangulation.h"

using namespace testing;

class S3BulkAction : public Test {
 protected:
  virtual void SetUp() {
    make_S3_triangulation(&T, number_of_simplices,
                          number_of_timeslices, no_output,
                          &three_one, &two_two, &one_three);
  }

  const bool output = true;
  const bool no_output = false;
  const int number_of_simplices = 2;
  const int number_of_timeslices = 2;
  Delaunay T;
  std::vector<Cell_handle> three_one;
  std::vector<Cell_handle> two_two;
  std::vector<Cell_handle> one_three;
};

TEST_F(S3BulkAction, GetN3Values) {
  int N3_31 = three_one.size() + one_three.size();
  int N3_22 = two_two.size();
  std::cout << "N3(3,1) = " << N3_31 << std::endl;
  std::cout << "N3(2,2) = " << N3_22 << std::endl;
  ASSERT_EQ(T.number_of_finite_cells(), N3_31 + N3_22)
    << "N3(3,1) + N3(2,2) should be total number of cells.";
}

TEST_F(S3BulkAction, GetN1Values) {
  int N1_SL{0};
  int N1_TL{0};

  // classify_edges(&T, &N1_TL, &N1_SL);
  std::pair<unsigned, unsigned> Links = classify_edges(&T);
  N1_SL = Links.first;
  N1_TL = Links.second;

  ASSERT_EQ(T.number_of_finite_edges(), N1_TL + N1_SL)
    << "N1_TL + N1_SL should be total number of edges.";
}

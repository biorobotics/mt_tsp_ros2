#include "mt_tsp_ros2/trajopt_through_convex_sets.h"
#include <limits>
#include <Clarabel>
#include <fstream>

MatrixXd openData(std::string fileToOpen)
{

	// the inspiration for creating this function was drawn from here (I did NOT copy and paste the code)
	// https://stackoverflow.com/questions/34247057/how-to-read-csv-file-and-assign-to-eigen-matrix

	// the input is the file: "fileToOpen.csv":
	// a,b,c
	// d,e,f
	// This function converts input file data into the Eigen matrix format



	// the matrix entries are stored in this variable row-wise. For example if we have the matrix:
	// M=[a b c
	//	  d e f]
	// the entries are stored as matrixEntries=[a,b,c,d,e,f], that is the variable "matrixEntries" is a row vector
	// later on, this vector is mapped into the Eigen matrix format
	std::vector<double> matrixEntries;

	// in this object we store the data from the matrix
	std::ifstream matrixDataFile(fileToOpen);

	// this variable is used to store the row of the matrix that contains commas
	std::string matrixRowString;

	// this variable is used to store the matrix entry;
	std::string matrixEntry;

	// this variable is used to track the number of rows
	int matrixRowNumber = 0;


	while (getline(matrixDataFile, matrixRowString)) // here we read a row by row of matrixDataFile and store every line into the string variable matrixRowString
	{
		std::stringstream matrixRowStringStream(matrixRowString); //convert matrixRowString that is a string to a stream variable.

		while (getline(matrixRowStringStream, matrixEntry, ' ')) // here we read pieces of the stream matrixRowStringStream until every comma, and store the resulting character into the matrixEntry
		{
			matrixEntries.push_back(stod(matrixEntry));   //here we convert the string to double and fill in the row vector storing all the matrix entries
		}
		matrixRowNumber++; //update the column numbers
	}

	// here we convet the vector variable into the matrix and return the resulting object,
	// note that matrixEntries.data() is the pointer to the first memory location at which the entries of the vector matrixEntries are stored;
	return Map<Matrix<double, Dynamic, Dynamic, RowMajor>>(matrixEntries.data(), matrixRowNumber, matrixEntries.size() / matrixRowNumber);
}

TrajoptThroughConvexSets::TrajoptThroughConvexSets(int dim_q, const Ref<const VectorXd> &set_lb, const Ref<const VectorXd> &set_ub, const Ref<const VectorXd> &start_pos_and_time, double vmax_agent) : dim_q(dim_q), set_lb(set_lb), set_ub(set_ub), start_pos_and_time(start_pos_and_time) {
  set_dim = set_lb.size();

  z_idx = 0;
  l_idx = set_dim;

  vars_per_step = set_dim + dim_q + 1; // Line segment endpoints and displacement/distance aux vars
  steps = 0;

  zero_cone_rows_per_step = VectorXi::Zero(dim_q*3);
  zero_cone_cols_per_step = VectorXi::Zero(dim_q*3);
  zero_cone_vals_per_step = VectorXd::Zero(dim_q*3);
  for (int q_idx = 0; q_idx < dim_q; ++q_idx) {
    int start_idx = 3*q_idx;
    zero_cone_rows_per_step(start_idx) = q_idx;
    zero_cone_rows_per_step(start_idx + 1) = q_idx;
    zero_cone_rows_per_step(start_idx + 2) = q_idx;

    zero_cone_cols_per_step(start_idx) = l_idx + q_idx;
    zero_cone_cols_per_step(start_idx + 1) = z_idx + dim_q + 1 + q_idx;
    zero_cone_cols_per_step(start_idx + 2) = q_idx;

    zero_cone_vals_per_step(start_idx) = 1.;
    zero_cone_vals_per_step(start_idx + 1) = -1.;
    zero_cone_vals_per_step(start_idx + 2) = 1.;
  }
  zero_cone_b_per_step = VectorXd::Zero(dim_q);
  num_zero_cone_per_step = dim_q;

  VectorXi dist_cons_cols(3);
  dist_cons_cols(0) = l_idx + dim_q;
  dist_cons_cols(1) = z_idx + dim_q + 1 + dim_q;
  dist_cons_cols(2) = z_idx + dim_q;
  VectorXd dist_cons_vals(3);
  dist_cons_vals(0) = 1.;
  dist_cons_vals(1) = -vmax_agent;
  dist_cons_vals(2) = vmax_agent;
  num_dist_cons_rows = 1;
  int dist_cons_nnz_per_row = 3;

  VectorXi dist_cons_rows = VectorXi::Zero(3);

  VectorXd dist_cons_b = VectorXd::Zero(num_dist_cons_rows);

  VectorXi var_ub_rows = num_dist_cons_rows*VectorXi::Ones(set_dim);
  VectorXi var_ub_cols = VectorXi::Zero(set_dim);
  VectorXd var_ub_vals = VectorXd::Ones(set_dim);
  for (int i = 1; i < set_dim; ++i) {
    var_ub_rows(i) += i;
    var_ub_cols(i) = i;
  }
  VectorXd var_ub_b = set_ub;

  VectorXi var_lb_rows = (num_dist_cons_rows + set_dim)*VectorXi::Ones(set_dim);
  VectorXi var_lb_cols = VectorXi::Zero(set_dim);
  VectorXd var_lb_vals = -VectorXd::Ones(set_dim);
  for (int i = 1; i < set_dim; ++i) {
    var_lb_rows(i) += i;
    var_lb_cols(i) = i;
  }
  VectorXd var_lb_b = -set_lb;

  linear_cone_rows_per_step = VectorXi(dist_cons_rows.size() + var_ub_rows.size() + var_lb_rows.size());
  linear_cone_rows_per_step.head(dist_cons_rows.size()) = dist_cons_rows;
  linear_cone_rows_per_step.segment(dist_cons_rows.size(), var_ub_rows.size()) = var_ub_rows;
  linear_cone_rows_per_step.tail(var_lb_rows.size()) = var_lb_rows;

  linear_cone_cols_per_step = VectorXi(dist_cons_cols.size() + var_ub_cols.size() + var_lb_cols.size());
  linear_cone_cols_per_step.head(dist_cons_cols.size()) = dist_cons_cols;
  linear_cone_cols_per_step.segment(dist_cons_cols.size(), var_ub_cols.size()) = var_ub_cols;
  linear_cone_cols_per_step.tail(var_lb_cols.size()) = var_lb_cols;

  linear_cone_vals_per_step = VectorXd(dist_cons_vals.size() + var_ub_vals.size() + var_lb_vals.size());
  linear_cone_vals_per_step.head(dist_cons_vals.size()) = dist_cons_vals;
  linear_cone_vals_per_step.segment(dist_cons_vals.size(), var_ub_vals.size()) = var_ub_vals;
  linear_cone_vals_per_step.tail(var_lb_vals.size()) = var_lb_vals;

  linear_cone_b_per_step = VectorXd(dist_cons_b.size() + var_ub_b.size() + var_lb_b.size());
  linear_cone_b_per_step.head(dist_cons_b.size()) = dist_cons_b;
  linear_cone_b_per_step.segment(dist_cons_b.size(), var_ub_b.size()) = var_ub_b;
  linear_cone_b_per_step.tail(var_lb_b.size()) = var_lb_b;
  num_linear_cone_per_step = num_dist_cons_rows + 2*set_dim;

  soc_rows_per_step = VectorXi::Zero(1 + dim_q);
  soc_rows_per_step(0) = 0;
  soc_cols_per_step = VectorXi::Zero(1 + dim_q);
  soc_cols_per_step(0) = l_idx + dim_q;
  for (int q_idx = 0; q_idx < dim_q; ++q_idx) {
    soc_rows_per_step(1 + q_idx) = 1 + q_idx;
    soc_cols_per_step(1 + q_idx) = l_idx + q_idx;
  }
  soc_vals_per_step = -VectorXd::Ones(1 + dim_q);
  soc_b_per_step = VectorXd::Zero(1 + dim_q);
  num_soc_per_step = 1 + dim_q;

  num_zero_cone = 0;
  num_linear_cone = 0;
  num_soc = 0;
}

void TrajoptThroughConvexSets::add_convex_set(const GCSNode &convex_set) {
  if (steps == 0) {
    zero_cone_rows.push_back(num_zero_cone*VectorXi::Ones(dim_q + 1));
    zero_cone_cols.push_back(vars_per_step*steps*VectorXi::Ones(dim_q + 1));
    for (int i = 1; i < dim_q + 1; ++i) {
      zero_cone_rows.back()(i) += i;
      zero_cone_cols.back()(i) += i;
    }
    zero_cone_vals.push_back(VectorXd::Ones(dim_q + 1));
    zero_cone_b.push_back(start_pos_and_time);
    num_zero_cone += dim_q + 1;
  }

  zero_cone_rows.push_back(num_zero_cone*VectorXi::Ones(zero_cone_rows_per_step.size()) + zero_cone_rows_per_step);
  zero_cone_cols.push_back(vars_per_step*steps*VectorXi::Ones(zero_cone_cols_per_step.size()) + zero_cone_cols_per_step);
  zero_cone_vals.push_back(zero_cone_vals_per_step);
  zero_cone_b.push_back(zero_cone_b_per_step);
  num_zero_cone += num_zero_cone_per_step;

  linear_cone_rows.push_back(num_linear_cone*VectorXi::Ones(linear_cone_rows_per_step.size()) + linear_cone_rows_per_step);
  linear_cone_cols.push_back(vars_per_step*steps*VectorXi::Ones(linear_cone_cols_per_step.size()) + linear_cone_cols_per_step);
  linear_cone_vals.push_back(linear_cone_vals_per_step);
  linear_cone_b.push_back(linear_cone_b_per_step);
  num_linear_cone += num_linear_cone_per_step;

  soc_rows.push_back(num_soc*VectorXi::Ones(soc_rows_per_step.size()) + soc_rows_per_step);
  soc_cols.push_back(vars_per_step*steps*VectorXi::Ones(soc_cols_per_step.size()) + soc_cols_per_step);
  soc_vals.push_back(soc_vals_per_step);
  soc_b.push_back(soc_b_per_step);
  num_soc += num_soc_per_step;
  soc_sizes.push_back(num_soc_per_step);

  if (steps != 0) {
    // Make line segment begin where previous line segment ends
    zero_cone_rows.push_back(num_zero_cone*VectorXi::Ones(dim_q + 1));
    for (int i = 1; i < dim_q + 1; ++i) {
      zero_cone_rows.back()(i) += i;
    }
    zero_cone_rows.push_back(num_zero_cone*VectorXi::Ones(dim_q + 1));
    for (int i = 1; i < dim_q + 1; ++i) {
      zero_cone_rows.back()(i) += i;
    }

    zero_cone_cols.push_back((vars_per_step*steps + z_idx)*VectorXi::Ones(dim_q + 1));
    for (int i = 1; i < dim_q + 1; ++i) {
      zero_cone_cols.back()(i) += i;
    }

    zero_cone_cols.push_back((vars_per_step*(steps - 1) + z_idx + dim_q + 1)*VectorXi::Ones(dim_q + 1));
    for (int i = 1; i < dim_q + 1; ++i) {
      zero_cone_cols.back()(i) += i;
    }

    zero_cone_vals.push_back(VectorXd::Ones(dim_q + 1));
    zero_cone_vals.push_back(-VectorXd::Ones(dim_q + 1));

    zero_cone_b.push_back(VectorXd::Zero(dim_q + 1));
    num_zero_cone += dim_q + 1;
  }

  if (!std::isnan(convex_set.A_eq(0, 0))) {
    int nnz = convex_set.A_eq.rows()*convex_set.A_eq.cols();
    zero_cone_rows.push_back(num_zero_cone*VectorXi::Ones(nnz));
    zero_cone_cols.push_back(vars_per_step*steps*VectorXi::Ones(nnz));
    zero_cone_vals.push_back(VectorXd::Zero(nnz));
    for (int row = 0; row < convex_set.A_eq.rows(); ++row) {
      for (int col = 0; col < convex_set.A_eq.cols(); ++col) {
        int flat_idx = row*convex_set.A_eq.cols() + col;
        zero_cone_rows.back()(flat_idx) += row;
        zero_cone_cols.back()(flat_idx) += col;
        zero_cone_vals.back()(flat_idx) = convex_set.A_eq(row, col);
      }
    }
    zero_cone_b.push_back(-convex_set.b_eq);
    num_zero_cone += convex_set.b_eq.rows();
  }

  if (!std::isnan(convex_set.A_ineq(0, 0))) {
    int nnz = convex_set.A_ineq.rows()*convex_set.A_ineq.cols();
    linear_cone_rows.push_back(num_linear_cone*VectorXi::Ones(nnz));
    linear_cone_cols.push_back(vars_per_step*steps*VectorXi::Ones(nnz));
    linear_cone_vals.push_back(VectorXd::Zero(nnz));
    for (int row = 0; row < convex_set.A_ineq.rows(); ++row) {
      for (int col = 0; col < convex_set.A_ineq.cols(); ++col) {
        int flat_idx = row*convex_set.A_ineq.cols() + col;
        linear_cone_rows.back()(flat_idx) += row;
        linear_cone_cols.back()(flat_idx) += col;
        linear_cone_vals.back()(flat_idx) = convex_set.A_ineq(row, col);
      }
    }
    linear_cone_b.push_back(-convex_set.b_ineq);
    num_linear_cone += convex_set.b_ineq.rows();
  }
  steps += 1;
}

double TrajoptThroughConvexSets::solve_trajopt(std::vector<VectorXd> &pos_seq, std::vector<double> &time_seq, bool min_time, bool verbose) {
  int num_decision_vars = vars_per_step*steps;

  int zero_cone_nnz_idx = 0;
  std::vector<Triplet<double>> A_triplets;
  for (int nnz_vec_idx = 0; nnz_vec_idx < zero_cone_rows.size(); ++nnz_vec_idx) {
    for (int nnz_idx = 0; nnz_idx < zero_cone_rows[nnz_vec_idx].size(); ++nnz_idx) {
      A_triplets.push_back(Triplet<double>(zero_cone_rows[nnz_vec_idx](nnz_idx),
                                           zero_cone_cols[nnz_vec_idx](nnz_idx),    
                                           zero_cone_vals[nnz_vec_idx](nnz_idx)));
    }
  }

  int linear_cone_nnz_idx = 0;
  for (int nnz_vec_idx = 0; nnz_vec_idx < linear_cone_rows.size(); ++nnz_vec_idx) {
    for (int nnz_idx = 0; nnz_idx < linear_cone_rows[nnz_vec_idx].size(); ++nnz_idx) {
      A_triplets.push_back(Triplet<double>(num_zero_cone + linear_cone_rows[nnz_vec_idx](nnz_idx),
                                           linear_cone_cols[nnz_vec_idx](nnz_idx),    
                                           linear_cone_vals[nnz_vec_idx](nnz_idx)));
    }
  }

  int soc_nnz_idx = 0;
  for (int nnz_vec_idx = 0; nnz_vec_idx < soc_rows.size(); ++nnz_vec_idx) {
    for (int nnz_idx = 0; nnz_idx < soc_rows[nnz_vec_idx].size(); ++nnz_idx) {
      A_triplets.push_back(Triplet<double>(num_zero_cone + num_linear_cone + soc_rows[nnz_vec_idx](nnz_idx),
                                           soc_cols[nnz_vec_idx](nnz_idx),    
                                           soc_vals[nnz_vec_idx](nnz_idx)));
    }
  }

  int total_row = 0;
  VectorXd b(num_zero_cone + num_linear_cone + num_soc);
  for (auto b_vec : zero_cone_b) {
    b.segment(total_row, b_vec.size()) = b_vec;
    total_row += b_vec.size();
  }

  for (auto b_vec : linear_cone_b) {
    b.segment(total_row, b_vec.size()) = b_vec;
    total_row += b_vec.size();
  }

  for (auto b_vec : soc_b) {
    b.segment(total_row, b_vec.size()) = b_vec;
    total_row += b_vec.size();
  }
  assert(total_row == num_zero_cone + num_linear_cone + num_soc);

  VectorXd gradient = VectorXd::Zero(num_decision_vars);
  if (min_time) {
    gradient[num_decision_vars - vars_per_step + dim_q + 1 + dim_q] = 1.;
    gradient[dim_q] = -1.;
  } else {
    // Minimize distance travelled
    for (int step = 0; step < steps; ++step) {
      gradient(step*vars_per_step + l_idx + dim_q) = 1.;
    }
  }

  SparseMatrix<double> P(num_decision_vars, num_decision_vars);

  SparseMatrix<double> A(total_row, num_decision_vars);
  A.setFromTriplets(A_triplets.begin(), A_triplets.end());
  A.makeCompressed();

  std::vector<clarabel::SupportedConeT<double>> cones
  {
      clarabel::ZeroConeT<double>(num_zero_cone), clarabel::NonnegativeConeT<double>(num_linear_cone), 
  };
  for (auto s : soc_sizes) {
    cones.push_back(clarabel::SecondOrderConeT(s));
  }

  // Settings
  clarabel::DefaultSettings<double> settings = clarabel::DefaultSettings<double>::default_settings();
  settings.verbose = verbose;

  // Build solver
  clarabel::DefaultSolver<double> solver(P, gradient, A, b, cones, settings);

  // Solve
  solver.solve();

  // Get solution
  clarabel::DefaultSolution<double> soln = solver.solution();

  if (std::isnan(soln.obj_val)) {
    return std::numeric_limits<double>::infinity();
  }

  pos_seq.push_back(soln.x.head(dim_q));
  time_seq.push_back(soln.x(dim_q));
  for (int step = 0; step < steps; ++step) {
    pos_seq.push_back(soln.x.segment(vars_per_step*step + dim_q + 1, dim_q));
    time_seq.push_back(soln.x(vars_per_step*step + dim_q + 1 + dim_q));
  }

  return soln.obj_val;
}

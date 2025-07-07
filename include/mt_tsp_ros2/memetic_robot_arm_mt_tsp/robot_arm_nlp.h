#pragma once

#include <Eigen/Dense>

#define HAVE_CSTDDEF
#include <IpTNLP.hpp>
#undef HAVE_CSTDDEF
#define HAVE_STDDEF
#include <IpTNLP.hpp>
#undef HAVE_STDDEF
#include "IpIpoptApplication.hpp"
#include "IpIpoptCalculatedQuantities.hpp"

#include "mt_tsp_ros2/SE3_spline.h"
#include "mt_tsp_ros2/memetic_robot_arm_mt_tsp/kuka_ik.h"

#include <chrono>

using namespace std::chrono;

using namespace Eigen;

typedef Matrix<Ipopt::Index, Dynamic, 1> VectorXIndex;
typedef Matrix<long, Dynamic, 1> VectorXl;
typedef Matrix<double, Dynamic, Dynamic, RowMajor> RowMatrixXd;

const double inf_val = 1e19;

const bool multiply_speed_constraint_by_delta_t = true;

const bool gauss_newton_hessian = false;

const bool quadratic_cost = false;

struct RobotArmNLPInfo {
  RobotArmNLPInfo(const Ref<const RowMatrixXd> &tw_per_target,
                  const std::vector<SE3Spline> &q_trj_per_target,
                  const Ref<const VectorXd> &q0,
                  const Ref<const VectorXd> &joint_limits,
                  const Ref<const VectorXd> &vmax,
                  const Ref<const Vector4d> &l,
                  const Ref<const RowMatrixXd> &dh) : tw_per_target(tw_per_target),
                                                      q_trj_per_target(q_trj_per_target),
                                                      q0(q0), 
                                                      joint_limits(joint_limits),
                                                      vmax(vmax),
                                                      l(l),
                                                      dh(dh) {
  }
  RowMatrixXd tw_per_target;
  std::vector<SE3Spline> q_trj_per_target;
  VectorXd q0;
  VectorXd joint_limits;
  VectorXd vmax;
  Vector4d l;
  RowMatrixXd dh;
};

class RobotArmNLP : public Ipopt::TNLP {
  public:
    /** Constructor */
    RobotArmNLP(std::shared_ptr<RobotArmNLPInfo> info,
                const Ref<const VectorXl> &target_seq) : info(info),
                                                         target_seq(target_seq) {
      finite_diff_gradient = false;
      num_targets = info->tw_per_target.rows();
      vars_per_step = 1 + dim_q; // t and q
      num_decision_vars = vars_per_step*num_targets;
      warm_start = VectorXd::Zero(0);
      int constraints_per_step;
      if (multiply_speed_constraint_by_delta_t) {
        constraints_per_step = 3 + 3 + dim_q + dim_q;
        J_nnz = 6*vars_per_step*num_targets + 2*2*dim_q + 2*4*dim_q*(num_targets - 1);
      } else {
        constraints_per_step = 3 + 3 + dim_q;
        J_nnz = 6*vars_per_step*num_targets + 2*dim_q + 4*dim_q*(num_targets - 1);
      }
      num_constraints = constraints_per_step*num_targets;
      // J_nnz = num_constraints*num_decision_vars; // Dense

      if (gauss_newton_hessian) {
        if (quadratic_cost) {
          H_nnz = dim_q + 2*dim_q*(num_targets - 1);
        } else {
          H_nnz = dim_q*(dim_q + 1)/2 + (dim_q*(dim_q + 1)/2 + dim_q*dim_q)*(num_targets - 1); // dim_q*(dim_q + 1)/2 is for step 0 and (dim_q*(dim_q + 1)/2 + dim_q*dim_q) is for subsequent steps
        }
      } else {
        H_nnz = vars_per_step*(vars_per_step + 1)/2 + (vars_per_step*(vars_per_step + 1)/2 + vars_per_step*vars_per_step)*(num_targets - 1); // vars_per_step*(vars_per_step + 1)/2 is for step 0 and (vars_per_step*(vars_per_step + 1)/2 + vars_per_step*vars_per_step) is for subsequent steps
      }
      // H_nnz = num_decision_vars*num_decision_vars; // Dense

      // Populate Jacobian sparsity
      Jrows = VectorXIndex(J_nnz);
      Jcols = VectorXIndex(J_nnz);
      
      int constraint_idx = 0;
      int nnz_idx = 0;
      for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
        // EE position error = 0
        for (int i = 0; i < 6; ++i) {
          for (int j = 0; j < vars_per_step; ++j) {
            Jrows(nnz_idx + i*vars_per_step + j) = constraint_idx + i; // ith element of position error, if i < 3, otherwise (i - 3)th element of rotation error
            Jcols(nnz_idx + i*vars_per_step + j) = vars_per_step*seq_idx + j; // t, or (j - 1)th element of q
          }
        }
        constraint_idx += 6;
        nnz_idx += 6*(1 + dim_q);

        // Joint velocity limits satisfied
        int nnz_per_q_component = seq_idx == 0 ? 2 : 4;
        for (int sign_idx = 0; sign_idx < (multiply_speed_constraint_by_delta_t ? 2 : 1); ++sign_idx) {
          for (int i = 0; i < dim_q; ++i) {
            Jrows(nnz_idx + nnz_per_q_component*i) = constraint_idx + i; // ith element of constraint
            Jcols(nnz_idx + nnz_per_q_component*i) = vars_per_step*seq_idx + 1 + i; // ith element of q

            Jrows(nnz_idx + nnz_per_q_component*i + 1) = constraint_idx + i; // ith element of constraint
            Jcols(nnz_idx + nnz_per_q_component*i + 1) = vars_per_step*seq_idx; // current t

            if (seq_idx == 0) {
              continue;
            }

            Jrows(nnz_idx + nnz_per_q_component*i + 2) = constraint_idx + i; // ith element of constraint
            Jcols(nnz_idx + nnz_per_q_component*i + 2) = vars_per_step*(seq_idx - 1) + 1 + i; // ith element of prev_q

            Jrows(nnz_idx + nnz_per_q_component*i + 3) = constraint_idx + i; // ith element of constraint
            Jcols(nnz_idx + nnz_per_q_component*i + 3) = vars_per_step*(seq_idx - 1); // prev_t
          }
          constraint_idx += dim_q;
          nnz_idx += nnz_per_q_component*dim_q;
        }
      }

      if (nnz_idx != J_nnz) {
        throw std::runtime_error("Number of added nonzeros not equal to J_nnz");
      }

      // Dense
      /*
      int nnz_idx = 0;
      for (int row = 0; row < num_constraints; ++row) {
        for (int col = 0; col < num_decision_vars; ++col) {
          Jrows(row*num_decision_vars + col) = row;
          Jcols(row*num_decision_vars + col) = col;
        }
      }
      */

      // Populate Hessian sparsity
      Hrows = VectorXIndex(H_nnz);
      Hcols = VectorXIndex(H_nnz);

      if (gauss_newton_hessian) {
        if (quadratic_cost) {
          nnz_idx = 0;
          for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
            for (int i = 1; i < vars_per_step; ++i) {
              Hrows(nnz_idx) = seq_idx*vars_per_step + i;
              Hcols(nnz_idx) = seq_idx*vars_per_step + i;
              ++nnz_idx;
            }

            if (seq_idx == 0) {
              continue;
            }

            for (int i = 1; i < vars_per_step; ++i) {
              Hrows(nnz_idx) = seq_idx*vars_per_step + i;
              Hcols(nnz_idx) = (seq_idx - 1)*vars_per_step + i;
              ++nnz_idx;
            }
          }
        } else {
          nnz_idx = 0;
          for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
            for (int i = 1; i < vars_per_step; ++i) {
              for (int j = 1; j <= i; ++j) { // Only filling lower triangle
                Hrows(nnz_idx) = seq_idx*vars_per_step + i;
                Hcols(nnz_idx) = seq_idx*vars_per_step + j;
                ++nnz_idx;
              }
            }

            if (seq_idx == 0) {
              continue;
            }

            for (int i = 1; i < vars_per_step; ++i) {
              for (int j = 1; j < vars_per_step; ++j) {
                Hrows(nnz_idx) = seq_idx*vars_per_step + i;
                Hcols(nnz_idx) = (seq_idx - 1)*vars_per_step + j;
                ++nnz_idx;
              }
            }
          }
        }
      } else {
        nnz_idx = 0;
        for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
          for (int i = 0; i < vars_per_step; ++i) {
            for (int j = 0; j <= i; ++j) { // Only filling lower triangle
              Hrows(nnz_idx) = seq_idx*vars_per_step + i;
              Hcols(nnz_idx) = seq_idx*vars_per_step + j;
              ++nnz_idx;
            }
          }

          if (seq_idx == 0) {
            continue;
          }

          for (int i = 0; i < vars_per_step; ++i) {
            for (int j = 0; j < vars_per_step; ++j) {
              Hrows(nnz_idx) = seq_idx*vars_per_step + i;
              Hcols(nnz_idx) = (seq_idx - 1)*vars_per_step + j;
              ++nnz_idx;
            }
          }
        }
      }
      if (nnz_idx != H_nnz) {
        throw std::runtime_error("Number of added nonzeros not equal to H_nnz");
      }

      // Dense
      /*
      for (int row = 0; row < num_decision_vars; ++row) {
        for (int col = 0; col < num_decision_vars; ++col) {
          Hrows(row*num_decision_vars + col) = row;
          Hcols(row*num_decision_vars + col) = col;
        }
      }
      */

      x_l = VectorXd::Zero(num_decision_vars);
      x_u = VectorXd::Zero(num_decision_vars);

      g_l = VectorXd::Zero(num_constraints);
      g_u = VectorXd::Zero(num_constraints);

      constraint_idx = 0;
      for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
        int target_idx = target_seq(seq_idx);
        // Time window
        x_l(vars_per_step*seq_idx) = info->tw_per_target(target_idx, 0);
        x_u(vars_per_step*seq_idx) = info->tw_per_target(target_idx, 1);

        // Joint limits
        x_l.segment(vars_per_step*seq_idx + 1, dim_q) = -info->joint_limits;
        x_u.segment(vars_per_step*seq_idx + 1, dim_q) = info->joint_limits;
        /*
        x_l.segment(vars_per_step*seq_idx + 1, dim_q).setConstant(-inf_val);
        x_u.segment(vars_per_step*seq_idx + 1, dim_q).setConstant(inf_val);
        */

        // EE position error = 0
        g_l.segment(constraint_idx, 3).setZero();
        g_u.segment(constraint_idx, 3).setZero();
        /*
        g_l.segment(constraint_idx, 3).setConstant(-inf_val);
        g_u.segment(constraint_idx, 3).setConstant(inf_val);
        */
        constraint_idx += 3;

        // EE rotation error = 0
        g_l.segment(constraint_idx, 3).setZero();
        g_u.segment(constraint_idx, 3).setZero();
        /*
        g_l.segment(constraint_idx, 3).setConstant(-inf_val);
        g_u.segment(constraint_idx, 3).setConstant(inf_val);
        */
        constraint_idx += 3;

        // Joint velocity limits
        if (multiply_speed_constraint_by_delta_t) {
          // vmax*delta_t - delta_q >= 0
          g_l.segment(constraint_idx, dim_q).setZero();
          g_u.segment(constraint_idx, dim_q).setConstant(inf_val);
          constraint_idx += dim_q;

          // delta_q - -vmax*delta_t >= 0
          g_l.segment(constraint_idx, dim_q).setZero();
          g_u.segment(constraint_idx, dim_q).setConstant(inf_val);
          constraint_idx += dim_q;
        } else {
          // -vmax <= delta_q/delta_t <= vmax
          g_l.segment(constraint_idx, dim_q) = -info->vmax;
          g_u.segment(constraint_idx, dim_q) = info->vmax;
          /*
          g_l.segment(constraint_idx, dim_q).setConstant(-inf_val);
          g_u.segment(constraint_idx, dim_q).setConstant(inf_val);
          */
          constraint_idx += dim_q;
        }
      }
    }

    /**@name Overloaded from TNLP */
    //@{
    /** Method to return some info about the NLP */
    virtual bool get_nlp_info(
       Ipopt::Index&          n,
       Ipopt::Index&          m,
       Ipopt::Index&          nnz_jac_g,
       Ipopt::Index&          nnz_h_lag,
       Ipopt::TNLP::IndexStyleEnum& index_style
    ) {
      n = num_decision_vars;
      m = num_constraints;
      nnz_jac_g = J_nnz;
      nnz_h_lag = H_nnz;
      index_style = TNLP::C_STYLE;
      return true;
    }

    /** Method to return the bounds for my problem */
    virtual bool get_bounds_info(
       Ipopt::Index   n,
       Ipopt::Number* x_l,
       Ipopt::Number* x_u,
       Ipopt::Index   m,
       Ipopt::Number* g_l,
       Ipopt::Number* g_u
    ) {
      Map<VectorXd>(x_l, num_decision_vars) = this->x_l;
      Map<VectorXd>(x_u, num_decision_vars) = this->x_u;

      Map<VectorXd>(g_l, num_constraints) = this->g_l;
      Map<VectorXd>(g_u, num_constraints) = this->g_u;
      return true;
    }

    /** Method to return the starting point for the algorithm */
    virtual bool get_starting_point(
       Ipopt::Index   n,
       bool    init_x,
       Ipopt::Number* x,
       bool    init_z,
       Ipopt::Number* z_L,
       Ipopt::Number* z_U,
       Ipopt::Index   m,
       bool    init_lambda,
       Ipopt::Number* lambda
    ) {
      if (init_x) {
        if (warm_start.size() == num_decision_vars) {
          Map<VectorXd>(x, num_decision_vars) = warm_start;
        } else {
          // Initialize qs at q0 and ts at center of time window
          for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
            int target_idx = target_seq(seq_idx);
            // Time window
            x[vars_per_step*seq_idx] = 0.5*(info->tw_per_target(target_idx, 0) + info->tw_per_target(target_idx, 1));

            // Joint limits
            Map<VectorXd>(x + vars_per_step*seq_idx + 1, dim_q) = info->q0;
          }
        }
      }
      if (init_z) {
        Map<VectorXd>(z_L, n).setZero();
        Map<VectorXd>(z_U, n).setZero();
      }
      if (init_lambda) {
        Map<VectorXd>(lambda, m).setZero();
      }
      return true;
    }

    /** Method to return the objective value */
    virtual bool eval_f(
       Ipopt::Index         n,
       const Ipopt::Number* x,
       bool          new_x,
       Ipopt::Number&       obj_value
    ) {
      obj_value = 0.;
      VectorXd prev_q = info->q0;
      for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
        Map<const VectorXd> q(x + seq_idx*vars_per_step + 1, dim_q);
        if (quadratic_cost) {
          obj_value += (q - prev_q).squaredNorm();
        } else {
          obj_value += sqrt((q - prev_q).squaredNorm() + 1e-8);
        }
        prev_q = q;
      }
      return true;
    }

    /** Method to return the gradient of the objective */
    virtual bool eval_grad_f(
       Ipopt::Index         n,
       const Ipopt::Number* x,
       bool          new_x,
       Ipopt::Number*       grad_f
    ) {
      if (finite_diff_gradient) {
        double eps = 1e-4;
        double o;
        eval_f(n, x, new_x, o);
        VectorXd xplus(Map<const VectorXd>(x, num_decision_vars));
        for (int i = 0; i < num_decision_vars; ++i) {
          xplus(i) += eps;
          double oplus;
          eval_f(n, xplus.data(), new_x, oplus);
          xplus(i) = x[i];
          grad_f[i] = (oplus - o)/eps;
        }
      } else {
        // Map<VectorXd>(grad_f, num_decision_vars).setZero();
        VectorXd prev_q = info->q0;
        for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
          Map<const VectorXd> q(x + seq_idx*vars_per_step + 1, dim_q);
          double dist = sqrt((q - prev_q).squaredNorm() + 1e-8);
          if (quadratic_cost) {
            Map<VectorXd>(grad_f + seq_idx*vars_per_step + 1, dim_q) = 2*(q - prev_q);
          } else {
            Map<VectorXd>(grad_f + seq_idx*vars_per_step + 1, dim_q) = 1/dist*(q - prev_q);
          }
          grad_f[seq_idx*vars_per_step] = 0.;
          if (seq_idx != 0) {
            if (quadratic_cost) {
              Map<VectorXd>(grad_f + (seq_idx - 1)*vars_per_step + 1, dim_q) += -2*(q - prev_q);
            } else {
              Map<VectorXd>(grad_f + (seq_idx - 1)*vars_per_step + 1, dim_q) += -1/dist*(q - prev_q);
            }
          }

          prev_q = q;
        }
      }
      return true;
    }

    /** Method to return the constraint residuals */
    virtual bool eval_g(
       Ipopt::Index         n,
       const Ipopt::Number* x,
       bool          new_x,
       Ipopt::Index         m,
       Ipopt::Number*       g
    ) {
      int constraint_idx = 0;
      VectorXd prev_q = info->q0;
      double prev_t = 0.;
      for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
        int target_idx = target_seq(seq_idx);
        double t = x[vars_per_step*seq_idx];
        Map<const VectorXd> q(x + seq_idx*vars_per_step + 1, dim_q);
        Matrix4d ee_pose = kuka_fk(q, info->l, info->dh);
        Matrix4d target_pose = info->q_trj_per_target[target_idx](t);

        // EE position error = 0
        Vector3d position_error = ee_pose.topRightCorner<3, 1>() - target_pose.topRightCorner<3, 1>();
        Map<Vector3d>(g + constraint_idx) = position_error;
        constraint_idx += 3;

        /*
        if (position_error.maxCoeff() > 1e-4) {
          std::cout << "ee position constraint violated" << std::endl;
          throw std::runtime_error("error");
        }

        if (position_error.minCoeff() < -1e-4) {
          std::cout << "ee position constraint violated" << std::endl;
          throw std::runtime_error("error");
        }
        */

        // EE rotation error = 0
        Matrix3d rotation_error_matrix = target_pose.topLeftCorner<3, 3>().transpose()*ee_pose.topLeftCorner<3, 3>();
        Matrix3d skew = 0.5*(rotation_error_matrix - rotation_error_matrix.transpose());
        Vector3d rotation_error(-skew(0, 1), skew(0, 2), -skew(1, 2));
        Map<Vector3d>(g + constraint_idx) = rotation_error;
        constraint_idx += 3;

        /*
        if (rotation_error.maxCoeff() > 1e-4) {
          std::cout << "ee rotation constraint violated" << std::endl;
          throw std::runtime_error("error");
        }

        if (rotation_error.minCoeff() < -1e-4) {
          std::cout << "ee rotation constraint violated" << std::endl;
          throw std::runtime_error("error");
        }
        */

        // Joint velocity limits
        if (multiply_speed_constraint_by_delta_t) {
          double delta_t = t - prev_t;
          Map<VectorXd>(g + constraint_idx, dim_q) = info->vmax*delta_t - (q - prev_q);
          constraint_idx += dim_q;

          Map<VectorXd>(g + constraint_idx, dim_q) = (q - prev_q) - -info->vmax*delta_t;
          constraint_idx += dim_q;
        } else {
          double delta_t = t - prev_t;
          Map<VectorXd>(g + constraint_idx, dim_q) = (q - prev_q)/delta_t;
          constraint_idx += dim_q;
        }

        prev_t = t;
        prev_q = q;
      }
      if (constraint_idx < num_constraints) {
        std::cout << "Did not populate all constraints" << std::endl;
        throw std::runtime_error("Did not populate all constraints");
      }
      return true;
    }

    /** Method to return:
     *   1) The structure of the jacobian (if "values" is NULL)
     *   2) The values of the jacobian (if "values" is not NULL)
     */
    virtual bool eval_jac_g(
       Ipopt::Index         n,
       const Ipopt::Number* x,
       bool          new_x,
       Ipopt::Index         m,
       Ipopt::Index         nele_jac,
       Ipopt::Index*        iRow,
       Ipopt::Index*        jCol,
       Ipopt::Number*       values
    ) {
      if (values == NULL) {
        Map<VectorXIndex>(iRow, nele_jac) = Jrows;
        Map<VectorXIndex>(jCol, nele_jac) = Jcols;
      } else {
        std::cout << "Did not implement exact constraint Jacobian" << std::endl;
        throw std::runtime_error("Did not implement exact constraint Jacobian");
      }
      return true;
    }

    /** Method to return:
     *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
     *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
     */
    virtual bool eval_h(
       Ipopt::Index         n,
       const Ipopt::Number* x,
       bool          new_x,
       Ipopt::Number        obj_factor,
       Ipopt::Index         m,
       const Ipopt::Number* lambda,
       bool          new_lambda,
       Ipopt::Index         nele_hess,
       Ipopt::Index*        iRow,
       Ipopt::Index*        jCol,
       Ipopt::Number*       values
    ) {
      if (values == NULL) {
        Map<VectorXIndex>(iRow, nele_hess) = Hrows;
        Map<VectorXIndex>(jCol, nele_hess) = Hcols;
      } else {
        if (!gauss_newton_hessian) {
          std::cout << "Did not implement exact constraint Hessian for non-Gauss-Newton case" << std::endl;
          throw std::runtime_error("Did not implement exact constraint Hessian for non-Gauss-Newton case");
        }

        // Gauss-Newton with quadratic cost
        if (quadratic_cost) {
          int nnz_idx = 0;
          VectorXd prev_q = info->q0;
          for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
            if (seq_idx == 0) {
              for (int i = 1; i < vars_per_step; ++i) {
                values[nnz_idx] = 4*obj_factor;
                ++nnz_idx;
              }
            } else if (seq_idx == num_targets - 1) {
              for (int i = 1; i < vars_per_step; ++i) {
                values[nnz_idx] = 2*obj_factor;
                ++nnz_idx;
              }
              for (int i = 1; i < vars_per_step; ++i) {
                values[nnz_idx] = -2*obj_factor;
                ++nnz_idx;
              }
            } else {
              for (int i = 1; i < vars_per_step; ++i) {
                values[nnz_idx] = 4*obj_factor;
                ++nnz_idx;
              }
              for (int i = 1; i < vars_per_step; ++i) {
                values[nnz_idx] = -2*obj_factor;
                ++nnz_idx;
              }
            }
          }
        } else {
          std::cout << "Did not implement Gauss-Newton hessian for distance cost" << std::endl;
          throw std::runtime_error("Did not implement Gauss-Newton hessian for distance cost");
        }
      }
      return true;
    }

    /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
    virtual void finalize_solution(
       Ipopt::SolverReturn               status,
       Ipopt::Index                      n,
       const Ipopt::Number*              x,
       const Ipopt::Number*              z_L,
       const Ipopt::Number*              z_U,
       Ipopt::Index                      m,
       const Ipopt::Number*              g,
       const Ipopt::Number*              lambda,
       Ipopt::Number                     obj_value,
       const Ipopt::IpoptData*           ip_data,
       Ipopt::IpoptCalculatedQuantities* ip_cq
    ) {
      soln = Map<const VectorXd>(x, n);
      cost = obj_value;
    }

    bool intermediate_callback(
       Ipopt::AlgorithmMode              mode,
       Ipopt::Index                      iter,
       Ipopt::Number                     obj_value,
       Ipopt::Number                     inf_pr,
       Ipopt::Number                     inf_du,
       Ipopt::Number                     mu,
       Ipopt::Number                     d_norm,
       Ipopt::Number                     regularization_size,
       Ipopt::Number                     alpha_du,
       Ipopt::Number                     alpha_pr,
       Ipopt::Index                      ls_trials,
       const Ipopt::IpoptData*           ip_data,
       Ipopt::IpoptCalculatedQuantities* ip_cq
    ) {
      if (mode == Ipopt::AlgorithmMode::RestorationPhaseMode) {
        // throw std::runtime_error("In restoration mode");
        restoration_invoked = true;
      }
      return true;
    }

    const VectorXd &get_soln() {
      return soln;
    }

    double get_cost() {
      return cost;
    }

    const VectorXIndex &get_Jrows() {
      return Jrows;
    }

    const VectorXIndex &get_Jcols() {
      return Jcols;
    }

    const VectorXIndex &get_Hrows() {
      return Hrows;
    }

    const VectorXIndex &get_Hcols() {
      return Hcols;
    }
    
    int get_num_constraints() {
      return num_constraints;
    }

    int get_num_decision_vars() {
      return num_decision_vars;
    }

    bool get_restoration_invoked() {
      return restoration_invoked;
    }

    void set_warm_start(const Ref<const VectorXd> &warm_start) {
      this->warm_start = warm_start;
      // this->x_l = warm_start;
      // this->x_u = warm_start;
    }

    void set_finite_diff_gradient(bool finite_diff_gradient) {
      this->finite_diff_gradient = finite_diff_gradient;
    }

  private:
    std::shared_ptr<RobotArmNLPInfo> info;
    VectorXl target_seq;

    int num_targets;

    int vars_per_step;
    int num_decision_vars;
    int J_nnz;
    int H_nnz;
    int num_constraints;

    VectorXIndex Jrows;
    VectorXIndex Jcols;

    VectorXIndex Hrows;
    VectorXIndex Hcols;

    VectorXd x_l;
    VectorXd x_u;

    VectorXd g_l;
    VectorXd g_u;

    VectorXd soln;

    double cost;

    VectorXd warm_start;

    bool finite_diff_gradient;

    bool restoration_invoked;
};

class RobotArmNLPSolver {
  public:
    void initialize(std::shared_ptr<RobotArmNLPInfo> info, const Ref<const VectorXl> &target_seq, int max_iter, int print_level) {
      num_targets = target_seq.size();
      nlp = new RobotArmNLP(info, target_seq);

      // Create a new instance of IpoptApplication
      //  (use a SmartPtr, not raw)
      // We are using the factory, since this allows us to compile this
      // example with an Ipopt Windows DLL
      app = new Ipopt::IpoptApplication();
      
      // Change some options
      // Note: The following choices are only examples, they might not be
      //       suitable for your optimization problem.
      app->Options()->SetNumericValue("tol", 1.0);
      app->Options()->SetNumericValue("compl_inf_tol", 1e-2);
      /*
      app->Options()->SetNumericValue("constr_vio_tol", 1e-4);
      */
      app->Options()->SetIntegerValue("max_iter", max_iter);
      app->Options()->SetStringValue("jacobian_approximation", "finite-difference-values");
      // app->Options()->SetStringValue("gradient_approximation", "finite-difference-values");
      if (!gauss_newton_hessian) {
        app->Options()->SetStringValue("hessian_approximation", "limited-memory");
      }

      app->Options()->SetIntegerValue("print_level", print_level);

      // Initialize the IpoptApplication and process the options
      Ipopt::ApplicationReturnStatus status;
      status = app->Initialize();
      if( status != Ipopt::Solve_Succeeded )
      {
         std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
         std::cout << "Status " << (int) status << std::endl;
         throw std::runtime_error("Ipopt initialization error");
      }
    }

    RobotArmNLPSolver(std::shared_ptr<RobotArmNLPInfo> info, const Ref<const VectorXl> &target_seq, int max_iter, int print_level) {
      initialize(info, target_seq, max_iter, print_level);
    }

    RobotArmNLPSolver(const Ref<const RowMatrixXd> &tw_per_target, const std::vector<SE3Spline> &q_trj_per_target, const Ref<const VectorXd> &q0, const Ref<const VectorXl> &target_seq, int max_iter, const Ref<const VectorXd> &joint_limits, const Ref<const VectorXd> &vmax, const Ref<const Vector4d> &l, const Ref<const RowMatrixXd> &dh, int print_level) {
      std::shared_ptr<RobotArmNLPInfo> info = std::make_shared<RobotArmNLPInfo>(tw_per_target, q_trj_per_target, q0, joint_limits, vmax, l, dh);
      initialize(info, target_seq, max_iter, print_level);
    }

    void set_finite_diff_gradient(bool finite_diff_gradient) {
      nlp->set_finite_diff_gradient(finite_diff_gradient);
    }

    VectorXIndex get_Jrows() {
      return nlp->get_Jrows();
    }

    VectorXIndex get_Jcols() {
      return nlp->get_Jcols();
    }

    VectorXIndex get_Hrows() {
      return nlp->get_Hrows();
    }

    VectorXIndex get_Hcols() {
      return nlp->get_Hcols();
    }

    int get_num_constraints() {
      return nlp->get_num_constraints();
    }

    int get_num_decision_vars() {
      return nlp->get_num_decision_vars();
    }

    void set_warm_start(const Ref<const VectorXd> &warm_start) {
      nlp->set_warm_start(warm_start);
    }

    MatrixXd get_dense_jacobian(const Ref<const VectorXd> &x) {
      MatrixXd J(nlp->get_num_constraints(), nlp->get_num_decision_vars());
      VectorXd g(nlp->get_num_constraints());
      VectorXd gplus(nlp->get_num_constraints());
      VectorXd xplus(x);
      nlp->eval_g(nlp->get_num_decision_vars(), xplus.data(), true, nlp->get_num_constraints(), g.data());
      double eps = 1e-4;
      for (int i = 0; i < nlp->get_num_decision_vars(); ++i) {
        xplus(i) += eps;
        nlp->eval_g(nlp->get_num_decision_vars(), xplus.data(), true, nlp->get_num_constraints(), gplus.data());
        xplus(i) = x(i);

        J.col(i) = (gplus - g)/eps;
      }
      return J;
    }

    VectorXd get_gradient(const Ref<const VectorXd> &x) {
      VectorXd grad(x.size());
      nlp->eval_grad_f(nlp->get_num_decision_vars(), x.data(), true, grad.data());
      return grad;
    }

    VectorXd get_hessian(const Ref<const VectorXd> &x) {
      VectorXd ret(nlp->get_Hrows().size());
      nlp->eval_h(
         nlp->get_num_decision_vars(),
         x.data(),
         true,
         1.,
         nlp->get_num_constraints(),
         NULL,
         true,
         nlp->get_Hrows().size(),
         NULL,
         NULL,
         ret.data()
      );
      return ret;
    }

    double solve(Ref<RowMatrixXd> trajectory, bool get_trajectory, bool &restoration_invoked, bool &restoration_failed) {
      // Ask Ipopt to solve the problem
      Ipopt::ApplicationReturnStatus status;
      status = app->OptimizeTNLP(nlp);

      if (get_trajectory) {
        if (trajectory.rows() != num_targets) {
          throw std::runtime_error("trajectory does not have num_targets rows");
        }
        if (trajectory.cols() != 1 + dim_q) {
          throw std::runtime_error("trajectory does not have 1 + dim_q columns");
        }
        const VectorXd &soln = nlp->get_soln();
        // Get arrival times and configurations
        for (int seq_idx = 0; seq_idx < num_targets; ++seq_idx) {
          trajectory.row(seq_idx) = soln.segment((1 + dim_q)*seq_idx, 1 + dim_q);
        }
      }

      restoration_invoked = nlp->get_restoration_invoked();
      restoration_failed = status == Ipopt::Restoration_Failed;

      if (status == Ipopt::Solve_Succeeded)
      {
         return nlp->get_cost();
      }
      else
      {
         // throw std::runtime_error("Ipopt failed");
         return std::numeric_limits<double>::infinity();
      }
      
      // As the SmartPtrs go out of scope, the reference count
      // will be decremented and the objects will automatically
      // be deleted.
    }
  private:
    Ipopt::SmartPtr<RobotArmNLP> nlp;
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app;
    int num_targets;
};

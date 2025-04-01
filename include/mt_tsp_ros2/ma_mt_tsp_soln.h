#pragma once
#include <Eigen/Dense>

using namespace Eigen;
typedef Matrix<long, Dynamic, Dynamic, RowMajor> RowMatrixXl;
typedef Matrix<long, Dynamic, 1> VectorXl;

class MAMTTSPSoln {
  public:
    MAMTTSPSoln(int num_agents, int num_targets) {
      // Needs to include a row for each depot and one row per target
      cell_seqs = RowMatrixXl::Zero(num_agents + num_targets, 2);
      time_seqs = VectorXl::Zero(num_agents + num_targets);
      agent_cell_seq_start_ptr = VectorXl::Zero(num_agents + 1);
    }

    Ref<RowMatrixXl> get_cell_seq(int agent) {
      int num_rows = agent_cell_seq_start_ptr(agent + 1) - agent_cell_seq_start_ptr(agent);
      return cell_seqs.block(agent_cell_seq_start_ptr(agent), 0, num_rows, 2);
    }

    Ref<VectorXl> get_time_seq(int agent) {
      int num_rows = agent_cell_seq_start_ptr(agent + 1) - agent_cell_seq_start_ptr(agent);
      return time_seqs.segment(agent_cell_seq_start_ptr(agent), num_rows);
    }

    RowMatrixXl &get_cell_seqs() {
      return cell_seqs;
    }

    VectorXl &get_time_seqs() {
      return time_seqs;
    }

    VectorXl &get_agent_cell_seq_start_ptr() {
      return agent_cell_seq_start_ptr;
    }

  private:
    RowMatrixXl cell_seqs;
    VectorXl time_seqs;
    VectorXl agent_cell_seq_start_ptr; // For convenience, has length num_agents + 1, where the final element equals the number of rows in cell_seqs
};

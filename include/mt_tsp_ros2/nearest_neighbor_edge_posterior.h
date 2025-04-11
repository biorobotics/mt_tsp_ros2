#include <ompl/datastructures/NearestNeighborsGNAT.h>
#include <omp.h>
#include <random>
#include <cmath>

using namespace Eigen;

typedef Matrix<long, Dynamic, 1> VectorXl;
typedef Matrix<bool, Dynamic, 1> VectorXb;

class NNItem {
  public:
    NNItem(VectorXd pt_pair, double cost) : pt_pair(pt_pair), cost(cost) {
    }

    const VectorXd &get_pt_pair() const {
      return pt_pair;
    }

    double get_cost() const {
      return cost;
    }

  private:
    VectorXd pt_pair;
    double cost;
};

typedef std::shared_ptr<NNItem> NNItemPtr;

class NearestNeighborEdgePosterior {
  public:
    NearestNeighborEdgePosterior(double eta, double cost_multiplier, bool ignore_last_elem, int num_threads) : eta(eta), cost_multiplier(cost_multiplier), ignore_last_elem(ignore_last_elem), uniform_dist(0, 1), num_threads(num_threads) {
      std::vector<std::uint32_t> seeds(num_threads);
      this_seed_seq.generate(seeds.begin(), seeds.end());
      for (int thread_idx = 0; thread_idx < num_threads; ++thread_idx) {
        rngs_per_thread.push_back(std::mt19937(seeds[thread_idx]));
      }
      nn.setDistanceFunction([this](const NNItemPtr a, const NNItemPtr b) { return (a->get_pt_pair() - b->get_pt_pair()).norm(); });
    }

    void add_edge(const Ref<const VectorXd> &pt1, const Ref<const VectorXd> &pt2, double cost) {
      VectorXd pt_pair(2*pt1.size() - 2*ignore_last_elem);
      pt_pair.head(pt1.size() - ignore_last_elem) = pt1.head(pt1.size() - ignore_last_elem);
      pt_pair.tail(pt2.size() - ignore_last_elem) = pt2.head(pt2.size() - ignore_last_elem);
      nn.add(std::make_shared<NNItem>(pt_pair, cost));
    }

    // I wrote the python code so that the gtsp_cost_mat_no_inf really is column major, which is why I don't say RowMajor for it
    void sample(const Ref<const VectorXl> &pt_to_target_ptr, 
                const Ref<const VectorXb> &evaluated_edge_mat,
                const Ref<const Matrix<bool, Dynamic, Dynamic, RowMajor>> &inf_mask,
                const Ref<const Matrix<double, Dynamic, Dynamic, RowMajor>> &all_pts,
                const Ref<const Matrix<double, Dynamic, Dynamic, RowMajor>> &gtsp_cost_mat,
                Ref<Matrix<long, Dynamic, Dynamic>> gtsp_cost_mat_no_inf,
                long inf_sub) {
      omp_set_num_threads(num_threads);
      int num_pts = all_pts.rows();
      #pragma omp parallel for collapse(2)
      for (int node_idx1 = 0; node_idx1 < num_pts; ++node_idx1) {
        for (int node_idx2 = 0; node_idx2 < num_pts; ++node_idx2) {
          if (pt_to_target_ptr(node_idx1) == pt_to_target_ptr(node_idx2) || 
              evaluated_edge_mat(node_idx1*num_pts + node_idx2) ||
              inf_mask(node_idx1, node_idx2) ||
              node_idx2 == 0) {
            continue;
          }

          VectorXd pt_pair(2*all_pts.cols() - 2*ignore_last_elem);
          pt_pair.head(all_pts.cols() - ignore_last_elem) = all_pts.row(node_idx1).transpose().head(all_pts.cols() - ignore_last_elem);
          pt_pair.tail(all_pts.cols() - ignore_last_elem) = all_pts.row(node_idx2).transpose().head(all_pts.cols() - ignore_last_elem);

          NNItemPtr nearest_neighbor = nn.nearest(std::make_shared<NNItem>(pt_pair, 0.));
          const VectorXd &nearest_neighbor_pt_pair = nearest_neighbor->get_pt_pair();
          double nearest_neighbor_cost = nearest_neighbor->get_cost();
          if (nearest_neighbor_cost < gtsp_cost_mat(node_idx1, node_idx2)) {
            continue;
          }
          double nearest_neighbor_dist = (nearest_neighbor_pt_pair - pt_pair).norm();
          double prob_sample_from_neighbor = exp(-eta*nearest_neighbor_dist);

          if (uniform_dist(rngs_per_thread[omp_get_thread_num()]) < prob_sample_from_neighbor) {
            if (std::isinf(nearest_neighbor_cost)) {
              gtsp_cost_mat_no_inf(node_idx1, node_idx2) = inf_sub;
              continue;
            }
            double scale = std::abs(nearest_neighbor_cost - gtsp_cost_mat(node_idx1, node_idx2))/3;
            gtsp_cost_mat_no_inf(node_idx1, node_idx2) = std::round(cost_multiplier*(nearest_neighbor_cost + scale*normal_dist(rngs_per_thread[omp_get_thread_num()])));

          } else {
            gtsp_cost_mat_no_inf(node_idx1, node_idx2) = std::round(cost_multiplier*gtsp_cost_mat(node_idx1, node_idx2));
          }
        }
      }
    }
  private:
    double eta;
    double cost_multiplier;
    bool ignore_last_elem;
    std::seed_seq this_seed_seq;
    int num_threads;
    std::vector<std::mt19937> rngs_per_thread;
    std::uniform_real_distribution<> uniform_dist;
    std::uniform_real_distribution<> normal_dist;
    ompl::NearestNeighborsGNAT<NNItemPtr> nn;
};

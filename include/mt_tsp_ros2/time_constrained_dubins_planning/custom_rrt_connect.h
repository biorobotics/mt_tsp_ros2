#include <ompl/geometric/planners/rrt/RRTConnect.h>
#include <ompl/base/goals/GoalSampleableRegion.h>
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_motion_validator.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_motion_validator_rects.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/NearestNeighborsSqrtApproxReturnDistance.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/NearestNeighborsSortByTime.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_time_state_space.h"
#include <ompl/tools/config/SelfConfig.h>
#include <unordered_set>
#include <chrono>

namespace og = ompl::geometric;
namespace ob = ompl::base;
using namespace og;

class DubinsSegmentChain : public og::PathGeometric {
  public:
    DubinsSegmentChain(const ob::SpaceInformationPtr &si);
    std::vector<RowMatrixXd> turns_chain;
};

class CustomRRTConnect : public og::RRTConnect {
  public:
    /** \brief Constructor */
    CustomRRTConnect(const ob::SpaceInformationPtr &si, bool nn_sort_by_time, bool use_approx_dist, bool addIntermediateStates = false);

    ob::PlannerStatus solve(const ob::PlannerTerminationCondition &ptc) override;

    double get_sampler_path_elongation_intervals_time();

    double get_sampling_time();

    int get_num_samples();

    double get_nearest_neighbor_time();

    int get_num_tree_nodes();

    double get_add_to_tree_time();

  protected:
    struct pair_hash {
      std::size_t operator()(const std::pair<Motion*, Motion*> &p) const;
    };

    // std::unordered_set<std::pair<Motion*, Motion*>, pair_hash> grow_pairs; // To stop advance from going in an infinite loop

    GrowState growTree(TreeData &tree, TreeGrowingInfo &tgi, Motion *rmotion, bool try_connect = false);

    double approx_distanceFunction(const Motion *a, const Motion *b) const;

    void setup() override;

    class DubinsMotion : public Motion
    {
    public:
        DubinsMotion();

        DubinsMotion(const ob::SpaceInformationPtr &si);

        RowMatrixXd turns;
    };

    double sampling_time;
    int num_samples;
    double add_to_tree_time;
    double grow_time;
    double motion_check_time;
    double nn_time;

    bool nn_sort_by_time;
    bool use_approx_dist;
};

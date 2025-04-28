#include <ompl/geometric/planners/rrt/RRTConnect.h>
#include <ompl/base/goals/GoalSampleableRegion.h>
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_motion_validator.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_motion_validator_rects.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/NearestNeighborsSqrtApproxReturnDistance.h"
//#include "mt_tsp_ros2/time_constrained_dubins_planning/NearestNeighborsSortByTime.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_time_state_space.h"
#include <ompl/tools/config/SelfConfig.h>
#include <unordered_set>
#include <chrono>

namespace og = ompl::geometric;
namespace ob = ompl::base;
using namespace og;

class DubinsSegmentChain : public og::PathGeometric {
  public:
    DubinsSegmentChain(const ob::SpaceInformationPtr &si) : og::PathGeometric(si) {
    }
    std::vector<RowMatrixXd> turns_chain;
};

class CustomRRTConnect : public og::RRTConnect {
  public:
    /** \brief Constructor */
    CustomRRTConnect(const ob::SpaceInformationPtr &si, bool addIntermediateStates = false) : og::RRTConnect(si, addIntermediateStates) {
      rng_ = ompl::RNG(3);
    }

    ob::PlannerStatus solve(const ob::PlannerTerminationCondition &ptc) override {
      checkValidity();
      auto *goal = dynamic_cast<ob::GoalSampleableRegion *>(pdef_->getGoal().get());

      if (goal == nullptr)
      {
          OMPL_ERROR("%s: Unknown type of goal", getName().c_str());
          throw std::runtime_error("Unknown type of goal");
          return ob::PlannerStatus::UNRECOGNIZED_GOAL_TYPE;
      }

      while (const ob::State *st = pis_.nextStart())
      {
          // Anoop: changed Motion to DubinsMotion
          auto *motion = new DubinsMotion(si_);
          si_->copyState(motion->state, st);
          motion->root = motion->state;
          tStart_->add(motion);
      }

      if (tStart_->size() == 0)
      {
          OMPL_ERROR("%s: Motion planning start tree could not be initialized!", getName().c_str());
          throw std::runtime_error("Motion planning start tree could not be initialized!");
          return ob::PlannerStatus::INVALID_START;
      }

      if (!goal->couldSample())
      {
          OMPL_ERROR("%s: Insufficient states in sampleable goal region", getName().c_str());
          throw std::runtime_error("Motion planning goal tree could not be initialized!");
          return ob::PlannerStatus::INVALID_GOAL;
      }

      if (!sampler_)
          sampler_ = si_->allocStateSampler();

      OMPL_INFORM("%s: Starting planning with %d states already in datastructure", getName().c_str(),
                  (int)(tStart_->size() + tGoal_->size()));

      TreeGrowingInfo tgi;
      tgi.xstate = si_->allocState();

      Motion *approxsol = nullptr;
      double approxdif = std::numeric_limits<double>::infinity();
      // Anoop: changed Motion to DubinsMotion
      auto *rmotion = new DubinsMotion(si_);
      ob::State *rstate = rmotion->state;
      bool solved = false;

      sampling_time = 0.;
      add_to_tree_time = 0.;
      motion_check_time = 0.;
      nn_time = 0.;

      while (!ptc)
      {
          TreeData &tree = startTree_ ? tStart_ : tGoal_;
          tgi.start = startTree_;
          startTree_ = !startTree_;
          TreeData &otherTree = startTree_ ? tStart_ : tGoal_;

          if (tGoal_->size() == 0 || pis_.getSampledGoalsCount() < tGoal_->size() / 2)
          {
              const ob::State *st = tGoal_->size() == 0 ? pis_.nextGoal(ptc) : pis_.nextGoal();
              if (st != nullptr)
              {
                  // Anoop: changed Motion to DubinsMotion
                  auto *motion = new DubinsMotion(si_);
                  si_->copyState(motion->state, st);
                  motion->root = motion->state;
                  tGoal_->add(motion);
              }

              if (tGoal_->size() == 0)
              {
                  OMPL_ERROR("%s: Unable to sample any valid states for goal tree", getName().c_str());
                  break;
              }
          }

          /* sample random state */
          auto timer_start = std::chrono::high_resolution_clock::now();
          sampler_->sampleUniform(rstate);
          auto timer_stop = std::chrono::high_resolution_clock::now();
          auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
          sampling_time += ((double)nanos)/1e9;

          GrowState gs = growTree(tree, tgi, rmotion);

          if (gs != TRAPPED)
          {
              /* remember which motion was just added */
              Motion *addedMotion = tgi.xmotion;

              /* attempt to connect trees */

              /* if reached, it means we used rstate directly, no need to copy again */
              if (gs != REACHED)
                  si_->copyState(rstate, tgi.xstate);

              tgi.start = startTree_;

              /* if initial progress cannot be done from the otherTree, restore tgi.start */

              GrowState gsc = growTree(otherTree, tgi, rmotion, true);
              if (gsc == TRAPPED)
                  tgi.start = !tgi.start;

              // Since we pass try_connect to growTree, we don't need to keep running EXTEND.
              // growTree does so internally. This relies on the fact that extending the nearest neighbor to rstate
              // can only produce closer states to rstate
              /*
              while (gsc == ADVANCED) {
                  gsc = growTree(otherTree, tgi, rmotion);
              }
              */

              /* update distance between trees */
              // const double newDist = tree->getDistanceFunction()(addedMotion, otherTree->nearest(addedMotion));
              // Anoop
              /*
              double newDist;
              std::static_pointer_cast<NearestNeighborsSqrtApproxReturnDistance<Motion*>>(otherTree)->nearest_and_distance(addedMotion, newDist, false);
              // std::static_pointer_cast<NearestNeighborsSortByTime<Motion*>>(otherTree)->nearest_and_distance(addedMotion, newDist, false);

              if (newDist < distanceBetweenTrees_)
              {
                  distanceBetweenTrees_ = newDist;
                  // OMPL_INFORM("Estimated distance to go: %f", distanceBetweenTrees_);
              }
              */

              Motion *startMotion = tgi.start ? tgi.xmotion : addedMotion;
              Motion *goalMotion = tgi.start ? addedMotion : tgi.xmotion;

              /* if we connected the trees in a valid way (start and goal pair is valid)*/
              if (gsc == REACHED && goal->isStartGoalPairValid(startMotion->root, goalMotion->root))
              {
                  // it must be the case that either the start tree or the goal tree has made some progress
                  // so one of the parents is not nullptr. We go one step 'back' to avoid having a duplicate state
                  // on the solution path
                  RowMatrixXd middle_turns;
                  if (startMotion->parent != nullptr) {
                      middle_turns = static_cast<DubinsMotion*>(startMotion)->turns;
                      startMotion = startMotion->parent;
                  } else {
                      middle_turns = static_cast<DubinsMotion*>(goalMotion)->turns;
                      goalMotion = goalMotion->parent;
                  }

                  connectionPoint_ = std::make_pair(startMotion->state, goalMotion->state);

                  /* construct the solution path */
                  Motion *solution = startMotion;
                  std::vector<Motion *> mpath1;
                  while (solution != nullptr)
                  {
                      mpath1.push_back(solution);
                      solution = solution->parent;
                  }

                  solution = goalMotion;
                  std::vector<Motion *> mpath2;
                  while (solution != nullptr)
                  {
                      mpath2.push_back(solution);
                      solution = solution->parent;
                  }

                  /*
                  auto path(std::make_shared<PathGeometric>(si_));
                  path->getStates().reserve(mpath1.size() + mpath2.size());
                  for (int i = mpath1.size() - 1; i >= 0; --i)
                      path->append(mpath1[i]->state);
                  for (auto &i : mpath2)
                      path->append(i->state);
                  */

                  // Anoop 
                  auto path(std::make_shared<DubinsSegmentChain>(si_));
                  path->getStates().reserve(mpath1.size() + mpath2.size());
                  for (int i = mpath1.size() - 1; i >= 0; --i) {
                      path->append(mpath1[i]->state);
                      if (static_cast<DubinsMotion*>(mpath1[i])->turns(0, 1) != -1) {
                        path->turns_chain.push_back(static_cast<DubinsMotion*>(mpath1[i])->turns);
                      }
                  }
                  path->turns_chain.push_back(middle_turns);
                  for (auto &i : mpath2) {
                      path->append(i->state);
                      if (static_cast<DubinsMotion*>(i)->turns(0, 1) != -1) {
                        path->turns_chain.push_back(static_cast<DubinsMotion*>(i)->turns);
                      }
                  }

                  pdef_->addSolutionPath(path, false, 0.0, getName());
                  solved = true;
                  break;
              }
              else
              {
                  // We didn't reach the goal, but if we were extending the start
                  // tree, then we can mark/improve the approximate path so far.
                  if (tgi.start)
                  {
                      // We were working from the startTree.
                      double dist = 0.0;
                      goal->isSatisfied(tgi.xmotion->state, &dist);
                      if (dist < approxdif)
                      {
                          approxdif = dist;
                          approxsol = tgi.xmotion;
                      }
                  }
              }
          }
      }

      si_->freeState(tgi.xstate);
      si_->freeState(rstate);
      delete rmotion;

      OMPL_INFORM("%s: Created %u states (%u start + %u goal)", getName().c_str(), tStart_->size() + tGoal_->size(),
                  tStart_->size(), tGoal_->size());

      if (approxsol && !solved)
      {
          /* construct the solution path */
          std::vector<Motion *> mpath;
          while (approxsol != nullptr)
          {
              mpath.push_back(approxsol);
              approxsol = approxsol->parent;
          }

          /*
          auto path(std::make_shared<PathGeometric>(si_));
          for (int i = mpath.size() - 1; i >= 0; --i)
              path->append(mpath[i]->state);
          */

          // Anoop
          auto path(std::make_shared<DubinsSegmentChain>(si_));
          for (int i = mpath.size() - 1; i >= 0; --i)
              path->append(mpath[i]->state);

          pdef_->addSolutionPath(path, true, approxdif, getName());
          // std::cout << "motion check time " << motion_check_time << std::endl;
          // std::cout << "nn time " << nn_time << std::endl;
          return ob::PlannerStatus::APPROXIMATE_SOLUTION;
      }

      return solved ? ob::PlannerStatus::EXACT_SOLUTION : ob::PlannerStatus::TIMEOUT;
    }

    double get_sampler_path_elongation_intervals_time() {
      return std::static_pointer_cast<DubinsTimeStateSampler>(sampler_)->get_path_elongation_intervals_time();
    }

    double get_sampling_time() {
      return sampling_time;
    }

    double get_nearest_neighbor_time() {
      return nn_time;
    }

    int get_num_tree_nodes() {
      return tStart_->size() + tGoal_->size();
    }

    double get_add_to_tree_time() {
      return add_to_tree_time;
    }

  protected:
    struct pair_hash {
      std::size_t operator()(const std::pair<Motion*, Motion*> &p) const {
        return (std::hash<Motion*>()(p.first) + 0x9e3779b9) ^ (std::hash<Motion*>()(p.second) + 0x9e3779b9);
      }
    };

    // std::unordered_set<std::pair<Motion*, Motion*>, pair_hash> grow_pairs; // To stop advance from going in an infinite loop

    GrowState growTree(TreeData &tree, TreeGrowingInfo &tgi, Motion *rmotion, bool try_connect = false) {
      /* find closest state in the tree */
      auto timer_start2 = std::chrono::high_resolution_clock::now();
      double dist;
      Motion *nmotion = std::static_pointer_cast<NearestNeighborsSqrtApproxReturnDistance<Motion*>>(tree)->nearest_and_distance(rmotion, dist, !try_connect);
      // Motion *nmotion = std::static_pointer_cast<NearestNeighborsSqrtApproxReturnDistance<Motion*>>(tree)->nearest_and_distance(rmotion, dist, false);
      // Motion *nmotion = std::static_pointer_cast<NearestNeighborsSortByTime<Motion*>>(tree)->nearest_and_distance(rmotion, dist, !try_connect);
      if (std::isinf(dist)) {
        return TRAPPED;
      }
      auto timer_stop2 = std::chrono::high_resolution_clock::now();
      auto nanos2 = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop2 - timer_start2).count();
      nn_time += ((double)nanos2)/1e9;

      /*
      if (grow_pairs.find(std::make_pair(nmotion, rmotion)) != grow_pairs.end()) {
        return TRAPPED;
      }
      grow_pairs.insert(std::make_pair(nmotion, rmotion));
      */

      /* assume we can reach the state we go towards */
      bool reach = true;

      /* find state to add */
      ob::State *dstate = rmotion->state;
      
      /*
      double d = si_->distance(nmotion->state, rmotion->state);
      if (d > maxDistance_)
      {
          si_->getStateSpace()->interpolate(nmotion->state, rmotion->state, maxDistance_ / d, tgi.xstate);

          // Check if we have moved at all. Due to some stranger state spaces (e.g., the constrained state spaces),
          // interpolate can fail and no progress is made. Without this check, the algorithm gets stuck in a loop as it
          // thinks it is making progress, when none is actually occurring.
          if (si_->equalStates(nmotion->state, tgi.xstate))
              return TRAPPED;

          dstate = tgi.xstate;
          reach = false;
      }

      bool validMotion = tgi.start ? si_->checkMotion(nmotion->state, dstate) :
                                     si_->isValid(dstate) && si_->checkMotion(dstate, nmotion->state);
      */

      // Anoop
      /*
      double d;
      if (tgi.start) {
        d = si_->distance(nmotion->state, rmotion->state);
      } else {
        d = si_->distance(rmotion->state, nmotion->state);
      }

      if (d > maxDistance_)
      {
          if (tgi.start) {
            si_->getStateSpace()->interpolate(nmotion->state, rmotion->state, maxDistance_ / d, tgi.xstate);
          } else {
            si_->getStateSpace()->interpolate(rmotion->state, nmotion->state, (1 - maxDistance_ / d), tgi.xstate);
          }

          // Check if we have moved at all. Due to some stranger state spaces (e.g., the constrained state spaces),
          // interpolate can fail and no progress is made. Without this check, the algorithm gets stuck in a loop as it
          // thinks it is making progress, when none is actually occurring.
          if (si_->equalStates(nmotion->state, tgi.xstate))
              return TRAPPED;

          dstate = tgi.xstate;
          reach = false;
          std::cout << "interpolated" << std::endl;
          std::cout << tgi.start << std::endl;
          double x1 = nmotion->state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getX();
          double y1 = nmotion->state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getY();
          double theta1 = nmotion->state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getYaw();
          double t1 = nmotion->state->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position;

          double x2 = rmotion->state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getX();
          double y2 = rmotion->state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getY();
          double theta2 = rmotion->state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getYaw();
          double t2 = rmotion->state->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position;

          double x3 = dstate->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getX();
          double y3 = dstate->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getY();
          double theta3 = dstate->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->getYaw();
          double t3 = dstate->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position;

          std::cout << x1 << " " << y1 << " " << theta1 << " " << t1 << std::endl;
          std::cout << x2 << " " << y2 << " " << theta2 << " " << t2 << std::endl;
          std::cout << x3 << " " << y3 << " " << theta3 << " " << t3 << std::endl;
      }

      std::cout << "checking valid" << std::endl;
      bool validMotion = tgi.start ? si_->checkMotion(nmotion->state, dstate) :
                                     si_->isValid(dstate) && si_->checkMotion(dstate, nmotion->state);
      std::cout << "done checking valid" << std::endl;
      */

      std::shared_ptr<DubinsMotionValidator> motionValidator = std::static_pointer_cast<DubinsMotionValidator>(si_->getMotionValidator());

      double minT = si_->getStateSpace()->as<DubinsTimeStateSpace>()->as<ob::TimeStateSpace>(1)->getMinTimeBound();
      double maxT = si_->getStateSpace()->as<DubinsTimeStateSpace>()->as<ob::TimeStateSpace>(1)->getMaxTimeBound();
      double maxDuration = 0.2*(maxT - minT);

      if (try_connect) {
        std::vector<RowMatrixXd> turns_vec;
        std::vector<Vector4d> states_vec;
        if (tgi.start) {
          std::static_pointer_cast<DubinsMotionValidatorRects>(motionValidator)->checkMotionForward_connect(turns_vec, states_vec, nmotion->state, dstate, maxDuration, reach);
        } else {
          std::static_pointer_cast<DubinsMotionValidatorRects>(motionValidator)->checkMotionBackward_connect(turns_vec, states_vec, dstate, nmotion->state, maxDuration, reach);
        }
        if (turns_vec.size() == 0) {
          return TRAPPED;
        }

        auto timer_start = std::chrono::high_resolution_clock::now();
        Motion *prev_motion = nmotion;
        for (int i = 0; i < turns_vec.size(); ++i) {
          auto *motion = new DubinsMotion(si_);
          motion->state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setX(states_vec[i](0));
          motion->state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setY(states_vec[i](1));
          motion->state->as<ob::CompoundState>()->as<ob::SE2StateSpace::StateType>(0)->setYaw(states_vec[i](2));
          motion->state->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position = states_vec[i](3);

          motion->parent = prev_motion;
          motion->root = nmotion->root;
          motion->turns = turns_vec[i];
          tree->add(motion);
          prev_motion = motion;
          /*
          if ((reach && turns_vec[i].col(1).sum() > maxDuration*5) ||
              (!reach && std::abs(turns_vec[i].col(1).sum() - maxDuration*5) > 1e-10)) {
            throw std::runtime_error("Incorrect segment dist");
          }
          */
        }

        // Need to use copyState as opposed to tgi.xstate = prev_motion->state,
        // or else we'll overwrite something in the tree
        // next time we update the elements of tgi.xstate
        si_->copyState(tgi.xstate, prev_motion->state);
        tgi.xmotion = prev_motion;

        auto timer_stop = std::chrono::high_resolution_clock::now();
        auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
        add_to_tree_time += ((double)nanos)/1e9;
      } else {
        bool validMotion;
        RowMatrixXd turns;
        // auto timer_start1 = std::chrono::high_resolution_clock::now();
        if (tgi.start) {
          turns = motionValidator->checkMotionForward(nmotion->state, dstate, maxDuration, tgi.xstate, reach, validMotion, try_connect);
        } else {
          turns = motionValidator->checkMotionBackward(dstate, nmotion->state, maxDuration, tgi.xstate, reach, validMotion, try_connect);
        }
        // auto timer_stop1 = std::chrono::high_resolution_clock::now();
        // auto nanos1 = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop1 - timer_start1).count();
        // motion_check_time += ((double)nanos1)/1e9;
        dstate = tgi.xstate;

        if (!validMotion)
            return TRAPPED;

        // assert(si_->isValid(dstate));

        // Anoop: took out addIntermediateStates_ case, and use DubinsMotion instead of Motion so I can add the turns
        assert(!addIntermediateStates_);
        auto timer_start = std::chrono::high_resolution_clock::now();
        auto *motion = new DubinsMotion(si_);
        si_->copyState(motion->state, dstate);
        motion->parent = nmotion;
        motion->root = nmotion->root;
        motion->turns = turns;
        tree->add(motion);
        auto timer_stop = std::chrono::high_resolution_clock::now();
        auto nanos = std::chrono::duration_cast<std::chrono::nanoseconds>(timer_stop - timer_start).count();
        add_to_tree_time += ((double)nanos)/1e9;

        tgi.xmotion = motion;
      }
      return reach ? REACHED : ADVANCED;
    }

    double approx_distanceFunction(const Motion *a, const Motion *b) const
    {
        return std::static_pointer_cast<DubinsTimeStateSpace>(si_->getStateSpace())->approx_distance(a->state, b->state);
    }

    void setup() override {
      // ompl::base::RRTConnect::setup();
      ompl::base::Planner::setup();
      ompl::tools::SelfConfig sc(si_, getName());
      sc.configurePlannerRange(maxDistance_);

      /*
      if (!tStart_)
          tStart_.reset(ompl::tools::SelfConfig::getDefaultNearestNeighbors<Motion *>(this));
      if (!tGoal_)
          tGoal_.reset(ompl::tools::SelfConfig::getDefaultNearestNeighbors<Motion *>(this));
      */
      /*
      if (!tStart_)
          tStart_.reset(new CustomNearestNeighborsSqrtApprox<Motion*>());
      if (!tGoal_)
          tGoal_.reset(new CustomNearestNeighborsSqrtApprox<Motion*>());
      */
      if (!tStart_)
          tStart_.reset(new NearestNeighborsSqrtApproxReturnDistance<Motion*>());
      if (!tGoal_)
          tGoal_.reset(new NearestNeighborsSqrtApproxReturnDistance<Motion*>());
      /*
      if (!tStart_)
          tStart_.reset(new NearestNeighborsSortByTime<Motion*>(true));
      if (!tGoal_)
          tGoal_.reset(new NearestNeighborsSortByTime<Motion*>(false));
      */
      tStart_->setDistanceFunction([this](const Motion *a, const Motion *b) { return distanceFunction(a, b); });
      tGoal_->setDistanceFunction([this](const Motion *a, const Motion *b) { return distanceFunction(b, a); });

      /*
      std::static_pointer_cast<NearestNeighborsSortByTime<Motion*>>(tStart_)->setApproxDistanceFunction([this](const Motion *a, const Motion *b) { return approx_distanceFunction(a, b); });
      std::static_pointer_cast<NearestNeighborsSortByTime<Motion*>>(tGoal_)->setApproxDistanceFunction([this](const Motion *a, const Motion *b) { return approx_distanceFunction(b, a); });
      */
      std::static_pointer_cast<NearestNeighborsSqrtApproxReturnDistance<Motion*>>(tStart_)->setApproxDistanceFunction([this](const Motion *a, const Motion *b) { return approx_distanceFunction(a, b); });
      std::static_pointer_cast<NearestNeighborsSqrtApproxReturnDistance<Motion*>>(tGoal_)->setApproxDistanceFunction([this](const Motion *a, const Motion *b) { return approx_distanceFunction(b, a); });
    }

  protected:
     class DubinsMotion : public Motion
     {
     public:
         DubinsMotion() {
           turns = -RowMatrixXd::Ones(1, 2);
         }

         DubinsMotion(const ob::SpaceInformationPtr &si) : Motion(si)
         {
           turns = -RowMatrixXd::Ones(1, 2);
         }

         RowMatrixXd turns;
     };

     double sampling_time;
     double add_to_tree_time;
     double grow_time;
     double motion_check_time;
     double nn_time;
};

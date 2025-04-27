#pragma once
 
#include <ompl/control/planners/rrt/RRT.h>
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_motion_validator_rects.h"
#include "mt_tsp_ros2/time_constrained_dubins_planning/dubins_control_sampler.h"
#include <ompl/control/spaces/RealVectorControlSpace.h>

using namespace ompl::control;
namespace ob = ompl::base;
 
class CustomControlRRT : public ompl::control::RRT
{
public:
    CustomControlRRT(const SpaceInformationPtr &si, std::shared_ptr<DubinsMotionValidatorRects> motion_validator) : ompl::control::RRT(si), motion_validator(motion_validator) {
      rng_ = ompl::RNG(1);
    }

    ob::PlannerStatus solve(const ob::PlannerTerminationCondition &ptc) override {
      checkValidity();
      ob::Goal *goal = pdef_->getGoal().get();
      auto *goal_s = dynamic_cast<ob::GoalSampleableRegion *>(goal);
  
      const ob::State *start_state;
      while (const ob::State *st = pis_.nextStart())
      {
          start_state = st; // I'm assuming there's only one start state so this is fine
          auto *motion = new Motion(siC_);
          si_->copyState(motion->state, st);
          siC_->nullControl(motion->control);
          nn_->add(motion);
      }
  
      if (nn_->size() == 0)
      {
          OMPL_ERROR("%s: There are no valid initial states!", getName().c_str());
          return ob::PlannerStatus::INVALID_START;
      }
  
      if (!sampler_)
          sampler_ = si_->allocStateSampler();
      if (!controlSampler_)
          controlSampler_ = siC_->allocDirectedControlSampler();
  
      OMPL_INFORM("%s: Starting planning with %u states already in datastructure", getName().c_str(), nn_->size());
  
      Motion *solution = nullptr;
      Motion *approxsol = nullptr;
      double approxdif = std::numeric_limits<double>::infinity();
  
      auto *rmotion = new Motion(siC_);
      ob::State *rstate = rmotion->state;
      Control *rctrl = rmotion->control;
      ob::State *xstate = si_->allocState();
  
      while (ptc == false)
      {
          /* sample random state (with goal biasing) */
          bool sample_goal = goal_s && rng_.uniform01() < goalBias_ && goal_s->canSample();
          if (sample_goal)
              goal_s->sampleGoal(rstate);
          else
              sampler_->sampleUniform(rstate);
  
          /* find closest state in the tree */
          Motion *nmotion = nn_->nearest(rmotion);
  
          /* sample a random control that attempts to go towards the random state, and also sample a control duration */

          bool goal_reached = false;
          unsigned int cd;
          if (sample_goal) {
            double start_t = start_state->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position;

            double t = nmotion->state->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position;

            double goal_t = rstate->as<ob::CompoundState>()->as<ob::TimeStateSpace::StateType>(1)->position;
            double delta_t = goal_t - t;
            if (0 <= delta_t && delta_t <= 0.2*(goal_t - start_t)) {
              bool reach;
              bool valid;
              motion_validator->checkMotionForward_internal(nmotion->state, rstate, std::numeric_limits<double>::infinity(), rmotion->state, reach, valid);
              if (valid) {
                cd = 1;
                // Just set the control to zero, we'll figure it out later
                RealVectorControlSpace::ControlType *real_vector_control = rmotion->control->as<RealVectorControlSpace::ControlType>();
                (*real_vector_control)[0] = 0.;
                (*real_vector_control)[1] = 0.;
                goal_reached = true;
              } else {
                cd = 0;
              }
            } else {
              cd = controlSampler_->sampleTo(rctrl, nmotion->control, nmotion->state, rmotion->state);
            }
          } else {
            cd = controlSampler_->sampleTo(rctrl, nmotion->control, nmotion->state, rmotion->state);
          }
  
          if (addIntermediateStates_)
          {
            throw std::runtime_error("addIntermediateStates not implemented");
          }
          else
          {
              if (cd >= siC_->getMinControlDuration())
              {
                  /* create a motion */
                  auto *motion = new Motion(siC_);
                  si_->copyState(motion->state, rmotion->state);
                  siC_->copyControl(motion->control, rctrl);
                  motion->steps = cd;
                  motion->parent = nmotion;
  
                  nn_->add(motion);
                  if (goal_reached)
                  {
                      solution = motion;
                      break;
                  }
              }
          }
      }
  
      bool solved = false;
      bool approximate = false;
      if (solution == nullptr)
      {
          solution = approxsol;
          approximate = true;
      }
  
      if (solution != nullptr)
      {
          lastGoalMotion_ = solution;
  
          /* construct the solution path */
          std::vector<Motion *> mpath;
          while (solution != nullptr)
          {
              mpath.push_back(solution);
              solution = solution->parent;
          }
  
          /* set the solution path */
          auto path(std::make_shared<PathControl>(si_));
          for (int i = mpath.size() - 1; i >= 0; --i)
              if (mpath[i]->parent)
                  path->append(mpath[i]->state, mpath[i]->control, mpath[i]->steps * siC_->getPropagationStepSize());
              else
                  path->append(mpath[i]->state);
          solved = true;
          pdef_->addSolutionPath(path, approximate, approxdif, getName());
      }
  
      if (rmotion->state)
          si_->freeState(rmotion->state);
      if (rmotion->control)
          siC_->freeControl(rmotion->control);
      delete rmotion;
      si_->freeState(xstate);
  
      OMPL_INFORM("%s: Created %u states", getName().c_str(), nn_->size());
  
      return {solved, approximate};
    }

    double get_collision_check_time() {
      return std::static_pointer_cast<DubinsControlSampler>(controlSampler_)->get_collision_check_time();
    }

protected:
    using ompl::control::RRT::rng_;
    std::shared_ptr<DubinsMotionValidatorRects> motion_validator;
};

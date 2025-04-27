#pragma once
#include "ompl/datastructures/NearestNeighborsSqrtApprox.h"
#include <algorithm>
#include <cmath>
#include <ompl/base/spaces/TimeStateSpace.h>
#include <ompl/base/StateSpace.h>
#include <ompl/base/State.h>

namespace ob = ompl::base;

template <typename _T>
class NearestNeighborsSortByTime : public ompl::NearestNeighborsSqrtApprox<_T>
{
public:
   NearestNeighborsSortByTime(bool reverse_sort) : sorted_data_(cmp(reverse_sort)) {
   }

   void add(const _T &data) override
   {
       sorted_data_.insert(data);
       updateCheckCount();
   }

   std::size_t size() const override
   {
       return sorted_data_.size();
   }

   virtual void setApproxDistanceFunction(const typename ompl::NearestNeighborsSqrtApprox<_T>::DistanceFunction &approx_distFun)
   {
       approx_distFun_ = approx_distFun;
   }

   _T nearest_and_distance(const _T &data, double &dmin, bool use_approx_dist) const
   {
       const std::size_t n = sorted_data_.size();

       if (checks_ > 0 && n > 0)
       {
           auto it = sorted_data_.begin();
           int k = 0;
           for (std::size_t j = 0; j < checks_; ++j)
           {
               std::size_t i = (j * checks_ + offset_) % n;
               for (; k < i; ++k) {
                 ++it;
               }

               double distance;
               if (use_approx_dist) {
                 distance = approx_distFun_(*it, data);
               } else {
                 distance = ompl::NearestNeighbors<_T>::distFun_(*it, data);
               }
               if (std::isfinite(distance)) {
                dmin = distance;
                return *it;
               }
           }
           offset_ = (offset_ + 1) % checks_;

           _T ret = *std::prev(sorted_data_.end());
           dmin = ompl::NearestNeighbors<_T>::distFun_(ret, data);
           return ret;
       }

       throw ompl::Exception("No elements found in nearest neighbors data structure");

       /*
       for (auto it = sorted_data_.begin(); it != sorted_data_.end(); ++it)
       {
           double distance = ompl::NearestNeighbors<_T>::distFun_(*it, data);
           if (std::isfinite(distance))
           {
               dmin = distance;
               return *it;
           }
       }
       if (sorted_data_.size()) {
         _T ret = *std::prev(sorted_data_.end());
         dmin = ompl::NearestNeighbors<_T>::distFun_(ret, data);
         return ret;
       }

       throw ompl::Exception("No elements found in nearest neighbors data structure");
       */
   }

protected:
  inline void updateCheckCount() 
  {
      checks_ = 1 + (std::size_t)floor(sqrt((double)sorted_data_.size()));
  }
  using ompl::NearestNeighborsSqrtApprox<_T>::checks_;
  using ompl::NearestNeighborsSqrtApprox<_T>::offset_;

  struct cmp {
    bool reverse_sort;
    cmp(bool reverse_sort) : reverse_sort(reverse_sort) {
    }
    bool operator() (_T a, _T b) const {
      double t_a = a->state->template as<ob::CompoundState>()->template as<ob::TimeStateSpace::StateType>(1)->position;
      double t_b = b->state->template as<ob::CompoundState>()->template as<ob::TimeStateSpace::StateType>(1)->position;
      return reverse_sort ? t_b < t_a : t_a < t_b;
    }
  };

  std::set<_T, cmp> sorted_data_;

  typename ompl::NearestNeighborsSqrtApprox<_T>::DistanceFunction approx_distFun_;
};

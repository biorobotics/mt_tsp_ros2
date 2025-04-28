#pragma once
#include "ompl/datastructures/NearestNeighborsSqrtApprox.h"
#include <algorithm>
#include <cmath>

template <typename _T>
class NearestNeighborsSqrtApproxReturnDistance : public ompl::NearestNeighborsSqrtApprox<_T>
{
public:
   virtual void setApproxDistanceFunction(const typename ompl::NearestNeighborsSqrtApprox<_T>::DistanceFunction &approx_distFun)
   {
       approx_distFun_ = approx_distFun;
   }

   _T nearest_and_distance(const _T &data, double &dmin, bool use_approx_dist) const
   {
       const std::size_t n = ompl::NearestNeighborsLinear<_T>::data_.size();
       std::size_t pos = n;

       if (checks_ > 0 && n > 0)
       {
           dmin = 0.0;
           for (std::size_t j = 0; j < checks_; ++j)
           {
               std::size_t i = (j * checks_ + offset_) % n;

               double distance;
               if (use_approx_dist) {
                 distance = approx_distFun_(ompl::NearestNeighborsLinear<_T>::data_[i], data);
               } else {
                 distance = ompl::NearestNeighbors<_T>::distFun_(ompl::NearestNeighborsLinear<_T>::data_[i], data);
               }
               if (pos == n || dmin > distance)
               {
                   pos = i;
                   dmin = distance;
               }
           }
           offset_ = (offset_ + 1) % checks_;
       }
       if (pos != n)
           return ompl::NearestNeighborsLinear<_T>::data_[pos];

       throw ompl::Exception("No elements found in nearest neighbors data structure");
       /*
       const std::size_t sz = data_.size();
       std::size_t pos = sz;
       dmin = 0.0;
       for (std::size_t i = 0; i < sz; ++i)
       {
           double distance;
           if (use_approx_dist) {
             distance = approx_distFun_(ompl::NearestNeighborsLinear<_T>::data_[i], data);
           } else {
             distance = ompl::NearestNeighbors<_T>::distFun_(data_[i], data);
           }

           if (pos == sz || dmin > distance)
           {
               pos = i;
               dmin = distance;
           }
       }
       if (pos != sz)
           return data_[pos];

       throw ompl::Exception("No elements found in nearest neighbors data structure");
       */
   }

protected:
  using ompl::NearestNeighborsSqrtApprox<_T>::checks_;
  using ompl::NearestNeighborsSqrtApprox<_T>::offset_;
  using ompl::NearestNeighborsLinear<_T>::data_;

  typename ompl::NearestNeighborsSqrtApprox<_T>::DistanceFunction approx_distFun_;
};

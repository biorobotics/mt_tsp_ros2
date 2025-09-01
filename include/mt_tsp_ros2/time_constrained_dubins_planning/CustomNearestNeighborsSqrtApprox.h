#include <ompl/datastructures/NearestNeighborsSqrtApprox.h>
#include <chrono>
#include <omp.h>
 
template <typename _T>
class CustomNearestNeighborsSqrtApprox : public ompl::NearestNeighborsSqrtApprox<_T>
{
public:
    _T nearest(const _T &data) const override {
      const std::size_t n = ompl::NearestNeighborsLinear<_T>::data_.size();
      std::size_t pos = n;

      if (checks_ > 0 && n > 0)
      {
          int num_threads = 8;
          std::vector<std::size_t> pos_per_thread(num_threads, 0);
          std::vector<double> dmin_per_thread(num_threads, std::numeric_limits<double>::infinity());
          omp_set_num_threads(num_threads);
          #pragma omp parallel for
          for (std::size_t j = 0; j < checks_; ++j)
          {
              std::size_t i = (j * checks_ + offset_) % n;

              double distance = ompl::NearestNeighbors<_T>::distFun_(ompl::NearestNeighborsLinear<_T>::data_[i], data);
              if (dmin_per_thread[omp_get_thread_num()] > distance)
              {
                  pos_per_thread[omp_get_thread_num()] = i;
                  dmin_per_thread[omp_get_thread_num()] = distance;
              }
          }
          pos = 0;
          double dmin = std::numeric_limits<double>::infinity();
          for (int i = 0; i < num_threads; ++i) {
            if (dmin_per_thread[i] < dmin) {
              pos = pos_per_thread[i];
              dmin = dmin_per_thread[i];
            }
          }
          offset_ = (offset_ + 1) % checks_;
      }
      if (pos != n)
          return ompl::NearestNeighborsLinear<_T>::data_[pos];

      throw ompl::Exception("No elements found in nearest neighbors data structure");
    }
protected:
  using ompl::NearestNeighborsSqrtApprox<_T>::checks_;
  using ompl::NearestNeighborsSqrtApprox<_T>::offset_;
};

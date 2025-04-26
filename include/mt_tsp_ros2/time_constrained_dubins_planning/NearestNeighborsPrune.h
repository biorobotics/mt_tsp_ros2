#include <ompl/datastructures/NearestNeighborsSqrtApprox.h>
#include <chrono>
 
template <typename _T>
class NearestNeighborsPrune : public ompl::NearestNeighborsSqrtApprox<_T>
{
public:
    _T nearest_and_distance(const _T &data, double &nearest_distance) const
    {
        const std::size_t n = ompl::NearestNeighborsLinear<_T>::data_.size();

        if (checks_ > 0 && n > 0)
        {
            nearest_distance = std::numeric_limits<double>::infinity();
            _T best = nullptr;
            std::deque<_T> queue;
            queue.push_back(ompl::NearestNeighborsLinear<_T>::data_[0]);
            for (std::size_t i = 0; i < checks_ && queue.size(); ++i) {
              std::size_t jlim = i == 0 ? offset_ : checks_;
              for (std::size_t j = 0; j < jlim && queue.size(); ++j) {
                _T pop = queue.front();
                queue.pop_front();
                for (_T child : pop->children) {
                  queue.push_back(child);
                }
              }
              if (queue.size() == 0) {
                break;
              }
              _T pop = queue.front();
              queue.pop_front();
              double distance = ompl::NearestNeighbors<_T>::distFun_(pop, data);
              if (std::isinf(distance)) {
                if (best == nullptr) {
                  best = pop;
                }
                continue;
              }
              if (distance < nearest_distance) {
                best = pop;
                nearest_distance = distance;
              }
              for (_T child : pop->children) {
                queue.push_back(child);
              }
            }
            offset_ = (offset_ + 1) % checks_;
            if (best == nullptr) {
              best = ompl::NearestNeighborsLinear<_T>::data_[0];
              nearest_distance = ompl::NearestNeighbors<_T>::distFun_(best, data);
            }
            return best;
        }

        throw ompl::Exception("No elements found in nearest neighbors data structure");
    }
protected:
  using ompl::NearestNeighborsSqrtApprox<_T>::checks_;
  using ompl::NearestNeighborsSqrtApprox<_T>::offset_;
};

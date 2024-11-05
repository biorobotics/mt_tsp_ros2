#pragma once
#include <Eigen/Dense>
#include <memory>

using namespace Eigen;

typedef VectorXi AStarCell;
typedef std::vector<AStarCell> AStarPath;

// Taken from https://jimmy-shen.medium.com/stl-map-unordered-map-with-a-vector-for-the-key-f30e5f670bae#:~:text=unordered_map%20uses%20vector%20as%20the%20key&text=You%20can%20use%20the%20following,make%20the%20best%20of%20STL.&text=%7D%3B,so%20that%20collisions%20are%20minimized
struct AStarCellHash {
  int operator()(const AStarCell &V) const {
    int hash = V.size();
    for(int i = 0; i < V.size(); ++i) {
      hash ^= V[i] + 0x9e3779b9 + (hash << 6) + (hash >> 2);
    }
    return hash;
  }
};

class AStarNode {
  private:
    // Cell
    AStarCell cell;

    // Cell of expanded predecessor with lowest g*
    AStarCell back;

    double g;
    double h;
    double eps;
    double eps_h; // Scaled heuristic value
    double f;

  public:
    /*
     * CONSTRUCTOR: initializes node object
     * ARGUMENTS
     * data: node data
     * back: cell of expanded predecessor with lowest g*
     * g: g* value of optimal predecessor so far
     * h: heuristic
     * eps: heuristic value
     */
    AStarNode(const AStarCell &back,
              double g, double h, double eps,
              const AStarCell &cell) : back(back),
                                       g(g), h(h), eps(eps), 
                                       eps_h(eps*h), f(g + eps_h),
                                       cell(cell) {}

    /*
     * get_cell: access cell
     * RETURN: cell
     */
    const AStarCell &get_cell() const {
      return cell;
    }

    /*
     * set_cell: set cell
     */
    void set_cell(const AStarCell &cell) {
      this->cell = cell;
    }

    /*
     * get_back: access cell of optimal predecessor
     * RETURN: cell of optimal predecessor
     */
    const AStarCell &get_back() {
      return back;
    }

    /*
     * get_g: access g value
     * RETURN: g value
     */
    double get_g() {
      return g;
    }

    /*
     * get_h: access h value
     * RETURN: h value
     */
    double get_h() {
      return h;
    }

    /*
     * get_eps_h: access h value
     * RETURN: eps_h value
     */
    double get_eps_h() {
      return eps_h;
    }

    /*
     * get_eps: access most recent inflation factor for heuristic
     * RETURN: eps value
     */
    double get_eps() {
      return eps;
    }

    /*
     * get_f: access f value
     * RETURN: f value
     */
    double get_f() {
      return f;
    }

    /*
     * inflate_h: update h value with new inflation factor
     * eps: new inflation factor for h
     */
    void inflate_h(double eps) {
      this->eps = eps;
      eps_h = eps*h;
      f = g + eps_h;
    }

    /*
     * set_h: update h value
     */
    void set_h(double h) {
      this->h = h;
      eps_h = eps*h;
      f = g + eps_h;
    }

    /*
     * update_path_to_node: updates the optimal predecessor and and cost. Effectively,
     * this corresponds to updating the path to the node
     * ARGUMENTS
     * back: new predecessor cell
     * g: g* for new predecessor
     */
    void update_path_to_node(const AStarCell &back, double g) {
      this->back = back;
      this->g = g;
      f = g + eps_h;
    }
};

typedef std::shared_ptr<AStarNode> AStarNodePtr;

struct compare_astar_nodes {
  bool operator() (AStarNodePtr &node1,
                   AStarNodePtr &node2) {
    return node1->get_f() > node2->get_f();
  }
};

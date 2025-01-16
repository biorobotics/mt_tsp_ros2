#pragma once
#include "mt_tsp_ros2/astar_problem.h"
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <vector>
#include <iostream>
#include <iomanip>
#include <memory>

class OpenListElement {
  private:
    AStarCell cell;
    double f;
  public:
    OpenListElement(const AStarCell &cell, double f) : cell(cell), f(f) {}

    const AStarCell& get_cell() const {
      return cell;
    }

    double get_f() const {
      return f;
    }
};

struct compare_open_list_elements {
  bool operator() (OpenListElement &elem1,
                   OpenListElement &elem2) {
    return elem1.get_f() > elem2.get_f();
  }
};

class ARAStarData {
  protected:
    std::priority_queue<OpenListElement, std::vector<OpenListElement>, compare_open_list_elements> open_list;
    std::unordered_map<AStarCell, AStarNodePtr, AStarCellHash> seen_nodes;
    std::unordered_set<AStarCell, AStarCellHash> closed_list;
    std::unordered_set<AStarCell, AStarCellHash> incons;

    std::shared_ptr<AStarProblem> problem;

    AStarCell start_cell;
    AStarNodePtr goal_node = nullptr;

    double eps;

    /* dust_off_open_list: when we find a better path to a node,
     * we insert a copy of it into the open list, which means 
     * there's a chance we run into an obsolete copy. Clear these out
     */
    void dust_off_open_list() {
      if (open_list.size() == 0) {
        return;
      }
      OpenListElement top = open_list.top();
      while (open_list.size() != 0 && is_closed(top.get_cell())) {
        open_list.pop();
        top = open_list.top();
      }
    }

  public:
    ARAStarData() {}

    // Adds start node to the open list
    ARAStarData(const std::shared_ptr<AStarProblem> &problem, const AStarCell &start_cell) : problem(problem), start_cell(start_cell), eps(1.0) {
      seen_nodes[start_cell] = std::make_shared<AStarNode>(start_cell,
                                                           0,
                                                           problem->heuristic(start_cell), eps,
                                                           start_cell);

      if (problem->is_goal(seen_nodes[start_cell])) {
        goal_node = seen_nodes[start_cell];
      }
      open(start_cell, seen_nodes[start_cell]->get_f());
    }

    ARAStarData(const AStarCell &start_cell) : start_cell(start_cell), eps(1.0) {
    }

    void print_open() {
      std::cout << "Printing open. Note that this empties the open list" << std::endl;
      while (open_list.size() != 0) {
        OpenListElement top = open_list.top();
        open_list.pop();
        std::cout << top.get_cell() << " " << top.get_f() << " " << seen_nodes[top.get_cell()]->get_g() << " " << seen_nodes[top.get_cell()]->get_h() << std::endl;
      }
    }

    std::shared_ptr<ARAStarData> copy() {
      std::shared_ptr<ARAStarData> copied_data = std::make_shared<ARAStarData>(start_cell);
      // std::shared_ptr<ARAStarData> copied_data = std::make_shared<ARAStarData>(problem, start_cell);
      copied_data->open_list = open_list;
      for (auto elem : seen_nodes) {
        assert(elem.second != nullptr);
        copied_data->seen_nodes[elem.first] = std::make_shared<AStarNode>(elem.second->get_back(),
                                                                          elem.second->get_g(),
                                                                          elem.second->get_h(),
                                                                          elem.second->get_eps(),
                                                                          elem.second->get_cell());
      }
      copied_data->closed_list = closed_list;
      copied_data->incons = incons;
      copied_data->goal_node = goal_node;
      return copied_data;
    }

    void nullify_problem() {
      problem = nullptr;
    }

    // Same start cell, but different goal and heuristic
    void update_problem(const std::shared_ptr<AStarProblem> &problem) {
      this->problem = problem;
      std::priority_queue<OpenListElement, std::vector<OpenListElement>, compare_open_list_elements> updated_open_list;
      while (open_list.size() != 0) {
        dust_off_open_list();
        if (open_list.size() == 0) {
          break;
        }
        AStarNodePtr node = seen_nodes[open_list.top().get_cell()];
        open_list.pop();
        node->set_h(problem->heuristic(node->get_cell()));
        updated_open_list.push(OpenListElement(node->get_cell(), node->get_f()));
      }
      assert(updated_open_list.size() != 0);
      open_list = updated_open_list;

      goal_node = nullptr;
    }

    void generate_successors(std::vector<AStarCell> &succ, std::vector<double> &transition_costs, const AStarNodePtr &node) {
      problem->generate_successors(succ, transition_costs, node);
    }

    const std::unordered_map<AStarCell, AStarNodePtr, AStarCellHash> &get_seen_nodes() {
      return seen_nodes;
    }

    /*
     * is_closed: tells whether this cell is closed
     * ARGUMENTS
     * cell: cell
     * RETURN: true if this cell is in the closed list, false if not
     */
    bool is_closed(AStarCell cell) {
      return closed_list.find(cell) != closed_list.end();
    }

     /*
     * reset: moves all inconsistent nodes back into the open list,
     * and updates heuristics of everything in the open list. Clears
     * the closed list
     * ARGUMENTS
     * eps: new heuristic inflation factor 
     */
    void reset(double eps) {
      this->eps = eps;
      while(open_list.size() != 0) {
        dust_off_open_list();
        make_incons(open_list.top().get_cell());
        open_list.pop();
      }

      closed_list.clear();
      for (std::unordered_set<AStarCell>::iterator it = incons.begin();
           it != incons.end(); ++it) {
        seen_nodes[*it]->inflate_h(eps);
        open(*it, seen_nodes[*it]->get_f());
      }
      incons.clear();
    }

    bool is_goal(const AStarNodePtr &node) {
      return problem->is_goal(node);
    }

    /*
     * update_path_to_node: updates the path to the node.
     * ARGUMENTS
     * problem: AStar problem
     * node: node
     */
    void update_path_to_node(const AStarNodePtr &back, const AStarCell &cell, double transition_cost) {
      double g = back->get_g() + transition_cost;
      if (!saw_cell(cell)) {
        // Since we haven't seen the node yet, its g value is infinity
        seen_nodes[cell] = std::make_shared<AStarNode>(back->get_cell(),
                                                       g,
                                                       problem->heuristic(cell), eps,
                                                       cell);

        if (problem->is_goal(seen_nodes[cell]) && (goal_node == nullptr || seen_nodes[cell]->get_g() < goal_node->get_g())) {
          goal_node = seen_nodes[cell];
        }
        open(cell, seen_nodes[cell]->get_f());
      } else if (g < seen_nodes[cell]->get_g()) {
        seen_nodes[cell]->update_path_to_node(back->get_cell(), g);
        // seen_nodes[cell] = node; // Do it this way if we want the updated data

        if (problem->is_goal(seen_nodes[cell]) && (goal_node == nullptr || seen_nodes[cell]->get_g() < goal_node->get_g())) {
          goal_node = seen_nodes[cell];
        }

        if (is_closed(cell)) {
          make_incons(cell);
        } else {
          open(cell, seen_nodes[cell]->get_f());
        }
      }
    }

    /*
     * open: pushes the node onto the open list
     * ARGUMENTS
     * node: node
     * REQUIRES: saw_cell(node->get_cell())
     */
    void open(const AStarCell &cell, double f) {
      open_list.push(OpenListElement(cell, f));
    }

    /*
     * make_incons: inserts the cell into the inconsistent set
     * ARGUMENTS
     * cell: cell to insert
     * REQUIRES: is_closed(cell)
     */
    void make_incons(const AStarCell &cell) {
      incons.insert(cell);
    }

    /*
     * get_f_goal: returns the f-value of the goal
     * RETURN: f-value of the goal
     */
    double get_f_goal() {
      if (goal_node == nullptr) {
        return std::numeric_limits<double>::infinity();
      }
      return goal_node->get_f();
    }

    /*
     * get_next: returns the node at the top of the open list
     * RETURN: node at the top of the open list
     * REQUIRES: !open_list_empty()
     */
    AStarNodePtr get_next() {
      dust_off_open_list();
      return seen_nodes[open_list.top().get_cell()];
    }

    /*
     * expand_next: pops off the node at the top of the open list
     * and returns it, inserting it into the closed list
     * RETURN: node previously at the top of the open list
     * REQUIRES: !open_list_empty()
     */
    AStarNodePtr expand_next() {
      dust_off_open_list();
      AStarNodePtr ret = seen_nodes[open_list.top().get_cell()];
      open_list.pop();
      closed_list.insert(ret->get_cell());
      return ret;
    }

    /*
     * open_list_empty: determines whether the open list is empty
     * RETURN: true if the open list is empty, false otherwise
     */
    bool open_list_empty() {
      dust_off_open_list();
      return open_list.size() == 0;
    }

    /*
     * saw_cell: return true if we've seen the node with the cell
     * ARGUMENTS
     * cell: the cell of the node
     * RETURN: true if we've seen the node, false otherwise
     */
    bool saw_cell(const AStarCell &cell) {
      return seen_nodes.find(cell) != seen_nodes.end();
    }

    /*
     * get_path: if we reached the goal, gets path to the goal.
     * Otherwise, gets a path to the node with the best h-value
     * ARGUMENTS
     * path: populated with path
     */
    void get_path(AStarPath &path) {
      int num_nodes = 1;

      AStarNodePtr terminal_node = goal_node;
      if (terminal_node == nullptr) {
        std::cout << "Did not reach goal" << std::endl;
        // If we didn't find the goal, get the path to the node with the best h-value
        double best_h = -1;
        for (auto it : seen_nodes) {
          if (best_h < 0 || it.second->get_h() < best_h) {
            best_h = it.second->get_h();
            terminal_node = it.second;
          }
        }
      } else {
        std::cout << "Goal node has path cost " << goal_node->get_g() << std::endl;
        assert(goal_node->get_g() == goal_node->get_f());
      }

      AStarCell cell = terminal_node->get_cell();
      AStarNodePtr node = terminal_node;
      while (!(node->get_cell() == start_cell)) {
        cell = node->get_back();
        node = seen_nodes[cell];
        ++num_nodes;
      }
      path.resize(num_nodes);

      cell = terminal_node->get_cell();
      node = terminal_node;
      int i = num_nodes - 1;
      while (!(node->get_cell() == start_cell)) {
        path[i] = node->get_cell();
        cell = node->get_back();
        node = seen_nodes[cell];
        --i;
      }
      path[i] = node->get_cell();
    }
};

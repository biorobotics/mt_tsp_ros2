#pragma once

#include <unordered_set>
#include <vector>
#include <array>
#include <queue>
#include <unordered_map>

using namespace std;

struct NodePos {
  int x;
  int y;
  int z;
  NodePos(int x, int y, int z) : x(x), y(y), z(z) {}
  
  bool operator==(const NodePos &other) const {
    return other.x == x && other.y == y && other.z == z;
  }
};

struct NodePosHash {
  size_t operator()(const NodePos& pos) const { 
    return (pos.x * 31) + (pos.y * 37) + (pos.z * 41);
  }
};

class Node {
  private:
    // Node position
    NodePos pos;

    // Position of expanded predecessor with lowest g*
    NodePos back;

    float g;
    float f;

    /*
     * set_predecessor: updates the optimal predecessor and cost
     * ARGUMENTS
     * bx, by: new predecessor location
     * g: g* for new predecessor
     */
    void set_predecessor(int bx, int by, int bz, float g) {
      f = f - this->g + g;
      this->g = g;
      back.x = bx;
      back.y = by;
      back.z = bz;
    }

    friend class AStarData3D;
    friend class GlennBFSData3D;

  public:
    /*
     * CONSTRUCTOR: initializes node object
     * ARGUMENTS
     * x, y, z: node location
     * bx, by, bz: location of expanded predecessor with lowest g*
     * g: g* value of optimal predecessor so far
     * h: heuristic
     */
    Node(int x, int y, int z, int bx, int by, int bz, float g, float h) : pos(x, y, z),
                                                                          back(bx, by, bz), g(g),
                                                                          f(g + h) {}

    Node() : pos(0, 0, 0),
             back(0, 0, 0), g(0),
             f(0) {}

    /*
     * get_pos: access position
     * RETURN: position 
     */
    NodePos get_pos() const {
      return pos;
    }

    /*
     * get_back: access position of optimal predecessor
     * RETURN: position of optimal predecessor
     */
    NodePos get_back() const {
      return back;
    }

    void set_back(const NodePos &back) {
      this->back = back;
    }

    /*
     * get_g: access g value
     * RETURN: g value
     */
    float get_g() const {
      return g;
    }

    /*
     * get_f: access f value
     * RETURN: f value
     */
    float get_f() const {
      return f;
    }

    /*
     * update_h: update h value
     * h: new heuristic
     */
    void update_h(float h) {
      f = h + g;
    }
};

struct compare_nodes {
  bool operator() (Node &node1,
                   Node &node2) {
    return node1.get_f() > node2.get_f();
  }
};

class AStarData3D {
  private:
    priority_queue<Node, vector<Node>, compare_nodes> open_list;
    unordered_map<NodePos, Node, NodePosHash> seen_nodes;
    unordered_set<NodePos, NodePosHash> closed_list;
    unordered_set<NodePos, NodePosHash> incons;

    int startx;
    int starty;
    int startz;
    int goalx;
    int goaly;
    int goalz;

    float eps;

    /* dust_off_open_list: when we find a better path to a node,
     * we insert a copy of it into the open list, which means 
     * there's a chance we run into an obsolete copy. Clear these out
     */
    void dust_off_open_list() {
      if (open_list.size() == 0) {
        return;
      }
      Node top = open_list.top();
      while (open_list.size() != 0 && is_closed(top.pos.x, top.pos.y, top.pos.z)) {
        open_list.pop();
        top = open_list.top();
      }
    }

    /*
     * heuristic: computes heuristic for a node
     * ARGUMENTS
     * goalposeX: goal x position
     * goalposeY: goal y position
     * goalposeZ: goal z position
     * x: node x position
     * y: node y position
     * z: node z position
     */
    float heuristic_8_connected_3d(int goalposeX, int goalposeY, int goalposeZ, int x, int y, int z) {
      float dx = (float)(std::abs(x - goalposeX));
      float dy = (float)(std::abs(y - goalposeY));
      float dz = (float)(std::abs(z - goalposeZ));
      float max = std::max(std::max(dx, dy), dz);
      float min = std::min(std::min(dx, dy), dz);
      float btw = dx + dy + dz - max - min;
      float d3 = min*sqrt(3);
      float d2 = (btw - min)*sqrt(2);
      float d1 = max - btw;
      return d1 + d2 + d3;
    }
  
  public:
    /*
     * CONSTRUCTOR: adds start position into the open list
     * ARGUMENTS
     * startx: initial robot x
     * starty: initial robot y
     * startz: initial robot z
     * goalx: target x
     * goaly: target y
     * goalz: target z
     * eps: inflation factor for heuristic
     */
    AStarData3D(int startx, int starty, int startz,
                int goalx, int goaly, int goalz, float eps) : startx(startx),
                                                              starty(starty),
                                                              startz(startz),
                                                              goalx(goalx),
                                                              goaly(goaly),
                                                              goalz(goalz),
                                                              eps(eps) {
      set_predecessor(startx, starty, startz, startx, starty, startz, 0);
      open(startx, starty, startz);
    }

    /*
     * is_closed: tells whether the node at this position is closed
     * ARGUMENTS
     * x, y, z: node position
     * RETURN: true if this node is in the closed list, false if not
     */
    bool is_closed(int x, int y, int z) {
      return closed_list.find(NodePos(x, y, z)) != closed_list.end();
    }

    /*
     * reset: moves all inconsistent states back into the open list,
     * and updates heuristics of everything in the open list. Clears
     * the closed list
     * ARGUMENTS
     * eps: new inflation factor to use when we update heuristics
     */
    void reset(float eps) {
      if (this->eps == eps) {
        return;
      }
      this->eps = eps;

      while(open_list.size() != 0) {
        dust_off_open_list();
        incons.insert(open_list.top().get_pos());
        open_list.pop();
      }

      closed_list.clear();
      for (unordered_set<NodePos, NodePosHash>::iterator it = incons.begin();
           it != incons.end(); ++it) {
        NodePos pos = *it;
        open(pos.x, pos.y, pos.z);
      }
      incons.clear();
    }

    /*
     * set_predecessor: sets the predecessor of the node at position
     * (x, y, z) to (bx, by, bz), and sets this node's g-value to g
     * ARGUMENTS
     * x, y, z: node position
     * bx, by: predecessor position
     * g: g-value of the node at (x, y, z) if we travel to it from (bx, by, bz)
     */
    void set_predecessor(int x, int y, int z, int bx, int by, int bz, float g) {
      NodePos pos(x, y, z);
      if (!saw_node(pos)) {
        seen_nodes[pos] = Node(x, y, z, bx, by, bz, g, eps*heuristic_8_connected_3d(goalx, goaly, goalz, x, y, z));
      } else {
        seen_nodes[pos].set_predecessor(bx, by, bz, g);
      }
    }

    /*
     * open: pushes the node at (x, y, z) into the open list
     * ARGUMENTS
     * x, y, z: node position
     * REQUIRES: saw_node(NodePos(x, y, z))
     */
    void open(int x, int y, int z) {
      NodePos pos(x, y, z);
      seen_nodes[pos].update_h(eps*heuristic_8_connected_3d(goalx, goaly, goalz, x, y, z));
      open_list.push(seen_nodes[pos]);
    }

    /*
     * make_incons: inserts the node at (x, y, z) into the inconsistent set
     * ARGUMENTS
     * x, y, z: node position
     * REQUIRES: is_closed(x, y, z)
     */
    void make_incons(int x, int y, int z) {
      incons.insert(NodePos(x, y, z));
    }

    /*
     * get_f: returns the f-value of the node at (x, y, z)
     * ARGUMENTS
     * x, y: node position
     * RETURN: f-value of the node
     */
    float get_f(int x, int y, int z) {
      NodePos pos(x, y, z);
      return saw_node(pos) ? seen_nodes[pos].get_f() : std::numeric_limits<float>::infinity();
    }

    /*
     * get_g: returns the g-value of the node at (x, y, z)
     * ARGUMENTS
     * x, y, z: node position
     * RETURN: g-value of the node
     */
    float get_g(int x, int y, int z) {
      NodePos pos(x, y, z);
      return saw_node(pos) ? seen_nodes[pos].get_g() : std::numeric_limits<float>::infinity();
    }

    /*
     * get_next: returns the node at the top of the open list
     * RETURN: node at the top of the open list
     * REQUIRES: !open_list_empty()
     */
    Node get_next() {
      dust_off_open_list();
      return open_list.top();
    }

    /*
     * expand_next: pops off the node at the top of the open list
     * and returns it, inserting it into the closed list
     * RETURN: node previously at the top of the open list
     * REQUIRES: !open_list_empty()
     */
    Node expand_next() {
      dust_off_open_list();
      Node ret = open_list.top();
      open_list.pop();
      closed_list.insert(ret.get_pos());
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
     * saw_node: return true if we've seen the node at pos
     * ARGUMENTS
     * pos: the position of the node
     * RETURN: true if we've seen the node, false otherwise
     */
    bool saw_node(const NodePos &pos) {
      return seen_nodes.find(pos) != seen_nodes.end();
    }

    /*
     * get_path: gets path
     * ARGUMENTS
     * path: populated with path
     */
    void get_path(vector<std::array<int, 3>> &path) {
      int num_nodes = 1;
      NodePos pos(goalx, goaly, goalz);
      Node node = seen_nodes[pos];
      while (pos.x != startx || pos.y != starty || pos.z != startz) {
        pos = node.get_back();
        node = seen_nodes[pos];
        ++num_nodes;
      }
      path.resize(num_nodes);

      pos = NodePos(goalx, goaly, goalz);
      node = seen_nodes[pos];
      int i = num_nodes - 1;
      while (pos.x != startx || pos.y != starty || pos.z != startz) {
        path[i][0] = pos.x;
        path[i][1] = pos.y;
        path[i][2] = pos.z;
        pos = node.get_back();
        node = seen_nodes[pos];
        --i;
      }
      path[i][0] = pos.x;
      path[i][1] = pos.y;
      path[i][2] = pos.z;
    }

    const std::unordered_map<NodePos, Node, NodePosHash> &get_seen_nodes() {
      return seen_nodes;
    }
};

class GlennBFSData3D {
  private:
    vector<Node> open_list1;
    vector<Node> open_list2;
    bool pop1;
    unordered_map<NodePos, Node, NodePosHash> seen_nodes;
    unordered_set<NodePos, NodePosHash> closed_list;

    int startx;
    int starty;
    int startz;
    int goalx;
    int goaly;
    int goalz;

  public:
    /*
     * CONSTRUCTOR: adds start position into the open list
     * ARGUMENTS
     * startx: initial robot x
     * starty: initial robot y
     * startz: initial robot z
     * goalx: target x
     * goaly: target y
     * goalz: target z
     */
    GlennBFSData3D(int startx, int starty, int startz,
              int goalx, int goaly, int goalz) : startx(startx),
                                                 starty(starty),
                                                 startz(startz),
                                                 goalx(goalx),
                                                 goaly(goaly),
                                                 goalz(goalz) {
      set_predecessor(startx, starty, startz, startx, starty, startz, 0);
      pop1 = true;
      open(startx, starty, startz);
    }

    /*
     * is_closed: tells whether the node at this position is closed
     * ARGUMENTS
     * x, y, z: node position
     * RETURN: true if this node is in the closed list, false if not
     */
    bool is_closed(int x, int y, int z) {
      return closed_list.find(NodePos(x, y, z)) != closed_list.end();
    }

    /*
     * set_predecessor: sets the predecessor of the node at position
     * (x, y, z) to (bx, by, bz), and sets this node's g-value to g
     * ARGUMENTS
     * x, y, z: node position
     * bx, by: predecessor position
     * g: g-value of the node at (x, y, z) if we travel to it from (bx, by, bz)
     */
    void set_predecessor(int x, int y, int z, int bx, int by, int bz, float g) {
      NodePos pos(x, y, z);
      if (!saw_node(pos)) {
        seen_nodes[pos] = Node(x, y, z, bx, by, bz, g, 0.);
      } else {
        seen_nodes[pos].set_predecessor(bx, by, bz, g);
      }
    }

    /*
     * open: pushes the node at (x, y, z) into the open list
     * ARGUMENTS
     * x, y, z: node position
     * REQUIRES: saw_node(NodePos(x, y, z))
     */
    void open(int x, int y, int z) {
      NodePos pos(x, y, z);
      if (pop1) {
        open_list2.push_back(seen_nodes[pos]);
      } else {
        open_list1.push_back(seen_nodes[pos]);
      }
    }

    /*
     * get_f: returns the f-value of the node at (x, y, z)
     * ARGUMENTS
     * x, y: node position
     * RETURN: f-value of the node
     */
    float get_f(int x, int y, int z) {
      NodePos pos(x, y, z);
      return saw_node(pos) ? seen_nodes[pos].get_f() : std::numeric_limits<float>::infinity();
    }

    /*
     * get_g: returns the g-value of the node at (x, y, z)
     * ARGUMENTS
     * x, y, z: node position
     * RETURN: g-value of the node
     */
    float get_g(int x, int y, int z) {
      NodePos pos(x, y, z);
      return saw_node(pos) ? seen_nodes[pos].get_g() : std::numeric_limits<float>::infinity();
    }

    /*
     * get_next: returns the node at the top of the open list
     * RETURN: node at the top of the open list
     * REQUIRES: !open_list_empty()
     */
    Node get_next() {
      if (pop1) {
        return open_list1.back();
      } else {
        return open_list2.back();
      }
    }

    /*
     * expand_next: pops off the node at the top of the open list
     * and returns it, inserting it into the closed list
     * RETURN: node previously at the top of the open list
     * REQUIRES: !open_list_empty()
     */
    Node expand_next() {
      if (pop1 && open_list1.size() == 0) {
        pop1 = false;
      }
      if (!pop1 && open_list2.size() == 0) {
        pop1 = true;
      }

      if (pop1) {
        Node ret = open_list1.back();
        open_list1.pop_back();
        closed_list.insert(ret.get_pos());
        return ret;
      } else {
        Node ret = open_list2.back();
        open_list2.pop_back();
        closed_list.insert(ret.get_pos());
        return ret;
      }
    }

    /*
     * open_list_empty: determines whether the open list is empty
     * RETURN: true if the open list is empty, false otherwise
     */
    bool open_list_empty() {
      return open_list1.size() == 0 && open_list2.size() == 0;
    }

    /*
     * saw_node: return true if we've seen the node at pos
     * ARGUMENTS
     * pos: the position of the node
     * RETURN: true if we've seen the node, false otherwise
     */
    bool saw_node(const NodePos &pos) {
      return seen_nodes.find(pos) != seen_nodes.end();
    }

    /*
     * get_path: gets path
     * ARGUMENTS
     * path: populated with path
     */
    void get_path(vector<std::array<int, 3>> &path) {
      int num_nodes = 1;
      NodePos pos(goalx, goaly, goalz);
      Node node = seen_nodes[pos];
      while (pos.x != startx || pos.y != starty || pos.z != startz) {
        pos = node.get_back();
        node = seen_nodes[pos];
        ++num_nodes;
      }
      path.resize(num_nodes);

      pos = NodePos(goalx, goaly, goalz);
      node = seen_nodes[pos];
      int i = num_nodes - 1;
      while (pos.x != startx || pos.y != starty || pos.z != startz) {
        path[i][0] = pos.x;
        path[i][1] = pos.y;
        path[i][2] = pos.z;
        pos = node.get_back();
        node = seen_nodes[pos];
        --i;
      }
      path[i][0] = pos.x;
      path[i][1] = pos.y;
      path[i][2] = pos.z;
    }

    const std::unordered_map<NodePos, Node, NodePosHash> &get_seen_nodes() {
      return seen_nodes;
    }
};

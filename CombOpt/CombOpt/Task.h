#pragma once

#include <vector>
#include <tuple>

typedef std::vector<char> TCHain;
typedef std::tuple<int, int, int> TCoord;

struct Direction 
  {
  enum EDir
    {
    UP,
    RIGHT,
    DOWN,
    LEFT,
    FRONT,
    BACK
    };
  };

class Task 
  {
  public:
    Task(const TCHain& i_chain) : m_chain(i_chain){}
    bool CheckDirections(const std::vector<Direction::EDir>& i_directions);
    std::vector<std::pair<TCoord, char>> GetCoords(const std::vector<Direction::EDir>& i_directions) const;
    std::vector<Direction::EDir> GRASP_Path() const;

  private:
    TCHain m_chain;
  };
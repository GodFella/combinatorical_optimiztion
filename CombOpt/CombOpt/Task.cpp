#include "Task.h"
#include <tuple>
#include <set>

std::vector<std::pair<TCoord, char>> Task::GetCoords(const std::vector<Direction::EDir>& i_directions) const
  {
  TCoord init(0, 0, 0);
  std::vector<std::pair<TCoord, char>> coords;
  size_t i = 0;
  coords.emplace_back(init, m_chain[i++]);
  
  for (const auto& dir : i_directions) {
    if (dir == Direction::EDir::UP)
      std::get<2>(init) += 1;
    else if (dir == Direction::EDir::RIGHT)
      std::get<0>(init) += 1;
    else if (dir == Direction::EDir::DOWN)
      std::get<2>(init) -= 1;
    else if (dir == Direction::EDir::LEFT)
      std::get<0>(init) -= 1;
    else if (dir == Direction::EDir::FRONT)
      std::get<1>(init) += 1;
    else if (dir == Direction::EDir::BACK)
      std::get<1>(init) -= 1;

    coords.emplace_back(init, m_chain[i++]);
    }

  return coords;
  }

bool Task::CheckDirections(const std::vector<Direction::EDir>& i_directions)
  {
  TCoord init(0, 0, 0);
  std::set<TCoord> visited;
  visited.insert(init);
  for (const auto& dir : i_directions) {
    if (dir == Direction::EDir::UP)
      std::get<2>(init) += 1;
    else if (dir == Direction::EDir::RIGHT)
      std::get<0>(init) += 1;
    else if (dir == Direction::EDir::DOWN)
      std::get<2>(init) -= 1;
    else if (dir == Direction::EDir::LEFT)
      std::get<0>(init) -= 1;
    else if (dir == Direction::EDir::FRONT)
      std::get<1>(init) += 1;
    else if (dir == Direction::EDir::BACK)
      std::get<1>(init) -= 1;
    
    if (visited.count(init))
      return false;
    else
      visited.insert(init);
    }

  return true;
  }
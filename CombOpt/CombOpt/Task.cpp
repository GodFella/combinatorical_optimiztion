#include "Task.h"
#include <tuple>
#include <set>
#include <array>
#include <chrono>

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

std::vector<Direction::EDir> Task::GRASP_Path() const {
  srand((std::chrono::system_clock::now().time_since_epoch().count(), 100));
  
  std::vector<Direction::EDir> result;
  std::vector<TCoord> coordinates;
  coordinates.emplace_back(0, 0, 0);

  std::vector<std::array<int, 3>> gaps = { { 1, 0, 0 },{ -1, 0, 0 },{ 0, 1, 0 },{ 0, -1, 0 },{ 0, 0, 1 },{ 0, 0, -1 } };
  std::set<TCoord> visited;
  visited.insert(coordinates.back());

  while (result.size() < m_chain.size() - 1) {
    
    std::vector<TCoord> cur_possible;
    for (const auto& g : gaps) {
      auto temp = coordinates.back();
      std::get<0>(temp) += g[0];
      std::get<1>(temp) += g[1];
      std::get<2>(temp) += g[2];
      if (visited.find(temp) == visited.end())
        cur_possible.push_back(temp);
      }

    if (cur_possible.empty()) {
      visited.erase(coordinates.back());
      coordinates.pop_back();
      result.pop_back();
      }
    coordinates.push_back(cur_possible[rand() % cur_possible.size()]);
    const auto& cur = coordinates.back();
    const auto& prev = coordinates[coordinates.size() - 2];
    if (std::get<0>(cur) - std::get<0>(prev) == -1)
      result.push_back(Direction::EDir::LEFT);
    else if(std::get<0>(cur) - std::get<0>(prev) == 1)
      result.push_back(Direction::EDir::RIGHT);
    else if (std::get<1>(cur) - std::get<1>(prev) == -1)
      result.push_back(Direction::EDir::BACK);
    else if (std::get<1>(cur) - std::get<1>(prev) == 1)
      result.push_back(Direction::EDir::FRONT);
    else if (std::get<2>(cur) - std::get<2>(prev) == -1)
      result.push_back(Direction::EDir::DOWN);
    else if (std::get<2>(cur) - std::get<2>(prev) == 1)
      result.push_back(Direction::EDir::UP);
    
    visited.insert(cur);
    }

  return result;
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
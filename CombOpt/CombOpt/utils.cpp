#include "utils.h"

#include <limits>
#include <array>

bool SelfIntersection(std::vector<std::vector<std::vector<char>>>& grid, 
  int shift, const std::vector<std::pair<TCoord, char>>& i_chain)
  {
  bool intersect = false;
  for (const auto& c : i_chain) {
    const auto& coord = c.first;
    if (grid[std::get<0>(coord) + shift][std::get<1>(coord) + shift][std::get<2>(coord) + shift] == 1)
      {
      intersect = true;
      break;
      }
    else
      grid[std::get<0>(coord) + shift][std::get<1>(coord) + shift][std::get<2>(coord) + shift] = 1;
    }
  
  for (const auto& c : i_chain) {
    const auto& coord = c.first;
    grid[std::get<0>(coord) + shift][std::get<1>(coord) + shift][std::get<2>(coord) + shift] = 0;
    }

  return intersect;
  }

int ComputeEnergy(std::vector<std::vector<std::vector<char>>>& grid, int shift, 
  const std::vector<std::pair<TCoord, char>>& i_chain) {

  for (const auto& c : i_chain) {
    const auto& coord = c.first;
    if(c.second == 1)
      grid[std::get<0>(coord) + shift][std::get<1>(coord) + shift][std::get<2>(coord) + shift] = 1;
    }

  int energy = 0;
  std::vector<std::array<int, 3>> gaps = { {1, 0, 0},{-1, 0, 0},{0, 1, 0},{0, -1, 0},{0, 0, 1}, {0, 0, -1} };
  
  for (const auto& c : i_chain) {
    const auto& coord = c.first;
    if (c.second == 1) {
      std::array<int, 3> temp = { 
        std::get<0>(coord) + shift, 
        std::get<1>(coord) + shift, 
        std::get<2>(coord) + shift };
      
      for (const auto& g : gaps) {
        auto cur = temp;
        cur[0] += g[0];
        cur[1] += g[1];
        cur[2] += g[2]; 
        if (cur[0] >= 0 && cur[0] < grid.size() &&
          cur[1] >= 0 && cur[1] < grid[0].size() && cur[2] >= 0 && cur[2] < grid[0][0].size())
          {
          if (grid[cur[0]][cur[1]][cur[2]] == 1)
            ++energy;
          }
        }
      grid[temp[0]][temp[1]][temp[2]] = 0;
      }
    }
  /*
  for (size_t i = 0; i < grid.size(); ++i) {
    for (size_t j = 0; j < grid[0].size(); ++j) {
      for (size_t k = 0; k < grid[0][0].size(); ++k) {
        if (grid[i][j][k] == 0)
          continue;
        for (const auto& g : gaps) {
          std::array<int, 3> cur = { i, j, k };
          cur[0] += g[0];
          cur[1] += g[1];
          cur[2] += g[2];
          if (cur[0] >= 0 && cur[0] < grid.size() &&
            cur[1] >= 0 && cur[1] < grid[0].size() && cur[2] >= 0 && cur[2] < grid[0][0].size())
            {
            if (grid[cur[0]][cur[1]][cur[2]] == 1)
              ++energy;
            }
          }
        grid[i][j][k] = 0;
        }
      }
    }*/

  return -energy;
  }

void _GetNeighbourDirections(std::vector<std::vector<Direction::EDir>>& result,
  std::vector<Direction::EDir>& dir, size_t pos, size_t radius) {
  if (radius == 0) {
    result.push_back(dir);
    return;
    }
    
  for (size_t i = pos; i < dir.size(); ++i) {
    for (int j = 0; j < 6; ++j) {
      if (j == dir[i])
        continue;
      auto temp = dir[i];
      dir[i] = static_cast<Direction::EDir>(j);
      _GetNeighbourDirections(result, dir, pos + 1, radius - 1);
      dir[i] = temp;
      }
    }
  }

std::vector<std::vector<Direction::EDir>> GetNeighbourDirections(const std::vector<Direction::EDir>& i_directions, size_t i_radius) {
  auto dir = i_directions;
  std::vector<std::vector<Direction::EDir>> result;
  _GetNeighbourDirections(result, dir, 0, i_radius);
  
  return result;
  }
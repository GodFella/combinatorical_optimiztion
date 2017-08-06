#include "utils.h"

#include <limits>
#include <array>

#include <iostream>

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

NeighbourProcessor& NeighbourProcessor::operator++() {
  ++m_counter;
  m_current[r_cur] = static_cast<Direction::EDir>((m_current[r_cur] + 1) % 6);
  if (m_current[r_cur] == 0) {
    m_current[r_cur] = m_initial[r_cur];
    ++r_cur;
    if (r_cur >= m_current.size()) {
      if (l_cur != -1) {
        m_current[l_cur] = static_cast<Direction::EDir>((m_current[l_cur] + 1) % 6);
        if (m_current[l_cur] == 0) {
          m_current[l_cur] = m_initial[l_cur];
          ++l_cur;
          if (l_cur == r_cur - 1) {
            l_cur = -1;
            r_cur = 0;
            m_current[r_cur] = static_cast<Direction::EDir>(0);
            }
          else {
            m_current[l_cur] = static_cast<Direction::EDir>(0);
            r_cur = l_cur + 1;
            m_current[r_cur] = static_cast<Direction::EDir>(0);
            }
          }
        else {
          r_cur = l_cur + 1;
          m_current[r_cur] = static_cast<Direction::EDir>(0);
          }
        }
      else {
        ++l_cur;
        m_current[l_cur] = static_cast<Direction::EDir>(0);
        r_cur = l_cur + 1;
        m_current[r_cur] = static_cast<Direction::EDir>(0);
        }
      }
    else {
      m_current[r_cur] = static_cast<Direction::EDir>(0);
      }
    }

  return *this;
  }

bool NeighbourProcessor::Finished() const{
  const int limit = ((m_initial.size() * (m_initial.size() - 1) * 36) / 2) + m_initial.size() * 6;
  return m_counter >= limit;
  }

const std::vector<Direction::EDir>& NeighbourProcessor::GetCurrentNeighbour() const {
  return m_current;
  }

void NeighbourProcessor::SetUp() {
  m_counter = 0;
  m_initial = m_current;
  }
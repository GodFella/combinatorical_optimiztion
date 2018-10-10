#include "utils.h"

#include <limits>
#include <array>
#include <string>
#include <sstream>
#include <numeric>
#include <set>

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

double ComputeEnergyPhysics(std::vector<std::vector<std::vector<char>>>& grid, int shift,
  const std::vector<std::pair<TCoord, char>>& i_chain) {

  double result = 0.0;
  for (size_t i = 0; i < i_chain.size(); ++i) {
    if (i_chain[i].second == 1) {
      for (size_t j = i + 1; j < i_chain.size(); ++j) {
        if (i_chain[j].second == 1) {
          const auto& f = i_chain[i].first;
          const auto& s = i_chain[j].first;
          const int x_diff = std::get<0>(f) - std::get<0>(s);
          const int y_diff = std::get<1>(f) - std::get<1>(s);
          const int z_diff = std::get<2>(f) - std::get<2>(s);
          const double cur_en = 1.0 / static_cast<double>(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);
          result += cur_en;
          }
        }
      }
    }
  
  return -result;
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
  
  for (size_t i = 0; i < i_chain.size(); ++i) {
    const auto& coord = i_chain[i].first;
    if (i_chain[i].second == 1) {
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
          std::array<int, 3> next = {};
          std::array<int, 3> prev = { -shift - 1, -shift - 1, -shift -1 };
          std::tie(next[0], next[1], next[2]) = i_chain[(i + 1) % i_chain.size()].first;
          if(i != 0)
            std::tie(prev[0], prev[1], prev[2]) = i_chain[(i - 1)].first;
          next[0] += shift; next[1] += shift; next[2] += shift;
          prev[0] += shift; prev[1] += shift; prev[2] += shift;
          
          if (grid[cur[0]][cur[1]][cur[2]] == 1 && cur != next && cur != prev)
            ++energy;
          }
        }
      grid[temp[0]][temp[1]][temp[2]] = 0;
      }
    }

  return -energy;
  }

double FitnessFunction(const std::vector<std::pair<TCoord, char>>& i_chain) {
  double x_min = std::numeric_limits<double>::infinity();
  double x_max = -x_min, y_min = x_min, y_max = -x_min, z_min = x_min, z_max = -x_min;

  for (const auto& c : i_chain) {
    const auto& coord = c.first;
    if (c.second == 1) {
      if (std::get<0>(coord) < z_min) {
        z_min = std::get<0>(coord);
        }
      if (std::get<0>(coord) > z_max) {
        z_max = std::get<0>(coord);
        }
      if (std::get<1>(coord) < y_min) {
        y_min = std::get<1>(coord);
        }
      if (std::get<1>(coord) > y_max) {
        y_max = std::get<1>(coord);
        }
      if (std::get<2>(coord) < x_min) {
        x_min = std::get<2>(coord);
        }
      if (std::get<2>(coord) > x_max) {
        x_max = std::get<2>(coord);
        }
      }
    }

  /*const double z_mid = (z_min + z_max) / 2.0;
  const double y_mid = (y_min + y_max) / 2.0;
  const double x_mid = (x_min + x_max) / 2.0;

  double res = 0.0;

  for (const auto& c : i_chain) {
    const auto& coord = c.first;
    if (c.second == 1) {
      double cur_z = static_cast<double>(std::get<0>(coord));
      double cur_y = static_cast<double>(std::get<1>(coord));
      double cur_x = static_cast<double>(std::get<2>(coord));
      res += (cur_z - z_mid) * (cur_z - z_mid) + (cur_y - y_mid) * (cur_y - y_mid) + (cur_x - x_mid) * (cur_x - x_mid);
      }
    }

  return res; */

  return (z_max - z_min + 1) * (y_max - y_min + 1) * (x_max - x_min + 1);
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
  return m_counter >= m_limit;
  }

int NeighbourProcessor::GetCounter() const {
  return m_counter;
  }

const std::vector<Direction::EDir>& NeighbourProcessor::GetCurrentNeighbour() const {
  return m_current;
  }

int NeighbourProcessor::GetRadiusSize() const {
  return m_limit;
  }

void NeighbourProcessor::print() const {
  std::cout << l_cur << " " << r_cur << std::endl;
  }

void NeighbourProcessor::SetUp() {
  m_counter = 0;
  m_initial = m_current;
  }

namespace {
  inline std::string serialize(const std::array<int, 3>& i_arr) {
    return std::to_string(i_arr[0]) + " " + std::to_string(i_arr[1]) + " " + std::to_string(i_arr[2]);
    }
  std::array<int, 3> deserialize(const std::string& i_serialized) {
    std::stringstream ss;
    ss << i_serialized;
    std::array<int, 3> result = { { 0, 0, 0 } };
    for (size_t j = 0; j < 3; ++j) {
      ss >> result[j];
      }

    return result;
    }

  }

namespace aco {
  
  int computeEnergyDiff(std::map<std::string, size_t>& info, const std::string& to_insert, int idx, const std::string& i_mask) {
    if (i_mask[idx] != '1') {
      return 0;
      }

    std::vector<std::array<int, 3>> gaps = { { 1, 0, 0 },{ -1, 0, 0 },{ 0, 1, 0 },{ 0, -1, 0 },{ 0, 0, 1 },{ 0, 0, -1 } };

    int diff = 0;
    for (const auto& g : gaps) {
      const auto& last = deserialize(to_insert);
      std::array<int, 3> cur = { { last[0] + g[0], last[1] + g[1], last[2] + g[2] } };
      const auto& ser_cur = serialize(cur);
      if (info.find(ser_cur) != std::end(info) && info[ser_cur] < idx - 1 && i_mask[info[ser_cur]] == '1') {
        ++diff;
        }
      }

    return diff;
    }
  
  std::vector<std::string> generatePath(const TMemory& i_memory, size_t i_size, const std::string& i_mask) {
    
    std::vector<std::string> result;
    if (i_size == 0) {
      return result;
      }
    result.push_back("0 0 0");
    std::vector<std::array<int, 3>> gaps = { { 1, 0, 0 },{ -1, 0, 0 },{ 0, 1, 0 },{ 0, -1, 0 },{ 0, 0, 1 },{ 0, 0, -1 } };
    std::map<std::string, size_t> visited_coords;
    std::vector<std::set<std::string>> visited_rotations;
    visited_rotations.emplace_back();
    visited_coords[result.front()] = 0;
    while (result.size() != i_size) {
      if (result.size() == 11) {
        int k = 0;
        k;
        }
      std::array<int, 3> last_pos = deserialize(result.back());
      std::vector<double> probs;
      std::vector<std::array<int, 3>> poss_pos;
      for (const auto& g : gaps) {
        std::array<int, 3> cur = { last_pos[0] + g[0], last_pos[1] + g[1], last_pos[2] + g[2] };
        const auto& checked = serialize(cur);
        if (visited_coords.find(checked) != std::end(visited_coords)
          || visited_rotations.back().find(checked) != std::end(visited_rotations.back())) {
          continue;
          }

        const auto& checked_pair = std::make_pair(std::to_string(result.size() - 1), serialize(g));
        double pherom_part = 0.2;
        if (i_memory.find(checked_pair) != i_memory.end()) {
          const double e = std::exp(i_memory.at(checked_pair));
          //const double pherom_part = e / (e + 1);
          pherom_part = e;
          }

        double eur_part = std::exp(computeEnergyDiff(visited_coords, checked, visited_coords.size(), i_mask));
        eur_part *= eur_part;
        probs.push_back(pherom_part * eur_part);

        poss_pos.push_back(cur);
        }

      if (poss_pos.empty()) {
        visited_coords.erase(result.back());
        visited_rotations.pop_back();
        result.pop_back();
        continue;
        }

      double sum = std::accumulate(std::begin(probs), std::end(probs), 0.0);
      probs[0] /= sum;
      for (unsigned short i = 1; i < probs.size(); ++i) {
        probs[i] /= sum;
        probs[i] += probs[i - 1];
        }

      probs.back() = 1.01;
      double rand_val = static_cast<double>(rand() % 10000) / 10000.0;
      auto idx = std::distance(
        std::begin(probs),
        std::lower_bound(std::begin(probs), std::end(probs), rand_val)
      );
      result.push_back(serialize(poss_pos[idx]));
      visited_coords[result.back()] = result.size() - 1;
      visited_rotations.back().insert(result.back());
      visited_rotations.emplace_back();
      }

    return result;
    }

  bool hasIntersections(const std::vector<std::string>& i_path) {
    std::set<std::string> visited;
    for (const auto& p : i_path) {
      if (visited.count(p)) {
        return true;
        }
      else {
        visited.insert(p);
        }
      }
    
    return false;
    }

  int computeEnergy(const std::vector<std::string>& i_path, const std::string& i_mask) {
    std::map<std::string, size_t> info;
    for (size_t i = 0; i < i_path.size(); ++i) {
      info[i_path[i]] = i;
      }

    int energy = 0;
    std::vector<std::array<int, 3>> gaps = { { 1, 0, 0 },{ -1, 0, 0 },{ 0, 1, 0 },{ 0, -1, 0 },{ 0, 0, 1 },{ 0, 0, -1 } };
    size_t idx = 0;
    for (const auto& path_comp : i_path) {
      const auto& last = deserialize(path_comp);
      if (i_mask[idx] != '1') {
        ++idx;
        continue;
        }
      for (const auto& g : gaps) {
        std::array<int, 3> cur = { {last[0] + g[0], last[1] + g[1], last[2] + g[2]} };
        const auto& ser_cur = serialize(cur);
        if (info.find(ser_cur) != std::end(info) && info[ser_cur] > idx + 1 && i_mask[info[ser_cur]] == '1') {
          ++energy;
          }
        }
      ++idx;
      }
    
    return -energy;
    }

  bool findRandNextBetter(std::vector<std::string>& better, int& better_energy, int cur_best_energy,
    size_t pos, int radius, const std::string& i_mask, std::vector<std::pair<size_t, std::array<int, 3>>>& changes) {
    if (radius == 0 || pos == better.size()) {
      if (changes.empty()) {
        return false;
        }
      auto possible_better = better;
      for (const auto& ch : changes) {
        for (size_t i = ch.first; i < better.size(); ++i) {
          auto n = deserialize(possible_better[i]);
          n[0] += ch.second[0];
          n[1] += ch.second[1];
          n[2] += ch.second[2];
          possible_better[i] = serialize(n);
          }
        }
      size_t idx_to_start = changes.back().first + 1;
      auto residue = generatePath(aco::TMemory(), better.size() - idx_to_start, i_mask.substr(idx_to_start));
      const auto& initial = idx_to_start < possible_better.size() ? 
        deserialize(possible_better[idx_to_start]) : std::array<int, 3>();
      for (size_t i = idx_to_start; i < better.size(); ++i) {
        const auto& cur_update = deserialize(residue[i - idx_to_start]);
        possible_better[i] = serialize({ cur_update[0] + initial[0], cur_update[1] + initial[1], cur_update[2] + initial[2] });
        }

      if (hasIntersections(possible_better)) {
        return false;
        }
      const auto possible_better_en = computeEnergy(possible_better, i_mask);
      if (possible_better_en < cur_best_energy) {
        better = possible_better;
        better_energy = possible_better_en;
        return true;
        }

      return false;
      }

    if (findRandNextBetter(better, better_energy, cur_best_energy, pos + 1, radius, i_mask, changes)) {
      return true;
      }

    std::vector<std::array<int, 3>> gaps = { { 1, 0, 0 },{ -1, 0, 0 },{ 0, 1, 0 },{ 0, -1, 0 },{ 0, 0, 1 },{ 0, 0, -1 } };
    for (const auto& g : gaps) {
      const auto& last = deserialize(better[pos]);
      const auto& pre_last = deserialize(better[pos - 1]);
      const std::array<int, 3> cur_rot = { last[0] - pre_last[0], last[1] - pre_last[1], last[2] - pre_last[2] };
      const std::array<int, 3> cur_change = { g[0] - cur_rot[0], g[1] - cur_rot[1], g[2] - cur_rot[2] };
      const std::array<int, 3> zero = { 0, 0, 0 };
      if (cur_change == zero) {
        continue;
        }
      changes.push_back({ pos, cur_change });
      if (findRandNextBetter(better, better_energy, cur_best_energy, pos + 1, radius - 1, i_mask, changes)) {
        return true;
        }
      changes.pop_back();
      }

    return false;
    }

  bool findNextBetter(std::vector<std::string>& better, int& better_energy, int cur_best_energy,
    size_t pos, int radius, const std::string& i_mask, std::vector<std::pair<size_t, std::array<int, 3>>>& changes) {
    if (radius == 0 || pos == better.size()) {
      auto possible_better = better;
      for (const auto& ch : changes) {
        for (size_t i = ch.first; i < better.size(); ++i) {
          auto n = deserialize(possible_better[i]);
          n[0] += ch.second[0];
          n[1] += ch.second[1];
          n[2] += ch.second[2];
          possible_better[i] = serialize(n);
          }
        }
      if (hasIntersections(possible_better)) {
        return false;
        }
      const auto possible_better_en = computeEnergy(possible_better, i_mask);
      if (possible_better_en < cur_best_energy) {
        better = possible_better;
        better_energy = possible_better_en;
        return true;
        }

      return false;
      }

    if (findNextBetter(better, better_energy, cur_best_energy, pos + 1, radius, i_mask, changes)) {
      return true;
      }

    std::vector<std::array<int, 3>> gaps = { { 1, 0, 0 },{ -1, 0, 0 },{ 0, 1, 0 },{ 0, -1, 0 },{ 0, 0, 1 },{ 0, 0, -1 } };
    for (const auto& g : gaps) {
      const auto& last = deserialize(better[pos]);
      const auto& pre_last = deserialize(better[pos - 1]);
      const std::array<int, 3> cur_rot = { last[0] - pre_last[0], last[1] - pre_last[1], last[2] - pre_last[2] };
      const std::array<int, 3> cur_change = { g[0] - cur_rot[0], g[1] - cur_rot[1], g[2] - cur_rot[2] };
      const std::array<int, 3> zero = { 0, 0, 0 };
      if (cur_change == zero) {
        continue;
        }
      changes.push_back({ pos, cur_change });
      if (findNextBetter(better, better_energy, cur_best_energy, pos + 1, radius - 1, i_mask, changes)) {
        return true;
        }
      changes.pop_back();
      }

    return false;
    }

  std::vector<std::string> dls(const std::vector<std::string>& i_initial, int radius, const std::string& i_mask) {
    int energy = computeEnergy(i_initial, i_mask);
    std::vector<std::string> better = i_initial, best = i_initial;
    int better_energy = 0;
    std::vector<std::pair<size_t, std::array<int, 3>>> changes;
    while (findNextBetter(better, better_energy, energy, 1, radius, i_mask, changes)) {
      changes.clear();
      energy = better_energy;
      best = better;
      std::cout << "   dls: " << energy << std::endl;
    }
    std::cout << "_____ dls finshed" << std::endl;
    return best;
    }

  std::vector<std::string> rdls(const std::vector<std::string>& i_initial, int radius, const std::string& i_mask) {
    int energy = computeEnergy(i_initial, i_mask);
    std::vector<std::string> better = i_initial, best = i_initial;
    int better_energy = 0;
    std::vector<std::pair<size_t, std::array<int, 3>>> changes;
    while (findRandNextBetter(better, better_energy, energy, 1, radius, i_mask, changes)) {
      changes.clear();
      energy = better_energy;
      best = better;
      //std::cout << "   dls: " << energy << std::endl;
      }
    //std::cout << "_____ dls finshed" << std::endl;
    return best;
    }

  void pheromon_expire(TMemory& io_memory, double rho) {
    for (auto& m : io_memory) {
      m.second = (1.0 - rho) * m.second;
      }
    }

  } // namespace aco
#pragma once
#include "task.h"

#include <map>
#include <vector>


bool SelfIntersection(std::vector<std::vector<std::vector<char>>>& grid, int shift, const std::vector<std::pair<TCoord, char>>& i_chain);
double ComputeEnergyPhysics(std::vector<std::vector<std::vector<char>>>& grid, int shift,
  const std::vector<std::pair<TCoord, char>>& i_chain);
int ComputeEnergy(std::vector<std::vector<std::vector<char>>>& grid, int shift, const std::vector<std::pair<TCoord, char>>& i_chain);
double FitnessFunction(const std::vector<std::pair<TCoord, char>>& i_chain);

class NeighbourProcessor {
private:
  std::vector<Direction::EDir> m_initial;
  std::vector<Direction::EDir> m_current;
  int l_cur = -1;
  int r_cur = 0;
  int m_counter = 0;
  int m_limit = -1;

public:
  NeighbourProcessor(const std::vector<Direction::EDir>& i_initial) : 
    m_initial(i_initial), m_current(i_initial), 
    m_limit(((m_initial.size() * (m_initial.size() - 1) * 36) / 2) + m_initial.size() * 6)
    {}
  NeighbourProcessor& operator++ ();
  bool Finished() const;
  int GetCounter() const;
  const std::vector<Direction::EDir>& GetCurrentNeighbour() const;
  int GetRadiusSize() const;
  void print() const;
  void SetUp();
  };

namespace aco {
  using TMemory = std::map<std::pair<std::string, std::string>, double>;
  
  std::vector<std::string> generatePath(const TMemory& i_memory, size_t i_size, const std::string& i_mask);
  std::vector<std::string> dls(const std::vector<std::string>& i_initial, int radius, const std::string& i_mask);
  std::vector<std::string> rdls(const std::vector<std::string>& i_initial, int radius, const std::string& i_mask);
  bool hasIntersections(const std::vector<std::string>& i_path);
  int computeEnergy(const std::vector<std::string>& i_path, const std::string& i_mask);
  void pheromon_expire(TMemory& io_memory, double rho);
  }
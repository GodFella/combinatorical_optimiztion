#pragma once
#include "task.h"

#include <vector>


bool SelfIntersection(std::vector<std::vector<std::vector<char>>>& grid, int shift, const std::vector<std::pair<TCoord, char>>& i_chain);
int ComputeEnergy(std::vector<std::vector<std::vector<char>>>& grid, int shift, const std::vector<std::pair<TCoord, char>>& i_chain);

class NeighbourProcessor {
private:
  std::vector<Direction::EDir> m_initial;
  std::vector<Direction::EDir> m_current;
  int l_cur = -1;
  int r_cur = 0;
  int m_counter = 0;

public:
  NeighbourProcessor(const std::vector<Direction::EDir>& i_initial) : 
    m_initial(i_initial), m_current(i_initial)
    {}
  NeighbourProcessor& operator++ ();
  bool Finished() const;
  const std::vector<Direction::EDir>& GetCurrentNeighbour() const;
  void SetUp();
  };
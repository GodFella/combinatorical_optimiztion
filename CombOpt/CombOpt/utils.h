#pragma once
#include "task.h"

#include <vector>

std::vector<std::vector<Direction::EDir>> GetNeighbourDirections(const std::vector<Direction::EDir>& i_directions, size_t i_radius);
bool SelfIntersection(std::vector<std::vector<std::vector<char>>>& grid, int shift, const std::vector<std::pair<TCoord, char>>& i_chain);
int ComputeEnergy(std::vector<std::vector<std::vector<char>>>& grid, int shift, const std::vector<std::pair<TCoord, char>>& i_chain);
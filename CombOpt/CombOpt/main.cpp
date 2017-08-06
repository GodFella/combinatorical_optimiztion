#include "Task.h"
#include "utils.h"

#include <fstream>
#include <iostream>
#include <random>
#include <algorithm>

int main()
  {
  std::ifstream fi("input.txt");
  size_t n = 0;
  fi >> n;
  std::vector<char> initial;

  for (size_t i = 0; i < n; ++i) {
    int c = 0;
    fi >> c;
    initial.push_back(c);
    }

  std::vector<std::vector<std::vector<char>>> grid(2 * n - 1,
    std::vector<std::vector<char>>(2 * n - 1,
      std::vector<char>(2 * n - 1, 0)));

  std::vector<Direction::EDir> d(n - 1);
  size_t ind = d.size() / 2;
  for (size_t i = 0; i < ind; ++i)
    d[i] = Direction::EDir::RIGHT;
  for (size_t i = ind; i < d.size(); ++i)
    d[i] = Direction::EDir::UP;

  Task problem(initial);
  NeighbourProcessor proc(d);
  int min_energy = 0;
  while (true) {
    //std::cout << "iter" << std::endl;
    auto coord0 = problem.GetCoords(d);
    min_energy = ComputeEnergy(grid, n - 1, coord0);
    auto min_dir = d;
    bool stop = true;
    while (!proc.Finished())
      {
      const auto& ne = proc.GetCurrentNeighbour();
      auto coord1 = problem.GetCoords(ne);
      if (SelfIntersection(grid, n - 1, coord1)) {
        ++proc;
        continue;
        }
        
      int energy1 = ComputeEnergy(grid, n - 1, coord1);
      if (energy1 < min_energy) {
        min_energy = energy1;
        min_dir = ne;
        stop = false;
        proc.SetUp();
        break;
        }

      ++proc;
      }

    if (stop)
      break;

    d = min_dir;
    }

for (const auto& dir : d) {
  if (dir == Direction::EDir::UP)
    std::cout << "UP" << std::endl;
  else if (dir == Direction::EDir::RIGHT)
    std::cout << "RIGHT" << std::endl;
  else if (dir == Direction::EDir::DOWN)
    std::cout << "DOWN" << std::endl;
  else if (dir == Direction::EDir::LEFT)
    std::cout << "LEFT" << std::endl;
  else if (dir == Direction::EDir::FRONT)
    std::cout << "FRONT" << std::endl;
  else if (dir == Direction::EDir::BACK)
    std::cout << "BACK" << std::endl;
  }

std::cout << min_energy << std::endl;
system("pause");
return 0;
  }
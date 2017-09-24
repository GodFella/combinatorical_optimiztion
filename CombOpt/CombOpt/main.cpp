#include "Task.h"
#include "utils.h"

#include <fstream>
#include <iostream>
#include <random>
#include <algorithm>
#include <chrono>

std::tuple<std::vector<Direction::EDir>, int> DetermLocalSearch(
  std::vector<std::vector<std::vector<char>>>& i_grid, int i_shift,
  const std::vector<Direction::EDir>& i_initial, const Task& i_problem) {
  NeighbourProcessor proc(i_initial);
  int min_energy = 0;
  auto d = i_initial;
  while (true) {
    //std::cout << "iter" << std::endl;
    auto coord0 = i_problem.GetCoords(d);
    min_energy = ComputeEnergy(i_grid, i_shift, coord0);
    auto min_dir = d;
    bool stop = true;
    while (!proc.Finished())
      {
      const auto& ne = proc.GetCurrentNeighbour();
      auto coord1 = i_problem.GetCoords(ne);
      if (SelfIntersection(i_grid, i_shift, coord1)) {
        ++proc;
        continue;
        }

      int energy1 = ComputeEnergy(i_grid, i_shift, coord1);
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

  return{ d, min_energy };
  }

std::tuple<std::vector<Direction::EDir>, int> SimulatedAnnealing(
  std::vector<std::vector<std::vector<char>>>& i_grid, int i_shift,
  const std::vector<Direction::EDir>& i_initial, const Task& i_problem) {
  
  srand(std::chrono::system_clock::now().time_since_epoch().count());
  NeighbourProcessor proc(i_initial);
  auto x = i_initial;
  auto x_res = x;
  double T = 10.0;
  const double alpha = 0.6;
  const double Tbound = 0.05;
  int x_energy = ComputeEnergy(i_grid, i_shift, i_problem.GetCoords(x));
  int min_energy = 0;

  const int switches = proc.GetRadiusSize();
  while (T > Tbound) {
    int i = 0;
    while (!proc.Finished() && i < switches) {
      const auto& y = proc.GetCurrentNeighbour();
      auto coord1 = i_problem.GetCoords(y);

      if (SelfIntersection(i_grid, i_shift, coord1)) {
        ++proc;
        continue;
        }

      int y_energy = ComputeEnergy(i_grid, i_shift, coord1);
      const int delta = y_energy - x_energy;
      const double p_thresh = static_cast<double>(rand() % 1000000) / 1000000.0;
      const double p = std::min(1.0, std::exp(-1.0 * (static_cast<double>(delta) / T)));
      if (p >= p_thresh) {
        x = y;
        ++i;
        x_energy = y_energy;
        if (x_energy < min_energy) {
          x_res = x;
          min_energy = x_energy;
          }
        proc.SetUp();
        }
      ++proc;
      }
    
    T *= alpha;
    }

  return{ x_res, min_energy };
  }

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

  Task problem(initial);
  auto&& d = problem.GRASP_Path();

  int min_energy = 0;
  //std::tie(d, min_energy) = DetermLocalSearch(grid, n - 1, d, problem);
  std::tie(d, min_energy) = SimulatedAnnealing(grid, n - 1, d, problem);

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
#include "Task.h"
#include "utils.h"

#include <fstream>
#include <iostream>
#include <random>
#include <algorithm>
#include <chrono>
#include <numeric>
#include <set>
#include <functional>
#include <string>

class PackStatistics {
public:
  double operator()(const std::vector<int>& i_pack) const {
    return std::accumulate(std::begin(i_pack), std::end(i_pack), 0) / static_cast<double>(i_pack.size());
    }
  };

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

namespace {
  struct Params {
    double T = 20.0;
    double alpha = 0.9999;
    double Tbound = 0.05;
    int pack_size = 300;
    double eps = 0.0001;
    };
  }

std::tuple<std::vector<Direction::EDir>, int> SimulatedAnnealing(
  std::vector<std::vector<std::vector<char>>>& i_grid, int i_shift,
  const std::vector<Direction::EDir>& i_initial, const Task& i_problem, const Params& i_pars) {
  
  srand((std::chrono::system_clock::now().time_since_epoch().count(), 100));
  NeighbourProcessor proc(i_initial);
  auto x = i_initial;
  auto x_res = x;
  double T = i_pars.T;
  const double alpha = i_pars.alpha;
  const double Tbound = i_pars.Tbound;
  int x_energy = ComputeEnergy(i_grid, i_shift, i_problem.GetCoords(x));
  int min_energy = 0;

  const int pack_size = i_pars.pack_size;
  std::vector<int> cur_pack;
  const double eps = i_pars.eps;

  auto less = [eps](double i_val0, double i_val1) -> bool {
    if (i_val0 < i_val1 - eps)
      return true;
    return false;
    };

  std::set<double, std::function<bool (double, double)>> statistics(less);
  std::cout << "start" << std::endl;
  while (T > Tbound) {
    std::cout << T << " " << min_energy << std::endl;
    while (true) {
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
        x_energy = y_energy;
        cur_pack.push_back(x_energy);
        if (x_energy < min_energy) {
          x_res = x;
          min_energy = x_energy;
          }
        proc.SetUp();
        }

      if (cur_pack.size() >= pack_size) {
        PackStatistics stat;
        double cur = stat(cur_pack);
        cur_pack.clear();
        if (statistics.empty() || statistics.find(cur) == statistics.end()) {
          statistics.insert(cur);
          }
        else {
          statistics.clear();
          break;
          }
        }

      ++proc;
      }
    
    T *= alpha;
    }

  std::cout << "finished" << std::endl << std::endl;

  return{ x_res, min_energy };
  }

int main()
  {
  std::ifstream fi("input.txt");
  std::ofstream fo("output.txt");
  size_t m = 0;
  fi >> m;
  
  for (size_t i = 0; i < m; ++i) {
    std::string str;
    fi >> str;
    std::vector<char> initial;

    for (size_t i = 0; i < str.size(); ++i) {
      initial.push_back(str[i] - '0');
      }

    const size_t n = str.size();
    std::vector<std::vector<std::vector<char>>> grid(2 * n - 1,
      std::vector<std::vector<char>>(2 * n - 1,
        std::vector<char>(2 * n - 1, 0)));

    Task problem(initial);
    auto&& d = problem.GRASP_Path();

    int min_energy = 0;
    Params pars;
    std::tie(d, min_energy) = DetermLocalSearch(grid, n - 1, d, problem);
    //std::tie(d, min_energy) = SimulatedAnnealing(grid, n - 1, d, problem, pars);

    fo << "For chain " << i << " : " << std::endl;
    for (const auto& dir : d) {
      if (dir == Direction::EDir::UP)
        fo << "UP" << std::endl;
      else if (dir == Direction::EDir::RIGHT)
        fo << "RIGHT" << std::endl;
      else if (dir == Direction::EDir::DOWN)
        fo << "DOWN" << std::endl;
      else if (dir == Direction::EDir::LEFT)
        fo << "LEFT" << std::endl;
      else if (dir == Direction::EDir::FRONT)
        fo << "FRONT" << std::endl;
      else if (dir == Direction::EDir::BACK)
        fo << "BACK" << std::endl;
      }

    fo << "min energy = " << min_energy << std::endl;
    fo << "________________________________________" << std::endl << std::endl;
    }
  
  system("pause");
  return 0;
  }
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
#include <map>

class PackStatistics {
public:
  double operator()(const std::vector<double>& i_pack) const {
    return std::accumulate(std::begin(i_pack), std::end(i_pack), 0.0) / static_cast<double>(i_pack.size());
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
    int pack_size = 400;
    double eps = 0.0001;
    };

  std::vector<std::function<double (int, int, double)>> GProbs =
    {
    [](int i_x, int i_y, double i_phi) -> double { 
      return 1.0 - (static_cast<double>(i_y - i_x) / (std::abs(i_x) * i_phi)); 
    }
    };

  std::vector<std::function<double (double, double)>> Gmugen = 
    {
    [](double i_x, double i_b) -> double {
      double d = (i_x - 3.0) * (i_x - 3.0) - 4 * (i_x - i_x * i_x + 1.0);
      double root_d = std::sqrt(d);
      double x_res = (3.0 - i_x - root_d) / 2.0;
      if (x_res <= i_x) {
        x_res = (3.0 - i_x + root_d) / 2.0;
        }
      return x_res;
    }
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
  double x_energy = ComputeEnergyPhysics(i_grid, i_shift, i_problem.GetCoords(x));
  double x_fit = FitnessFunction(i_problem.GetCoords(x));
  double min_energy = 0.0;
  double min_fit_sigma = std::numeric_limits<double>::infinity();
  double max_fit_sigma = -std::numeric_limits<double>::infinity();

  const int pack_size = i_pars.pack_size;
  std::vector<double> cur_pack;
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
    while (!proc.Finished()) {
      const auto& y = proc.GetCurrentNeighbour();
      auto coord1 = i_problem.GetCoords(y);

      if (SelfIntersection(i_grid, i_shift, coord1)) {
        ++proc;
        continue;
        }

      double y_energy = ComputeEnergyPhysics(i_grid, i_shift, coord1);
      double y_fit = FitnessFunction(coord1);
      const double delta = y_energy - x_energy;
      double fit_sigma = 1.0;
        /*y_fit < x_fit ? y_fit / (x_fit * std::sqrt(x_res.size())) : (y_fit * std::sqrt(x_res.size())) / x_fit;
      fit_sigma *= std::sqrt(x_res.size());*/
      
      min_fit_sigma = std::min(min_fit_sigma, fit_sigma);
      max_fit_sigma = std::max(max_fit_sigma, fit_sigma);
      //std::cout << fit_sigma << std::endl;
      const double p_thresh = static_cast<double>(rand() % 1000000) / 1000000.0;
      const double p = std::min(1.0, std::exp(-1.0 * (static_cast<double>(delta) / T)) / fit_sigma);
      if (p >= p_thresh) {
        x = y;
        x_energy = y_energy;
        x_fit = y_fit;
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

/*std::tuple<std::vector<Direction::EDir>, int> GAlgo(
  std::vector<std::vector<std::vector<char>>>& i_grid, int i_shift,
  const std::vector<Direction::EDir>& i_initial, const Task& i_problem, const Params& i_pars) {
  auto x = i_initial;
  auto x_res = x;
  int x_energy = ComputeEnergy(i_grid, i_shift, i_problem.GetCoords(x));
  int min_energy = 0;
  NeighbourProcessor proc(i_initial);

  double mu = 0.0;
  
  std::vector<double> golden;
  golden.push_back(0.0);
  for (int i = 0; i < 10; ++i)
    golden.push_back(Gmugen[0](golden.back(), 0.1));
  const int magic_number = 8000;
  const int number_per_segment = magic_number / 10;

  std::vector<double> mus;
  for (int i = 0; i < golden.size() - 1; ++i) {
    const double delta = (golden[i + 1] - golden[i]) / static_cast<double>(number_per_segment);
    mus.push_back(golden[i]);
    for (int j = 0; j < number_per_segment; ++j) {
      mus.push_back(mus.back() + delta);
      }
    }

  std::vector<int> cur_pack;
  const double eps = i_pars.eps;

  auto less = [eps](double i_val0, double i_val1) -> bool {
    if (i_val0 < i_val1 - eps)
      return true;
    return false;
    };
  std::set<double, std::function<bool(double, double)>> statistics(less);

  std::cout << "G start" << std::endl;
  const double mu_eps = 0.00001;
  for (const auto& mu : mus) {
    std::cout << mu << " " << proc.GetCounter() << " " << min_energy << std::endl;
    while (!proc.Finished()) {
      const auto& y = proc.GetCurrentNeighbour();
      auto coord1 = i_problem.GetCoords(y);

      if (SelfIntersection(i_grid, i_shift, coord1)) {
        ++proc;
        continue;
        }

      int y_energy = ComputeEnergy(i_grid, i_shift, coord1);
      const double p = std::min(1.0, GProbs[0](x_energy, y_energy, 1.0));
      const double p_thresh = mu + (static_cast<double>(rand() % 1000000) / 1000000.0) * (1.0 - mu);

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

      if (cur_pack.size() >= i_pars.pack_size) {
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

    //mu = Gmugen[0](mu, 0.1);
    }

  std::cout << "G finished" << std::endl;
  return{ x_res, min_energy };
  }*/

void Compare()
  {
  std::ifstream fi("D:/CppProjects/projects/Git/CombOpt/CombOpt/Debug/input.txt");
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
    auto d = problem.GRASP_Path();

    int min_energy = 0;
    Params pars;
    //std::tie(d, min_energy) = DetermLocalSearch(grid, n - 1, d, problem);
    std::cout << "initial energy was: " << ComputeEnergyPhysics(grid, n - 1, problem.GetCoords(d)) << std::endl;
    std::tie(d, min_energy) = SimulatedAnnealing(grid, n - 1, d, problem, pars);
    //std::tie(d, min_energy) = GAlgo(grid, n - 1, d, problem, pars);

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

    fo << "min energy phyics= " << min_energy << ", min_energy_classical= " 
      << ComputeEnergy(grid, n - 1, problem.GetCoords(d)) << std::endl;
    fo << "________________________________________" << std::endl << std::endl;
    }
  
  }

void Train() {
  std::ifstream fi("Train.txt");
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
    auto d = problem.GRASP_Path();

    double alph_add = 0.9;
    int min_energy = 0;
    Params arg_min;
    fo << "initial energy was: " << ComputeEnergy(grid, n - 1, problem.GetCoords(d)) << std::endl;
    for (double cur_alph = 0.9; cur_alph < 0.99999; cur_alph += alph_add) {
      alph_add *= 0.1;
      for (int pack = 50; pack <= 400; pack += 50) {
        for (double ep = 0.1; ep > 0.00001; ep *= 0.1) {
          Params params;
          params.alpha = cur_alph;
          params.pack_size = pack;
          params.eps = ep;
          int cur_energy = 0;
          fo << "trying " << params.alpha << " " << params.pack_size << " " << params.eps << std::endl;
          std::tie(d, cur_energy) = SimulatedAnnealing(grid, n - 1, d, problem, params);
          fo << "      result: " << cur_energy << std::endl;
          if (cur_energy < min_energy) {
            min_energy = cur_energy;
            arg_min = params;
            fo << "found new better params: " << arg_min.alpha << " " << arg_min.pack_size << " " << arg_min.eps <<
              "with energy " << min_energy << std::endl;
            }
          }
        }
      }

    fo << std::endl;
    fo << std::endl;
    fo << "so the best params: " << arg_min.alpha << " " << arg_min.pack_size << " " << arg_min.eps <<
      "with energy " << min_energy << std::endl;

    }
  }

struct ACOParams {
  size_t generation_am_;
  size_t ants_in_generation_;
  std::function<double(double)> pheremon_aspiration_;
  };

std::tuple<std::vector<Direction::EDir>, int> ACO(const Task& i_problem, const ACOParams& i_aco) {
  aco::TMemory memory;
  int min_energy = 0;
  std::vector<std::string> min_path;

  srand((std::chrono::system_clock::now().time_since_epoch().count(), 100));
  for (size_t g_idx = 0; g_idx < i_aco.generation_am_; ++g_idx) {
    for (size_t a_idx = 0; a_idx < i_aco.ants_in_generation_; ++a_idx) {
      aco::TMemory local_memory;
      const auto& path = aco::generatePath(memory, i_problem.getSize(), local_memory);
      if (aco::hasIntersections(path)) {
        continue;
        }
      for (const auto& l_m : local_memory) {
          memory[l_m.first] += l_m.second;
        }
      int cur_energy = aco::computeEnergy(path);
      if (cur_energy < min_energy) {
        min_energy = cur_energy;
        min_path = path;
        }
      }
    aco::pheromon_expire(memory);
    }
  }

int main() {
  //Train();
  Compare();

  system("pause");
  return 0;
  }
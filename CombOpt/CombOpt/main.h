#pragma once

#include "Task.h"
#include <fstream>
#include <random>

int main()
  {
  std::ifstream fi("input.txt");
  size_t n = 0;
  fi >> n;
  std::vector<char> initial;
  
  for (size_t i = 0; i < n; ++i) {
    char c;
    fi >> c;
    initial.push_back(c);
    }

  srand(0);
  std::vector<Direction::EDir> d(n - 1, rand() % 6);
  Task problem(initial);
  while (true) {

    }
  return 0;
  }
/*
 * @Author: 董泰宏 2396203400@qq.com
 * @Date: 2022-12-07 12:10:57
 * @LastEditors: 董泰宏 2396203400@qq.com
 * @LastEditTime: 2023-04-18 12:56:45
 * @FilePath: /SplineSmoother/src/main.cpp
 * @Description:
 * Copyright (c) 2022 by 董泰宏 email: 2396203400@qq.com, All Rights Reserved.
 */
#include "SplineSmoother.hpp"

namespace plt = matplotlibcpp;
int main() {
  const char* filename = "../source/path.csv";
  SplineSmoother spline(filename);
  return 0;
}
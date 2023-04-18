/*
 * @Author: 董泰宏 2396203400@qq.com
 * @Date: 2022-12-07 12:10:33
 * @LastEditors: 董泰宏 2396203400@qq.com
 * @LastEditTime: 2023-04-18 12:23:37
 * @FilePath: /SplineSmoother/src/SplineSmoother.cpp
 * @Description:输入离散的全局路径点(要求每两点之间间隔不高于2.5m)，得到平滑后的全局参考线
 * Copyright (c) 2022 by 董泰宏 email: 2396203400@qq.com, All Rights Reserved.
 */
#include "SplineSmoother.hpp"

/**
 * @description: 构造函数
 * @param {char*} filename
 * @return {*}
 */
SplineSmoother::SplineSmoother(const char* filename) {
  GetSourceData(filename);
  Q.emplace_back(CostFunction());
  ConstriantFunction();
  Solution();
  ParametricToCartesian();
}

/**
 * @description: 采样获取原始路径点以及鞍点,每5米取1个input_point,
 * 每2.5米取1个anchor_point
 *               TODO: 测试不同间隔下算法的表现
 * @param {char*} filename
 * @return {*}
 */
void SplineSmoother::GetSourceData(const char* filename) {
  ifstream readCSV(filename, ios::in);
  string line;
  int input_index = 0;
  int anchor_index = 0;
  int line_number = 0;
  double distance = 0;
  while (getline(readCSV, line)) {
    stringstream ss(line);
    string s;
    vector<string> temp_line;
    pair<double, double> temp_point;
    while (getline(ss, s, ',')) {
      temp_line.emplace_back(s);
    }
    temp_point.first = stod(temp_line[0]);
    temp_point.second = stod(temp_line[1]);
    //将起点放入
    if (line_number == 0) {
      input_points.emplace_back(temp_point);
    }

    //检查路径点是鞍点或者input_points?
    distance =
        sqrt(pow(temp_point.first - input_points[input_index].first, 2) +
             pow(temp_point.second - input_points[input_index].second, 2));

    //每隔5米左右采一个参考点（间隔可以自己去调，比如1米就把这里判断条件改成1左右）
    if (distance > 4.7 && distance < 5.3) {
      input_points.emplace_back(temp_point);
      input_index++;
    }
    //在两个参考点的中间也就是2.5米左右采一个anchor point
    if (distance > 2.2 && distance < 2.8) {
      anchor_points.emplace_back(temp_point);
      anchor_index++;
    }
    line_number++;
  }
  this->number = input_points.size();
  //初始化起点和终点的朝向
  this->heading_start = (input_points[1].second - input_points[0].second) /
                        (input_points[1].first - input_points[0].first);
  this->heading_end =
      (input_points[number - 1].second - input_points[number - 2].second) /
      (input_points[number - 1].first - input_points[number - 2].first);
}

/**
 * @description: 目标函数的Q矩阵，包含正则化惩罚项
 *               目标变量: [a10,a11,a12,a13,a14,a15,...,an0,an1,an2,an3,an4,an5,
 *                              b10,b11,b12,b13,b14,b15,...,bn0,bn1,bn2,bn3,bn4,bn5]
 * 12xn
 * @return QP-the cost matrix
 */
Eigen::SparseMatrix<double> SplineSmoother::CostFunction() {
  //目标函数的q矩阵--舒适项
  Eigen::SparseMatrix<double> q(number * 12, number * 12);
  for (int i = 0; i < number * 2; i++) {
    q.insert(6 * i + 3, 6 * i + 3) = 36;
    q.insert(6 * i + 3, 6 * i + 4) = 72;
    q.insert(6 * i + 3, 6 * i + 5) = 120;
    q.insert(6 * i + 4, 6 * i + 3) = 72;
    q.insert(6 * i + 4, 6 * i + 4) = 192;
    q.insert(6 * i + 4, 6 * i + 5) = 360;
    q.insert(6 * i + 5, 6 * i + 3) = 120;
    q.insert(6 * i + 5, 6 * i + 4) = 360;
    q.insert(6 * i + 5, 6 * i + 5) = 720;
  }
  //正则化惩罚项L矩阵
  Eigen::SparseMatrix<double> L(number * 12, number * 12);
  for (int i = 0; i < number * 12; i++) {
    L.insert(i, i) = 1.0e-5;
  }
  //二次矩阵
  Eigen::SparseMatrix<double> QP(number * 12, number * 12);
  QP = q + L;

  //一次矩阵, 这里没有，全部设置为0
  Eigen::VectorXd grad = Eigen::VectorXd::Zero(number * 12);
  this->gradient.emplace_back(grad);
  return QP;
}

/**
 * @description: 构建A矩阵以及约束边界
 * @return {*}
 */
void SplineSmoother::ConstriantFunction() {
  //光滑性约束的约束：假设输入6个点，则有6段（终点后面还有一段），则有5个连接点，x对应4*5，y对应4*5,
  //则总共4*5*2
  int smooth_constriant_number = 4 * (number - 1) * 2;
  int position_constriant_number = (number - 1) * 2;
  int start_end_constriant_number = 2;
  this->total_constriant_number = smooth_constriant_number +
                                  position_constriant_number +
                                  start_end_constriant_number;

  Eigen::SparseMatrix<double> a(total_constriant_number, number * 12);

  Eigen::VectorXd LowerBound(total_constriant_number);
  Eigen::VectorXd UpperBound(total_constriant_number);

  Eigen::VectorXd T0(12);     // 0阶导数向量
  Eigen::VectorXd T0_05(12);  // 0阶导数向量参数0.5时
  Eigen::VectorXd T1(12);     // 1阶导数向量
  Eigen::VectorXd T2(12);     // 2阶导数向量
  Eigen::VectorXd T3(12);     // 3阶导数向量
  T0 << 1, 1, 1, 1, 1, 1, -1, 0, 0, 0, 0, 0;
  T1 << 0, 1, 2, 3, 4, 5, 0, -1, 0, 0, 0, 0;
  T2 << 0, 0, 2, 6, 12, 20, 0, 0, -2, 0, 0, 0;
  T3 << 0, 0, 0, 6, 24, 60, 0, 0, 0, -6, 0, 0;
  T0_05 << 1, 0.5, pow(0.5, 2), pow(0.5, 3), pow(0.5, 4), pow(0.5, 5), 0, 0, 0,
      0, 0, 0;

  //光滑性约束 4 * (number - 1) * 2
  for (int i = 0; i < number - 1; i++) {
    // fi(ti) = fi+1(t0) && gi(ti) = gi+1(t0)
    a.insert(i, 6 * i + 0) = T0[0];
    a.insert(i, 6 * i + 1) = T0[1];
    a.insert(i, 6 * i + 2) = T0[2];
    a.insert(i, 6 * i + 3) = T0[3];
    a.insert(i, 6 * i + 4) = T0[4];
    a.insert(i, 6 * i + 5) = T0[5];
    a.insert(i, 6 * i + 6) = T0[6];
    LowerBound(i) = 0;
    UpperBound(i) = 0;
    a.insert(i + total_constriant_number / 2, 6 * i + 6 * number + 0) = T0[0];
    a.insert(i + total_constriant_number / 2, 6 * i + 6 * number + 1) = T0[1];
    a.insert(i + total_constriant_number / 2, 6 * i + 6 * number + 2) = T0[2];
    a.insert(i + total_constriant_number / 2, 6 * i + 6 * number + 3) = T0[3];
    a.insert(i + total_constriant_number / 2, 6 * i + 6 * number + 4) = T0[4];
    a.insert(i + total_constriant_number / 2, 6 * i + 6 * number + 5) = T0[5];
    a.insert(i + total_constriant_number / 2, 6 * i + 6 * number + 6) = T0[6];
    LowerBound(i + total_constriant_number / 2) = 0;
    UpperBound(i + total_constriant_number / 2) = 0;
    // fi'(ti) = fi+1'(t0) && gi'(ti) = gi+1'(t0)
    a.insert(i + number - 1, 6 * i + 0) = T1[0];
    a.insert(i + number - 1, 6 * i + 1) = T1[1];
    a.insert(i + number - 1, 6 * i + 2) = T1[2];
    a.insert(i + number - 1, 6 * i + 3) = T1[3];
    a.insert(i + number - 1, 6 * i + 4) = T1[4];
    a.insert(i + number - 1, 6 * i + 5) = T1[5];
    a.insert(i + number - 1, 6 * i + 7) = T1[7];
    LowerBound(i + number - 1) = 0;
    UpperBound(i + number - 1) = 0;
    a.insert(i + total_constriant_number / 2 + number - 1,
             6 * i + 6 * number + 0) = T1[0];
    a.insert(i + total_constriant_number / 2 + number - 1,
             6 * i + 6 * number + 1) = T1[1];
    a.insert(i + total_constriant_number / 2 + number - 1,
             6 * i + 6 * number + 2) = T1[2];
    a.insert(i + total_constriant_number / 2 + number - 1,
             6 * i + 6 * number + 3) = T1[3];
    a.insert(i + total_constriant_number / 2 + number - 1,
             6 * i + 6 * number + 4) = T1[4];
    a.insert(i + total_constriant_number / 2 + number - 1,
             6 * i + 6 * number + 5) = T1[5];
    a.insert(i + total_constriant_number / 2 + number - 1,
             6 * i + 6 * number + 7) = T1[7];
    LowerBound(i + total_constriant_number / 2 + number - 1) = 0;
    UpperBound(i + total_constriant_number / 2 + number - 1) = 0;

    // fi''(ti) = fi+1''(t0) && gi''(ti) = gi+1''(t0)
    a.insert(i + (number - 1) * 2, 6 * i + 0) = T2[0];
    a.insert(i + (number - 1) * 2, 6 * i + 1) = T2[1];
    a.insert(i + (number - 1) * 2, 6 * i + 2) = T2[2];
    a.insert(i + (number - 1) * 2, 6 * i + 3) = T2[3];
    a.insert(i + (number - 1) * 2, 6 * i + 4) = T2[4];
    a.insert(i + (number - 1) * 2, 6 * i + 5) = T2[5];
    a.insert(i + (number - 1) * 2, 6 * i + 8) = T2[8];
    LowerBound(i + (number - 1) * 2) = 0;
    UpperBound(i + (number - 1) * 2) = 0;
    a.insert(i + total_constriant_number / 2 + (number - 1) * 2,
             6 * i + 6 * number + 0) = T2[0];
    a.insert(i + total_constriant_number / 2 + (number - 1) * 2,
             6 * i + 6 * number + 1) = T2[1];
    a.insert(i + total_constriant_number / 2 + (number - 1) * 2,
             6 * i + 6 * number + 2) = T2[2];
    a.insert(i + total_constriant_number / 2 + (number - 1) * 2,
             6 * i + 6 * number + 3) = T2[3];
    a.insert(i + total_constriant_number / 2 + (number - 1) * 2,
             6 * i + 6 * number + 4) = T2[4];
    a.insert(i + total_constriant_number / 2 + (number - 1) * 2,
             6 * i + 6 * number + 5) = T2[5];
    a.insert(i + total_constriant_number / 2 + (number - 1) * 2,
             6 * i + 6 * number + 8) = T2[8];
    LowerBound(i + total_constriant_number / 2 + (number - 1) * 2) = 0;
    UpperBound(i + total_constriant_number / 2 + (number - 1) * 2) = 0;

    // fi'''(ti) = fi+1'''(t0) && gi'''(ti) = gi+1'''(t0)
    a.insert(i + (number - 1) * 3, 6 * i + 0) = T3[0];
    a.insert(i + (number - 1) * 3, 6 * i + 1) = T3[1];
    a.insert(i + (number - 1) * 3, 6 * i + 2) = T3[2];
    a.insert(i + (number - 1) * 3, 6 * i + 3) = T3[3];
    a.insert(i + (number - 1) * 3, 6 * i + 4) = T3[4];
    a.insert(i + (number - 1) * 3, 6 * i + 5) = T3[5];
    a.insert(i + (number - 1) * 3, 6 * i + 9) = T3[9];
    LowerBound(i + (number - 1) * 3) = 0;
    UpperBound(i + (number - 1) * 3) = 0;
    a.insert(i + total_constriant_number / 2 + (number - 1) * 3,
             6 * i + 6 * number + 0) = T2[0];
    a.insert(i + total_constriant_number / 2 + (number - 1) * 3,
             6 * i + 6 * number + 1) = T2[1];
    a.insert(i + total_constriant_number / 2 + (number - 1) * 3,
             6 * i + 6 * number + 2) = T2[2];
    a.insert(i + total_constriant_number / 2 + (number - 1) * 3,
             6 * i + 6 * number + 3) = T2[3];
    a.insert(i + total_constriant_number / 2 + (number - 1) * 3,
             6 * i + 6 * number + 4) = T2[4];
    a.insert(i + total_constriant_number / 2 + (number - 1) * 3,
             6 * i + 6 * number + 5) = T2[5];
    a.insert(i + total_constriant_number / 2 + (number - 1) * 3,
             6 * i + 6 * number + 9) = T2[9];
    LowerBound(i + total_constriant_number / 2 + (number - 1) * 3) = 0;
    UpperBound(i + total_constriant_number / 2 + (number - 1) * 3) = 0;
  }

  //起止点约束
  a.insert(smooth_constriant_number / 2, 1) = -heading_start;
  a.insert(smooth_constriant_number / 2, 6 * number + 1) = 1;
  LowerBound(smooth_constriant_number / 2) = 0;
  UpperBound(smooth_constriant_number / 2) = 0;
  a.insert((smooth_constriant_number + total_constriant_number) / 2,
           6 * (number - 1) + 1) = -heading_end * 1;
  a.insert((smooth_constriant_number + total_constriant_number) / 2,
           6 * (number - 1) + 2) = -heading_end * 2;
  a.insert((smooth_constriant_number + total_constriant_number) / 2,
           6 * (number - 1) + 3) = -heading_end * 3;
  a.insert((smooth_constriant_number + total_constriant_number) / 2,
           6 * (number - 1) + 4) = -heading_end * 4;
  a.insert((smooth_constriant_number + total_constriant_number) / 2,
           6 * (number - 1) + 5) = -heading_end * 5;
  a.insert((smooth_constriant_number + total_constriant_number) / 2,
           6 * (2 * number - 1) + 1) = 1;
  a.insert((smooth_constriant_number + total_constriant_number) / 2,
           6 * (2 * number - 1) + 2) = 2;
  a.insert((smooth_constriant_number + total_constriant_number) / 2,
           6 * (2 * number - 1) + 3) = 3;
  a.insert((smooth_constriant_number + total_constriant_number) / 2,
           6 * (2 * number - 1) + 4) = 4;
  a.insert((smooth_constriant_number + total_constriant_number) / 2,
           6 * (2 * number - 1) + 5) = 5;
  LowerBound((smooth_constriant_number + total_constriant_number) / 2) = 0;
  UpperBound((smooth_constriant_number + total_constriant_number) / 2) = 0;
  //位置约束
  for (int i = 1; i <= number - 1; i++) {
    // x boundary && y boundary
    a.insert(smooth_constriant_number / 2 + i, 6 * (i - 1) + 0) = T0_05[0];
    a.insert(smooth_constriant_number / 2 + i, 6 * (i - 1) + 1) = T0_05[1];
    a.insert(smooth_constriant_number / 2 + i, 6 * (i - 1) + 2) = T0_05[2];
    a.insert(smooth_constriant_number / 2 + i, 6 * (i - 1) + 3) = T0_05[3];
    a.insert(smooth_constriant_number / 2 + i, 6 * (i - 1) + 4) = T0_05[4];
    a.insert(smooth_constriant_number / 2 + i, 6 * (i - 1) + 5) = T0_05[5];
    LowerBound(smooth_constriant_number / 2 + i) =
        anchor_points[i - 1].first - boundary;
    UpperBound(smooth_constriant_number / 2 + i) =
        anchor_points[i - 1].first + boundary;
    a.insert((smooth_constriant_number + total_constriant_number) / 2 + i,
             6 * (i - 1) + 6 * number + 0) = T0_05[0];
    a.insert((smooth_constriant_number + total_constriant_number) / 2 + i,
             6 * (i - 1) + 6 * number + 1) = T0_05[1];
    a.insert((smooth_constriant_number + total_constriant_number) / 2 + i,
             6 * (i - 1) + 6 * number + 2) = T0_05[2];
    a.insert((smooth_constriant_number + total_constriant_number) / 2 + i,
             6 * (i - 1) + 6 * number + 3) = T0_05[3];
    a.insert((smooth_constriant_number + total_constriant_number) / 2 + i,
             6 * (i - 1) + 6 * number + 4) = T0_05[4];
    a.insert((smooth_constriant_number + total_constriant_number) / 2 + i,
             6 * (i - 1) + 6 * number + 5) = T0_05[5];
    LowerBound((smooth_constriant_number + total_constriant_number) / 2 + i) =
        anchor_points[i - 1].second - boundary;
    UpperBound((smooth_constriant_number + total_constriant_number) / 2 + i) =
        anchor_points[i - 1].second + boundary;
  }
  this->A.emplace_back(a);
  this->bound.emplace_back(LowerBound);
  this->bound.emplace_back(UpperBound);
}

/**
 * @description:
 * @return {*}
 */
bool SplineSmoother::Solution() {
  OsqpEigen::Solver solver;
  solver.settings()->setVerbosity(false);
  solver.settings()->setWarmStart(true);
  solver.data()->setNumberOfVariables(number * 12);
  solver.data()->setNumberOfConstraints(total_constriant_number);

  if (!solver.data()->setHessianMatrix(Q[0])) return false;
  if (!solver.data()->setGradient(gradient[0])) return false;
  if (!solver.data()->setLinearConstraintsMatrix(A[0])) return false;
  if (!solver.data()->setLowerBound(bound[0])) return false;
  if (!solver.data()->setUpperBound(bound[1])) return false;
  if (!solver.initSolver()) return false;

  Eigen::VectorXd Solution;
  if (!solver.solve()) return false;
  Solution = solver.getSolution();
  smooth_parametric.emplace_back(Solution);
  return true;
}

/**
 * @description: 将参数方程转换为笛卡尔坐标系点列，并且画出实时图像
 * @return {*}
 */
void SplineSmoother::ParametricToCartesian() {
  double x, y;
  pair<double, double> temp_point;
  vector<double> X;
  vector<double> Y;
  vector<double> X1;
  vector<double> Y1;
  vector<double> X2;
  vector<double> Y2;
  for (auto& p : anchor_points) {
    X1.emplace_back(p.first);
    Y1.emplace_back(p.second);
  }
  for (auto& p : input_points) {
    X2.emplace_back(p.first);
    Y2.emplace_back(p.second);
  }
  for (int i = 0; i < number * 6; i += 6) {
    for (double j = 0; j <= 1; j += 0.01) {
      x = smooth_parametric[0][i] + smooth_parametric[0][i + 1] * j +
          smooth_parametric[0][i + 2] * pow(j, 2) +
          smooth_parametric[0][i + 3] * pow(j, 3) +
          smooth_parametric[0][i + 4] * pow(j, 4) +
          smooth_parametric[0][i + 5] * pow(j, 5);
      y = smooth_parametric[0][i + 6 * number] +
          smooth_parametric[0][i + 6 * number + 1] * j +
          smooth_parametric[0][i + 6 * number + 2] * pow(j, 2) +
          smooth_parametric[0][i + 6 * number + 3] * pow(j, 3) +
          smooth_parametric[0][i + 6 * number + 4] * pow(j, 4) +
          smooth_parametric[0][i + 6 * number + 5] * pow(j, 5);
      temp_point.first = x;
      temp_point.second = y;
      smooth_points.emplace_back(temp_point);
      X.push_back(x);
      Y.push_back(y);
      plt::clf();
      std::map<std::string, std::string> keywords1;
      keywords1.insert(pair<string, string>("color", "r"));
      keywords1.insert(pair<string, string>("label", "Smooth Path"));
      plt::plot(X, Y, keywords1);
      std::map<std::string, std::string> keywords2;
      keywords2.insert(pair<string, string>("marker", "^"));
      keywords2.insert(pair<string, string>("label", "Path Point"));
      plt::scatter(X1, Y1, 55.0, keywords2);
      std::map<std::string, std::string> keywords3;
      keywords3.insert(pair<string, string>("marker", "p"));
      keywords3.insert(pair<string, string>("label", "Anchor Point"));
      plt::scatter(X2, Y2, 35.0, keywords3);
      plt::legend();
      plt::pause(0.001);
    }
  }
}
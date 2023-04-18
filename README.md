# SplineSmoother
基于osqp复现apollo的分段spline参考线平滑算法 
使用matplotlib完成过程的可视化
## 依赖
--osqp-eigen
求解二次规划的c++库
```shell
安装参考：https://blog.csdn.net/zk1639712996/article/details/126568314
```
--matplotlib-cpp
    这个是基础的画图库，是很基本很重要经常会用到的画图库，安装也较为简单
```shell
pip3 install numpy
pip3 install matplotlib
git clone https://github.com/lava/matplotlib-cpp
CMakeLists.txt里面需要添加python相关的依赖，具体参考我的CMakeLists.txt
```
## 编译运行
```shell
git clone https://github.com/dongtaihong/SplineSmoother.git
cd SplineSmoother
mkdir build && cd build
cmake .. && make
./SplineSmoother
```
## 效果
说明：算法原理及如何调试请参考pdf文件
![image](./result_image/spline.png)

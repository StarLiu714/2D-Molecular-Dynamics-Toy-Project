# 2D-Molecular-Dynamics-Toy-Project

<center class="half">
  <image src= "https://github.com/StarLiu714/2D-Molecular-Dynamics-Toy-Project/blob/20230408/Lab-Report/assets/Hydro-20K-900.jpg" width="600"\>
</center>


## Quick guide

[Project Codes](https://github.com/StarLiu714/2D-Molecular-Dynamics-Toy-Project/tree/20230408/Proj)  

[Lab Report](https://github.com/StarLiu714/2D-Molecular-Dynamics-Toy-Project/blob/20230408/Lab-Report/Lab%20report.md)   

[Preview in Jupyter Notebook](https://github.com/StarLiu714/2D-Molecular-Dynamics-Toy-Project/blob/20230408/Notebook_MD_2d_Argon_Proj.ipynb)   

The output results are stored in the ./plot/ folder by default, and will be created if it does not exist


## Environment & Packages

`Python` 3.7+; 

Packages:  numpy, numba, os, tqdm, matplotlib


## Default Parameters
1. Argon atoms 
2. 20\*20 block
3. Desired temperature follows the list: 5~100 Kelvin, every 5 Kelvin interval outputs a series of results.
4. Thermostat is of kappa: $$\kappa = (\frac{T_{desired}}{T_{now}})^{\frac{1}{2}} $$ with a tolerance of 2e-4






# 二维分子动力学模拟

## 快速指南

[项目代码](https://github.com/StarLiu714/2D-Molecular-Dynamics-Toy-Project/tree/20230408/Proj)  

[实验报告](https://github.com/StarLiu714/2D-Molecular-Dynamics-Toy-Project/blob/20230408/Lab-Report/Lab%20report.md)  

[逐条结果可视化-ipynb](https://github.com/StarLiu714/2D-Molecular-Dynamics-Toy-Project/blob/20230408/Notebook_MD_2d_Argon_Proj.ipynb)  

输出结果默认存放于./plot/文件夹下，若不存在将创建

## 环境要求

Python 3.7及以上版本; 

相关包:  numpy, numba, os, tqdm, matplotlib

## 默认参数
1. 纯组分氩气单原子分子
2. 正方形2D块体，长\*宽=20\*20
3. 预设温度范围为5~100开尔文；每间隔5开尔文，输出一组结果
4. 恒温器检测环境温度波动，当偏离设定温度2e-4时自动回调。 
  $$\kappa = (\frac{T_{desired}}{T_{now}})^{\frac{1}{2}} $$

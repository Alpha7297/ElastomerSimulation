对于单个四面体

假设初始位置$x_i$,形变后位置$y_i$

为了消除平动可以取差值，以下$i=1,2,3$

$$
X'_i=x_i-x_0;Y_i=y_i-y_0
$$

由线性形变假设

$$
F x_i+b=y_i\\
F X_i=Y_i
$$

记

$$
R_s=(X_1,X_2,X_3)\\
R_m=(Y_1,Y_2,Y_3)
$$

此时

$$
F=R_mR_s^{-1}
$$

得到形变向量之后，假设应变张量可以写为

$$
\epsilon=\frac{1}{2}(F+F^T)-I
$$

由自由能公式，总能量为

$$
E=V_0\cdot(\mu \mathbf{tr}(\epsilon^T \epsilon)+\frac{\lambda}{2}(\mathbf{tr}\epsilon)^2)
$$

下一步进行求导，由链式法则

$$
f_i=\frac{\partial E}{\partial \epsilon}\frac{\partial \epsilon}{\partial F}\frac{\partial F}{\partial R_m}\frac{\partial R_m}{\partial y_i}
$$

记应力

$$
P=2\mu\epsilon+\lambda\mathbf{tr}\epsilon
$$

可以得到

$$
(f_1,f_2,f_3)=PR_s^{-T}
$$

而$f_0=-f_1-f_2-f_3$

下一步，考虑约束和外力，每个点质量一定，直接用牛顿第二定律即可

为了实现隐式欧拉，需要对力再对位置求一阶导。先对应力张量$P$对形变梯度$F$求导。由

$$
\epsilon=\frac{1}{2}(F+F^T)-I
$$

和

$$
P=2\mu\epsilon+\lambda\mathbf{tr}(\epsilon)I
$$

可写为

$$
P=\mu(F+F^T-2I)+\lambda(\mathbf{tr}(F)-3)I
$$

因此对其取微分有

$$
dP=\mu(dF+dF^T)+\lambda\mathbf{tr}(dF)I
$$

也就是四阶张量

$$
C_{ijkl}=\frac{\partial P_{ij}}{\partial F_{kl}}
=\mu(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk})+\lambda\delta_{ij}\delta_{kl}
$$

记

$$
B=R_s^{-1}
$$

四面体的形函数梯度可由$B$表示：

$$
g_1=(B_{1,0},B_{1,1},B_{1,2})^T\\
g_2=(B_{2,0},B_{2,1},B_{2,2})^T\\
g_3=(B_{3,0},B_{3,1},B_{3,2})^T\\
g_0=-(g_1+g_2+g_3)
$$

此时对任意一个顶点$b$的位置变分，有

$$
dF=\sum_b dy_b\otimes g_b
$$

即

$$
dF_{kl}=\sum_b dy_{b,k}g_{b,l}
$$

单元中顶点$a$的内力为

$$
f_a=-V_0Pg_a
$$

因此

$$
df_a=-V_0dPg_a
$$

从而得到力对位置的雅可比块

$$
\frac{\partial f_{a,i}}{\partial y_{b,k}}
=-V_0\sum_{j,l}C_{ijkl}g_{a,j}g_{b,l}
$$

如果定义刚度矩阵为

$$
K=-\frac{\partial f}{\partial y}
$$

则单元刚度块为

$$
K_{ab,ik}=V_0\sum_{j,l}C_{ijkl}g_{a,j}g_{b,l}
$$

对于当前这个线性弹性模型，$C$不依赖当前$F$，因此$K$可以在初始形状下预计算。

记时间步长为$h$，它对应代码中每个物理子步的$\Delta t$，也就是

$$
h=\Delta t=\frac{dt}{\text{substeps}}
$$

隐式欧拉中，不考虑阻尼时，一步可以写成

$$
(M+h^2K)\Delta v=h f
$$

然后

$$
v_{n+1}=v_n+\Delta v\\
y_{n+1}=y_n+hv_{n+1}
$$

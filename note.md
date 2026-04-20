### 傅立叶变换

$$f(x)=\frac{a_0}{2}+\sum_{n=1}^\infty \left[a_n\cos(\frac{2\pi n}{T}x)+b_n\sin(\frac{2\pi n}{T}x)\right]$$

$$a_n=\frac{2}{T}\int_{-\frac{T}{2}}^{\frac{T}{2}}f(x)\cos(\frac{2\pi n}{T}x)$$

$$b_n=\frac{2}{T}\int_{-\frac{T}{2}}^{\frac{T}{2}}f(x)\sin(\frac{2\pi n}{T}x)$$

$\omega_n=\frac{2\pi n}{T}$

$e^{i\theta}=\cos\theta+i\sin\theta$

$e^{i\omega_nx}=\cos(\omega_n x)+i\sin(\omega_n x)$

$e^{-i\omega_nx}=\cos(\omega_n x)-i\sin(\omega_n x)$

$\cos(w_n x)=\frac{e^{i\omega_n x}+e^{-i\omega_nx}}{2}$

$\sin(w_n x)=-i\frac{e^{i\omega_n x}-e^{-i\omega_nx}}{2}$

$a_n\cos(\omega_n x)+b_n\sin(\omega_n x)=c_n(e^{i\omega_n x}+e^{-i\omega_n x})$

$$f(x)=\sum_{n=-\infty}^\infty c_n e^{i\omega_n x}$$

$$c_n=\frac{\int_{-\frac{T}{2}}^{\frac{T}{2}}f(x) e^{-i\omega_n x}}{\int_{-\frac{T}{2}}^{\frac{T}{2}}e^{i\omega_n x} e^{-i\omega_n x}}=\frac{1}{T}\int_{-\frac{T}{2}}^{\frac{T}{2}}f(x) e^{-i\omega_n x}$$

$$f(x)=\sum_{n=-\infty}^\infty \left(\frac{1}{T}\int_{-\frac{T}{2}}^{\frac{T}{2}}f(t) e^{-i\omega_n t}dt \right)e^{i\omega_n x}$$

$\omega_n$ 为角频率，$\Delta \omega=\frac{2\pi}{T}$，$T=\frac{2\pi}{\Delta \omega}$

所以有：

$$f(x)=\sum_{n=-\infty}^\infty \left(\frac{\Delta \omega}{2\pi}\int_{-\frac{T}{2}}^{\frac{T}{2}}f(t) e^{-i\omega_n t}dt \right)e^{i\omega_n x}$$

$T\rightarrow \infty,\Delta\omega\rightarrow 0$，$\omega_n$ 连续

此时有：

$$f(x)=\int_{-\infty}^\infty \left(\frac{1}{2\pi}\int_{-\frac{T}{2}}^{\frac{T}{2}}f(t) e^{-i\omega t}dt \right)e^{i\omega x}d\omega$$

即：

$$f(x)=\frac{1}{2\pi}\int_{-\infty}^\infty F(\omega)e^{i\omega x}d\omega$$

$$F(\omega)=\int_{-\infty}^\infty f(x)e^{-i\omega x}dx$$

如果是 $x$ 是离散的，则有：

$$f(x)=\sum_{n=-\infty}^\infty c_n e^{i\omega_n x}$$

$$X_k=\sum_{n=0}^{N-1}x_n e^{-i2\pi kn/N}$$

连续情形里写的是：$e^{-i\omega x}$

离散化后，$x$ 变成采样点编号 $n$，$\omega$ 变成离散频率 $\frac{2\pi k}{N}$

---

### FFT 求多项式乘法：$A\cdot B=C$

用系数表示法暴力求两个项数为 $n$ 的多项式相乘的复杂度为 $O(n^2)$，这很慢。

考虑用点值表示法来优化。

比如我们知道多项式在很多点上的值，比如：

$$A(x_0),A(x_1),A(x_2),...,A(x_n)\\ B(x_0),B(x_1),B(x_2),...,B(x_m)$$

那么我们可以得到

$$C(x_i)=A(x_i)B(x_i)$$

如果我们能根据点值表示法 $C(x_i)$ 反推出系数表示法，就完成了。

即：$系数表示法\rightarrow 点值表示法\rightarrow 系数表示法$

傅立叶变换本质上是时域与频域的转换，我们也可以将其应用到系数表示法与点值表示法的转换上。

如果我们选择的 $k$ 个点为单位根 $\omega_N^k=e^{\frac{2\pi k}{N}i}$，那么对应的函数值为：

$$A(\omega_N^k)=\sum_{n=0}^{N-1}a_n(\omega_N^k)^n=\sum_{n=0}^{N-1}a_ne^{\frac{2\pi kn}{N}i}\\ B(\omega_N^k)=\sum_{n=0}^{N-1}b_n(\omega_N^k)^n=\sum_{n=0}^{N-1}b_ne^{\frac{2\pi kn}{N}i}$$

对应点值相乘可得：

$$C(\omega_N^k)=\sum_{n=0}^{N-1}c_n(\omega_N^k)^n=\sum_{n=0}^{N-1}c_ne^{\frac{2\pi kn}{N}i}$$

其中 $c_n$ 就是我们要求的东西。

我们发现，这和傅立叶变换的形式高度一致，于是 $c_n$ 有：

$$c_n=\frac{\sum_{n=0}^{N-1}C(\omega_N^k) e^{-\frac{2\pi kn}{N}i}}{\sum_{n=0}^{N-1}e^{-\frac{2\pi kn}{N}i}e^{\frac{2\pi kn}{N}i}}=\frac{1}{N}\sum_{n=0}^{N-1}C(\omega_N^k) e^{-\frac{2\pi kn}{N}i}$$

所以我们只需要把 $C(\omega_N^k)$ 当作系数 $t_n$，令 $\omega_N^k=e^{-\frac{2\pi k}{N}i}$ 再求一遍即可。 

即：

$$c_n=\frac{1}{N}\sum_{n=0}^{N-1}t_n(\omega_N^k)^n$$

现在我们发现，$A(\omega_N^k),A(\omega_N^k),c_n$ 暴力求的复杂度为 $O(n^2)$，而 FFT 算法则利用单位根 $\omega_N^k$ 的一些性质将复杂度降低到 $O(n\log n)$，具体如下：

FFT 要求多项式项数为 $2$ 的整数幂，考虑用 $0$ 补齐。

假如多项式为：

$$f(x) = a_0 + a_1x + a_2x^2 + a_3x^3 + a_4x^4 + a_5x^5 + a_6x^6 + a_7x^7$$

按照次数的奇偶来分成两组，然后右边提出来一个 $x$：

$$
\begin{aligned}
f(x) &= (a_0 + a_2x^2 + a_4x^4 + a_6x^6) + (a_1x + a_3x^3 + a_5x^5 + a_7x^7) \\
     &= (a_0 + a_2x^2 + a_4x^4 + a_6x^6) + x(a_1 + a_3x^2 + a_5x^4 + a_7x^6)
\end{aligned}
$$

分别用奇偶次次项数建立新的函数：

$$
\begin{aligned}
G(x) &= a_0 + a_2x + a_4x^2 + a_6x^3 \\
H(x) &= a_1 + a_3x + a_5x^2 + a_7x^3
\end{aligned}
$$

那么原来的 $f(x)$ 用新函数表示为：

$$f(x) = G \left( x^2 \right) + x \times H \left( x^2 \right)$$

利用偶数次单位根的性质 $\omega_n^i = -\omega_n^{i+n/2}$，和 $G \left( x^2 \right)$ 和 $H \left( x^2 \right)$ 是偶函数，我们知道在复平面上 $\omega_n^i$ 和 $\omega_n^{i+n/2}$ 的 $G(x^2)$ 的 $H(x^2)$ 对应的值相同. 得到：

$$
\begin{aligned}
f(\omega_n^k) &= G((\omega_n^k)^2) + \omega_n^k \times H((\omega_n^k)^2) \\
              &= G(\omega_n^{2k}) + \omega_n^k \times H(\omega_n^{2k}) \\
              &= G(\omega_{n/2}^k) + \omega_n^k \times H(\omega_{n/2}^k)
\end{aligned}
$$

和：

$$
\begin{aligned}
f(\omega_n^{k+n/2}) &= G(\omega_n^{2k+n}) + \omega_n^{k+n/2} \times H(\omega_n^{2k+n}) \\
                    &= G(\omega_n^{2k}) - \omega_n^k \times H(\omega_n^{2k}) \\
                    &= G(\omega_{n/2}^k) - \omega_n^k \times H(\omega_{n/2}^k)
\end{aligned}
$$

因此我们求出了 $G(\omega_{n/2}^k)$ 和 $H(\omega_{n/2}^k)$ 后，就可以同时求出 $f(\omega_n^k)$ 和 $f(\omega_n^{k+n/2})$. 于是对 $G$ 和 $H$ 分别递归 DFT 即可.

FFT.cpp:

```cpp
#include<bits/stdc++.h>
using namespace std;
typedef double db;
const int N = 5e6 + 5;
const db pi = acos(-1);
struct node{
    db x, y;
    node(db x = 0.0, db y = 0.0):x(x), y(y){}
};
node operator + (node a, node b) {return node(a.x + b.x, a.y + b.y);}
node operator - (node a, node b) {return node(a.x - b.x, a.y - b.y);}
node operator * (node a, node b) {return node(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);}
void FFT(node *a, int n, int tp) {
    if (n == 1) return;
    node even[n >> 1], odd[n >> 1];
    for (int i = 0; i < n; i += 2) {
        even[i >> 1] = a[i];
        odd[i >> 1] = a[i + 1];
    }
    FFT(even, n >> 1, tp);
    FFT(odd, n >> 1, tp);
    node Wn = node(cos(2 * pi / n), tp * sin(2 * pi / n)), w = node(1, 0);
    for (int i = 0; i < (n >> 1); ++i, w = w * Wn) {
        node E = even[i], O = w * odd[i];
        a[i] = E + O;
        a[i + (n >> 1)] = E - O;
    }
}
int n, m, len = 1, cnt;
node a[N], b[N], c[N];
int main() {
    ios::sync_with_stdio(false);
    cin >> n >> m;
    for (int i = 0; i <= n; ++i) cin >> a[i].x;
    for (int i = 0; i <= m; ++i) cin >> b[i].x;
    while(len <= n + m) len <<= 1, ++cnt;
    FFT(a, len, 1);
    FFT(b, len, 1);
    for (int i = 0; i < len; ++i) c[i] = a[i] * b[i];
    FFT(c, len, -1);
    for (int i = 0; i <= n + m; ++i) cout << (long long)(c[i].x / len + 0.5) << " ";
    return 0;
}
```

这种实现方式比较浪费空间，因为每一层都需要创建 `even` 和 `odd` 去递归。我们考虑将序列一开始就变成最终形态，然后从小往大合并，这个过程称为蝴蝶变换。

具体的，以 8 项多项式为例，模拟拆分的过程：

* 初始序列为 $\{x_0, x_1, x_2, x_3, x_4, x_5, x_6, x_7\}$
* 一次二分之后 $\{x_0, x_2, x_4, x_6\}, \{x_1, x_3, x_5, x_7\}$
* 两次二分之后 $\{x_0, x_4\}\{x_2, x_6\}, \{x_1, x_5\}, \{x_3, x_7\}$
* 三次二分之后 $\{x_0\}\{x_4\}\{x_2\}\{x_6\}\{x_1\}\{x_5\}\{x_3\}\{x_7\}$

规律：其实就是原来的那个序列，每个数用二进制表示，然后把二进制翻转对称一下，就是最终那个位置的下标. 比如 $x_1$ 是 001，翻转是 100，也就是 4，而且最后那个位置确实是 4. 

根据它的定义，我们可以在 $O(n)$ 的时间内求出每个数变换后的结果，

FFT_butterfly.cpp:

```cpp
#include<bits/stdc++.h>
using namespace std;
typedef double db;
const int N = 5e6 + 5;
const db pi = acos(-1);
struct node{
    db x, y;
    node(db x = 0.0, db y = 0.0):x(x), y(y){}
};
node operator + (node a, node b) {return node(a.x + b.x, a.y + b.y);}
node operator - (node a, node b) {return node(a.x - b.x, a.y - b.y);}
node operator * (node a, node b) {return node(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);}
int to[N];
void FFT(node *a, int n, int tp) {
    for (int i = 0; i < n; ++i) 
        if (i < to[i])
            swap(a[i], a[to[i]]);
    for (int t = 1; t < n; t <<= 1) {
        node Wn = node(cos(pi / t), tp * sin(pi / t)), w, E, O;
        for (int len = t << 1, j = 0; j < n; j += len) {
            w = node(1, 0);
            for (int k = 0; k < t; ++k, w = w * Wn) {
                E = a[j + k];
                O = w * a[t + j + k];
                a[j + k] = E + O;
                a[t + j + k] = E - O;
            }
        }
    }
}
int n, m, len = 1, cnt;
node a[N], b[N], c[N];
int main() {
    ios::sync_with_stdio(false);
    cin >> n >> m;
    for (int i = 0; i <= n; ++i) cin >> a[i].x;
    for (int i = 0; i <= m; ++i) cin >> b[i].x;
    while(len <= n + m) len <<= 1, ++cnt;
    for (int i = 0; i < len; ++i) to[i] = (to[i >> 1] >> 1) | ((i & 1) << (cnt - 1));
    FFT(a, len, 1);
    FFT(b, len, 1);
    for (int i = 0; i < len; ++i) c[i] = a[i] * b[i];
    FFT(c, len, -1);
    for (int i = 0; i <= n + m; ++i) cout << (long long)(c[i].x / len + 0.5) << " ";
    return 0;
}
```
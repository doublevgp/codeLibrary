## 快速傅里叶变换FFT(Fast Fourier Transform)多用于多项式乘法以及大整数乘法

### 多项式乘法

给定A(n + 1), B(m + 1)数组，一般是系数表达式，也可能会是点值表达式。

求C(n + m + 1) = A(x) * B(x)

常规的做法是O(nm)级别的暴力枚举

复杂度过高，所以我们需要降低其复杂度

FFT的思路其实就是要确定一个n + m 次的多项式，则需要使用n + m + 1个不同的点来确定，所以我们可以从A中选出n + m + 1个点，再从B中选出n + m + 1个点，再将这些点对应相乘，生成n+m+1个点，这就是C的点值表达，再要做的就是将点值表达式转回系数表达式。

但是这样选点求点值的话复杂度也是比较高的，所以我们需要再次降，我们考虑，$y=x^2$的时候，当$x=1$时，$y(1) = y(-1) = 1$，$y=x^3$时，$y(-1)=-y(1)=-1$，我们发现，我们只需要选择一个点值，那么其相反数的点值我们也能得到.

我们采用单位根的形式来进行分治 $W^{2j}_{n} = W^{j}_{n / 2}$

如果常数过大

可以采用蝴蝶操作来将分治后的顺序保存下来

```cpp
m += n;
n = 1;
int cnt = 0;
while(n <= m) {
  n <<= 1;
  cnt += 1;
}

for (int i = 0; i < n; i++) {
  r[i] = (r[i >> 1] >> 1) | ((i & 1) << (cnt - 1));
}
```



FFT操作

```cpp
void FFT(cp *t, const int &n, int inv) {
  for (int i = 0; i < n; i++) {
    if (i < r[i]) swap(t[i], t[r[i]]);
  }
  
  for (int len = 1; len <= (n >> 1); len <<= 1) {
    cp w1(cos(PI / len), inv * sin(PI / len));
    for (int i = 0; i <= n - (len << 1); i += (len << 1)) {
      cp w(1.0, 0.0);
      for (int j = 0; j <= len - 1; j++) {
        cp x = t[i + j], y = w * t[i + j + len];
        t[i + j] = x + y;
        t[i + j + len] = x - y;
        w = w * w1;
      }
    }
  }
}
```



Main

```cpp
FFT(a, n, 1);
FFT(b, n, 1);
for (int i = 0; i <= n; i++) {
  c[i] = a[i] * b[i];
}
FFT(c, n, -1);
for (int i = 0; i <= m; i++) {
  cout << (int)(c[i].r / x + 0.5) << " ";
}
cout << endl;
```


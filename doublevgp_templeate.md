## 倍增ST表

## 乘法逆元

$A·X\equiv1\pmod{p}$

其中 当A 与 p互质时，X有解，即存在A的逆元X，称X为A的逆元

求逆元 可采用扩展欧几里德算法或费马小定理

### 1、扩展欧几里德算法

```cpp
// cpp
inline ll extend_gcd(ll a, ll b, ll x, ll y) {
  if (a == 0 && b == 0) {
    return -1ll;
  }
  if (b == 0) {
    x = 1ll;
    y = 0ll;
    return a;
  }
  ll d = extend_gcd(b, a % b, y, x);
  y -= a / b * x;
  return d;
}
inline ll mod_reverse(ll a, ll n) {
  ll x, y, d = extend_gcd(a, n, x, y);
  if (d == 1) {
    if (x % y <= 0) {
      return x % n + n;
    }
    return x % n;
  }
  return -1ll;
}
```

### 2、费马小定理 

$a^{(p-1)}\equiv1(mod\space p)$
$$
前提a\%p=1,p为质数,令arr[i]=i(0<i<p),则a*arr[i]\%p=1,2,...,(p-2),(p-1) \\
例如a=2,p=5,arr=[1,2,3,4],a*arr[i]\%p=[2,4,1,3]
$$

$$
\begin{align*}
&\because 1\%p*2\%p*...*(p-1)\%p=a*1\%p*a*2\%p*...*a*(p-1)\%p \\
&\therefore(p-1)!\%p=a^{(p-1)}*(p-1)!\%p \\
&a^{(p-1)}\equiv1(mod\;p) \\
&a^{(p-2)}*a\equiv1(mod\;p) \qquad两边同除a \\
&a^{(p-2)}\equiv a^{-1}(mod\;p)
\end{align*}
$$

```cpp
// cpp
ll ksm(ll a, ll k, ll mod) {
  ll ans = 1;
  while (k) {
    if (k & 1) {
      ans = ans * a % mod;
    }
    a = a * a % mod;
    k >>= 1ll;
  }
  return ans;
}
ll mod_reverse(ll a, ll mod) {
  return ksm(a, mod - 2, mod);
}
```

### 欧拉定理求逆

观察式子 $ax\equiv1(mod\;p)$ ，我们发现它和欧拉定理 $a^{\varphi(p)}\equiv1(mod\;p)$ 很像……

而欧拉定理又是 $a,p$ 互质时成立，这和逆元存在的条件相同……

故将两个式子合在一起，得：

$ax\equiv a^{\varphi(p)}\quad (mod\;p)$

两遍同除$a$ 得：

$x\equiv a^{\varphi(p)-1} \quad(mod\;p)$

我们可以在 $O(\sqrt{p})$ 的时间内求出单个数的欧拉函数。

> 引理：设 $n=p^{k_1}_1\times p^{k_2}_2\times \cdots \times p^{k_m}_m=\prod^{m}_{i=1}p^{k_i}_i$为正整数 $n$ 的质数幂乘积表示式，则：
> $\varphi(n)=n\times (1−\frac{1}{p_1})\times (1−\frac{1}{p_2})\times \cdots \times (1−\frac{1}{p_m})=n \times \prod^{m}_{i=1}(1−\frac{1}{p_i})$

证明：

由于各质数幂之间是互质的，根据欧拉函数的性质可得：

$\varphi(n)=\prod\limits_{i=1}^m \varphi(p^{k_i}_i)=\prod\limits_{i=1}^m p^{k_i}_i \times \prod\limits_{i=1}^m (1-\frac{1}{p_i})=n\times \prod\limits_{i=1}^m(1-\frac{1}{p_i})$

进而得证。

在求一个数的欧拉函数时，不妨将上式后面的分数部分通分化作 $\varphi(n)=n\times \prod^m_{i=1}(\frac{p_i-1}{p_i})$ ，我们可以在$ O(\sqrt{n})$ 的时间内枚举 $n$ 的每个质因子，初始化 $ans$为 $n$ ，每枚举到一个质因子 $p_i$ 就让 $ans\times \frac{p_i-1}{p_i}$ ，发现每个质因子只对答案贡献一次，故在统计完该质因子的答案后，这个质因子就无用了，为了方便最后判断我们是否将 $n$ 的所有质因子全都枚举到，我们可以通过不断除该质因子来消去它，防止其再被计入答案。

求出欧拉函数后，可以使用快速幂求出答案。

故这样做的时间复杂度为 $O(log\,\varphi(p)+\sqrt{p})$ 。

模板题参考代码：

### 分数取模

要求$\frac{1}{a}\;mod \;p$，就是求$a^{p-2}mod\;p$，如果求$\frac{m}{a}\;mod\;p$，则可以使用$(m*a^{-1})\;mod\;p$再将乘积的模转换成模的乘积$((m\%p)*(a^{-1}\%p))\;mod\;p$

```cpp
// cpp
ll frac_mod(ll a, ll b, ll mod) { // a/b % p
    ll ans = 1ll * a * ksm(b, mod - 2, mod) % mod;
    return ans;
}
```

### 中国剩余定理



### 扩展Lucas

扩展 $Lucas$ 一般是用来解决 $C^m_n\;mod\;p\;(p\;\in composite)$ 的组合数取模问题。

我们对 $p$ 进行质因数分解得到 $p$ 的质因数分解形式 $p=\prod^k_{i=1}p^{m_i}_i$ 。如果我们对每个质因子 $p_i$ 解决了$C^m_n\;mod\;p^{m_i}_i$ ，再通过中国剩余定理我们就能解出 $C^m_n\;mod\;p$。

对于 $C^m_n\;mod\;p^{m_i}_i$  ，因为 $C^m_n=\frac{n!}{m!(n-m)!}$ ，所以下面的阶乘部分由于可能含有$p_i$的幂次因此不能直接求逆元得到，我们考虑先提出质因子 $p_i$ 。

我们可以得到 $\frac{n!}{m!(n-m)!}=\frac{\frac{n!}{p^{k_1}_i}}{\frac{m!}{p^{k_2}_i}\frac{(n-m)!}{p^{k_3}_i}}p^{k_1-k_2-k_3}_i$。在提出$p_i$的幂次之后我们就可以对分母求逆元计算出在模$p^{m_i}_i$意义下的答案。现在只需要解决 形如$\frac{n!}{p^{k_i}_i}\;mod\;p^{m_i}_i$这样的问题就行了。

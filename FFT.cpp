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
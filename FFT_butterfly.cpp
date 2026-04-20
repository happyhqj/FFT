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
#include <bits/stdc++.h>

using namespace std;

typedef long double ld;
typedef complex<ld> c;
// polynomial
typedef vector<c> poly;

const ld pi = 3.1415926535897932384626433832795;

// fft of polynomial f, processing n entries starting from a with common difference d
// inv: whether it's inverse fft
void fft(poly &f, int n, int a, int d, bool inv) {
    if(n==1) return;

    // even and odd
    fft(f, n/2, a, d*2, inv);
    fft(f, n/2, a+d, d*2, inv);

    // combine
    poly res(f.size(), 0);
    // w = n-th root of unity
    ld arg = 2*pi*(inv?1:-1) / n;
    c wn = c(cos(arg), sin(arg));
    c w = c(1, 0);
    for(int i=0; i<n/2; i++) {
        c even = f[a+(2*i)*d];
        c odd = w * f[a+(2*i+1)*d];
        res[a+i*d] = even + odd;
        res[a+(i+n/2)*d] = even - odd;
        // in inverse fft, we do the a=even+odd, b=even-odd step in reverse
        // which yields even=(a+b)/2 and odd=(a-b)/2
        if(inv) {
            res[a+i*d] /= 2;
            res[a+(i+n/2)*d] /= 2;
        }
        w *= wn;
    }
    for(int i=0; i<n; i++){
        f[a+d*i] = res[a+d*i];
    }
}

poly mult(poly &f, poly &g) {
    // length must be 2^n
    int len = 1 << (__lg(f.size() + g.size()) + 1);
    f.resize(len);
    g.resize(len);
    fft(f, len, 0, 1, false);
    fft(g, len, 0, 1, false);
    poly res;
    for(int i=0; i<len; i++) {
        res.push_back(f[i]*g[i]);
    }
    fft(res, len, 0, 1, true);
    return res;
}

// input: two lines
//     first line: an integer m, then m numbers f_i
//     denoting f = \sum f_i x^i
//     second line: an integer n, then n numbers g_i
//     denoting g = \sum g_i x^i
// output: one line containing m+n-1 numbers p_i, representing the product p = \sum p_i x^i
// time complexity: O(nlogn), where n is the greater of f and g's degrees
int main() {
    int m, n;
    ld x;
    poly f, g;
    cin >> m;
    for(int i=0; i<m; i++) {
        cin >> x;
        f.push_back(x);
    }
    cin >> n;
    for(int i=0; i<n; i++) {
        cin >> x;
        g.push_back(x);
    }
    poly res = mult(f, g);
    cout << fixed << setprecision(4);
    for(int i=0; i<m+n-1; i++) {
        if(i) cout << ' ';
        cout << res[i].real();
    }
    cout << '\n';
    return 0;
}

typedef long double ld;
typedef complex<ld> c;
typedef vector<c> poly;
const ld pi = 3.141592653589793238;

int n;

void fft(poly &f, int a, int d, bool inv) {
    if(d==n) return;
    fft(f, a, d*2, inv);
    fft(f, a+d, d*2, inv);
    poly res(f.size());
    ld r = 2*pi*(inv?d:-d)/n;
    c w = c(1, 0);
    for(int i=a; i<n/2; i+=d) {
        c even = f[2*i-a];
        c odd = w * f[2*i+d-a];
        res[i] = even + odd;
        res[i+n/2] = even - odd;
        w *= c(cos(r), sin(r));
    }
    for(int i=a; i<n; i+=d)
        f[i] = res[i] / (c)(inv?2:1);
}

poly mult(poly &f, poly &g) {
    n = 1<<(__lg(f.size()+g.size())+1);
    f.resize(n);
    g.resize(n);
    fft(f, 0, 1, 0);
    fft(g, 0, 1, 0);
    poly res;
    for(int i=0; i<n; i++)
        res.push_back(f[i]*g[i]);
    fft(res, 0, 1, 1);
    return res;
}

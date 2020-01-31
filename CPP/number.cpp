// ТЕОРИЯ ЧИСЕЛ
typedef long long ll;
// делители числа
void nd(ll n){
    for (ll i = 2; i <= sqrt(n); i++) {
        if (n % i == 0) {
            // i, n / i - делители числа n
        }
    }
}
// простые делители числа
void pnd(ll n){
    for (ll i = 2; i <= n; i++) {
        if (n % i == 0) {
            while (n % i == 0) {
                n /= i;
            }
            // i - простой делитель
        }
    }
}

//НОД алгоритм Евклида
ll gcd(ll a, ll b) {
    while (a > 0 && b > 0) {
        if (a > b) {
            a %= b;
        }
        else {
            b %= a;
        }
    }
    return max(a, b);
}

// НОК\
НОК(a, b) = a*b / НОД(a, b)
ll lcm(ll a, ll b){
    return a * b / gcd(a, b);
}

// Решето Эратосфена
vector<bool> Eratos_Sieve(ll n) {
    vector<bool> vis(n + 1, true);
    vis[0] = vis[1] = false;
    for (ll i = 2; i <= n; i++) {
        if (vis[i]) {
            for (int j = 2 * i; j <= n; j += i) {
                vis[j] = false;
            }
        }
    }
    return vis;
}

// if(abs(a - b) < EPS) => то числа равны с определенной точностью

// Bin pow
ll bp(ll a, ll n, ll mod) {
    if (n == 1) {
        return a;
    }
    ll z = bp(a, n / 2, mod);
    z = (z * z) % mod;
    if (n % 2 == 0) {
        return z;
    }
    else {
        return (z * a) % mod;
    }
}


// i-ое число Фибоначи за log
struct matrix {
    ll a[2][2];
};
matrix multi(matrix a, matrix b) {
    matrix c;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            c.a[i][j] = a.a[i][0] * b.a[0][j] + a.a[i][1] * b.a[1][j];
        }
    }
    return c;
}
matrix binPow(matrix a, ll n) {
    if (n == 1) {
        return a;
    }
    matrix z = binPow(a, n / 2);
    z = multi(z, z);
    if (n % 2 == 0) {
        return z;
    }
    else {
        return multi(z, a);
    }
}
ll Fib(ll n) {
    matrix a = { 0, 1, 1, 1 };
    a = binPow(a, n);
    return a.a[0][0];
}

// Для вывода double(чисел с плавающей точкой) с определенной точностью используют\
cout << fixed << setprecision(точность) << value;

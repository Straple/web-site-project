// Вычислительная Геометрия
#define x first
#define y second

typedef long double ld;
typedef long long ll;

const ld eps = 1e-9;

struct dot { // структура точка
    ll x, y;
    dot() {// обнуление
        x = y = 0;
    }
    dot(ll _x, ll _y) {
        x = _x;
        y = _y;
    }
    dot operator + (dot p) {
        return dot(x + p.x, y + p.y);
    }
    dot operator - (dot p) {
        return dot(x - p.x, y - p.y);
    }
    ll operator % (dot p) { // векторное произведение
        return x * p.y - y * p.x;
    }// от знака векторного произведения зависит расположение точки p

    ll operator * (dot p) { // скалярное произведение
        return x * p.x + y * p.y;
    }

    dot operator * (ll k) { // умножение вектора на число k
        return dot{ x * k, y * k };
    }

    ll getlen() { // длина вектора
        return sqrt(x * x + y * y);
    }

    ll getsquarelen() { // длина вектора^2
        return x * x + y * y;
    }

    bool operator == (dot p) {
        return (p.x == x) && (p.y == y);
    }

    pair<ld, ld> normalize() { // предает x и y значения 0 или 1
        return make_pair(x / getlen(), y / getlen());
    }

    bool iscent(dot p1, dot p2) { // точка находиться между p1 и p2
        if ((p2.x <= x && x <= p1.x) || (p1.x <= x && x <= p2.x)) {
            if ((p1.y <= y && y <= p2.y) || (p2.y <= y && y <= p1.y)) {
                return true;
            }
        }
        return false;
    }
};

bool iscent(dot p1, dot p2, pair<ld, ld> c) // тоже что и чуть выше, только с eps
{
    if ((c.x <= p1.x + eps && c.x + eps >= p2.x) || (c.x + eps >= p1.x && c.x <= p2.x + eps)) {
        if ((c.y <= p1.y + eps && c.y + eps >= p2.y) || (c.y + eps >= p1.y && c.y <= p2.y + eps)) {
            return 1;
        }
    }
    return 0;
}

dot infv = dot(LLONG_MAX, LLONG_MAX);

const ld pi = acos(-1);

ld getangle(dot a, dot b) { // возвращает угол между векторами
    return atan2(a % b, a * b);
}

ld getgoodangle(dot a, dot b) { // возвращает не отрицательный угол
    ld res = atan2(a % b, a * b);
    if (res < 0) {
        res += 2 * pi;
    }
    return res;
}

double getverygoodangle(dot a, dot b) { // возвращает не отрицательный угол меньше 180
    double res = atan2(a % b, a * b);
    if (res < 0) {
        res += 2 * pi;
    }
    if (res > pi) {
        res = 2 * pi - res;
    }
    return res;
}

struct line { // структура линия; задается a, b, c, где A*x + B*y + C = 0 \
            и двумя точками лежащими на прямой
    ll a, b, c;
    dot on1, on2;
    line(dot p, dot q) { // по двум точкам инициализируем значения line
        a = p.y - q.y;
        b = q.x - p.x;
        c = -a * p.x - b * p.y; // C = -A*x - B*y
        on1 = p;
        on2 = q;
    }
    line(ll _a, ll _c) { // инициализация с помощью A и C
        a = _a;
        b = -1;
        c = _c;
        on1 = dot(1, a + c);
        on2 = dot(2, 2 * a + c);
    }
    line(ll _a, ll _b, ll _c) { // инициализация A,B,C 
        a = _a;
        b = _b;
        c = _c;
    }
    bool ison(dot p) { // если точка p находится примерно на линии
        if (abs(a * p.x + b * p.y + c) < eps) {
            return true;
        }
        return false;
    }
    pair<ld, ld> intersect(line l) { // возвращают точку пересечения линий
        ld x, y;
        if (l.b != 0) {
            x = ((((ld)b) * l.c / l.b - c) / (a - ((ld)b) * l.a / l.b));
            y = (-x * l.a - l.c) / l.b;
        }
        else {
            x = -((ld)l.c) / l.a;
            y = (-x * a - c) / b;
        }
        return make_pair(x, y);
    }
};

ld getdist(line l, dot a) { // расстояние между линией и точкой
    if (l.ison(a)) {
        return 1;
    }
    dot parr = l.on1 - l.on2; // берет вектор от l
    dot norm = dot(-parr.y, parr.x);// берет перпендикулярный вектор 
    dot r = l.on1 - a;
    double ang = getverygoodangle(r, norm);
    if (ang > pi / 2) {
        norm = norm * -1;
        ang = pi - ang;
    }
    ang = pi / 2 - ang;
    return r.getlen() * sin(ang);
}

ld getdist2(line l, dot p) {
    ld A = l.on1.y - l.on2.y,
        B = l.on1.x - l.on2.x,
        C = l.on1.x * l.on2.y - l.on1.y * l.on2.x;
    ld d = A * p.x - B * p.y + C;

    return abs(d / sqrt(A * A + B * B));
}

dot perp(line l, dot a) { // почти тоже что и чуть выше, но возвращает точку пересечения
    dot parr = l.on1 - l.on2;
    dot norm = dot(-parr.y, parr.x);
    dot r = l.on1 - a;
    double ang = getverygoodangle(r, norm);
    if (ang > pi / 2) {
        norm = norm * -1;
        ang = pi - ang;
    }
    ang = pi / 2 - ang;
    double dd = r.getlen() * sin(ang);
    dot ddd = norm * dd;
    return a + ddd;
}

bool isParallel(line a, line b) { // if A1*B2 - A2*B1 == 0 - line is parallel
    return a.a * b.b - a.b * b.a == 0;
}
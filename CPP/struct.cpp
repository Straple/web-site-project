// ДД (Декартово дерево, Декартач, Курево, trep(tree + heap), Дерамида - бинарное дерево по x и куча по y)
    // Если y - random number, то ДД превращается в сбалансированное или рандомизированное бинарное дерево поиска
#include <random>
//#define x first\
#define y second

struct node;
typedef node* pnode;
mt19937_64 rnd(42);
pnode root = nullptr;

struct node { // ДД
    int cnt;
    int sz;
    int val;
    long long prior;
    pnode l, r; // ветки
    node() {}
    node(int k)
    {
        cnt = 1;
        val = k;
        prior = rnd();
        sz = 1;
        l = r = nullptr;
    }
};
int getsz(pnode t) // возвращает кол-во элементов в t
{
    if (t == nullptr) {
        return 0;
    }
    return t->sz;
}
void upd(pnode t) // обновление размера в t
{
    if (t != nullptr) {
        t->sz = getsz(t->l) + getsz(t->r) + 1;
    }
}

//Merge(T1, T2) - объединяет два поддерева T1 и T2, и возвращает новое дерево.\
    Эта операция также реализуется за O(log N).Она работает в предположении, \
    что T1 и T2 обладают соответствующим порядком(все значения X в первом меньше значений X во втором).\
    Таким образом, нам нужно объединить их так, чтобы не нарушить порядок по приоритетам Y.\
    Для этого просто выбираем в качестве корня то дерево, у которого Y в корне больше,\
    и рекурсивно вызываем себя от другого дерева и соответствующего сына выбранного дерева.
pnode merge(pnode t, pnode v)
{
    if (t == nullptr) {
        return v;
    }
    else if (v == nullptr) {
        return t;
    }
    else if (t->prior > v->prior) {
        t->r = merge(t->r, v);
        upd(t);
        return t;
    }
    else {
        v->l = merge(t, v->l);
        upd(v);
        return v;
    }
}

//Split(T, X) - разделяет дерево T на два дерева L и R(которые являются возвращаемым значением)\
    таким образом, что L содержит все элементы, меньшие по ключу X, а R содержит все элементы, большие X.Эта операция выполняется за O(log N).
pair<pnode, pnode> split(pnode t, int x)
{
    if (t == nullptr) {
        return make_pair(nullptr, nullptr);
    }
    else if (getsz(t->l) >= x) {
        pair<pnode, pnode> v = split(t->l, x);
        t->l = v.y;
        upd(t);
        return make_pair(v.x, t);
    }
    else {
        pair<pnode, pnode> v = split(t->r, x - 1 - getsz(t->l));
        t->r = v.x;
        upd(t);
        return make_pair(t, v.y);
    }
}

//Теперь очевидна реализация Insert (X, Y). Сначала спускаемся по дереву (как в обычном бинарном дереве поиска по X),\
    но останавливаемся на первом элементе, в котором значение приоритета оказалось меньше Y.\
    Мы нашли позицию, куда будем вставлять наш элемент. Теперь вызываем Split (X) от найденного элемента\
    (от элемента вместе со всем его поддеревом), и возвращаемые ею L и R записываем в качестве левого и правого сына добавляемого элемента.\
    \
    Также понятна и реализация Erase (X). Спускаемся по дереву (как в обычном бинарном дереве поиска по X), ища удаляемый элемент.\
    Найдя элемент, мы просто вызываем Merge от его левого и правого сыновей, и возвращаемое ею значение ставим на место удаляемого элемента.\

pair<pnode, pnode> split(pnode t, int val) // сплит по ключу
{
    if (t == nullptr) {
        return make_pair(nullptr, nullptr);
    }
    pnode l;
    pnode r;
    upd(t);
    if (t->val > val) {
        pair<pnode, pnode> pr = split(t->l, val);
        t->l = pr.y;
        upd(t);
        return make_pair(pr.x, t);
    }
    else {
        pair<pnode, pnode> pr = split(t->r, val);
        t->r = pr.x;
        upd(t);
        return make_pair(t, pr.y);
    }
}

int getkth(pnode t, int k)
{
    if (getsz(t->l) < k) {
        return getkth(t->r, k - getsz(t->l) - 1);
    }
    else if (getsz(t->l) == k) {
        return t->val;
    }
    else
        return getkth(t->l, k);
}

pnode findkey(pnode t, int x)
{
    if (t == nullptr) {
        return nullptr;
    }
    else if (t->val == x) {
        return t;
    }
    else if (t->val > x) {
        return findkey(t->l, x);
    }
    else return findkey(t->r, x);
}

pnode ins(pnode t, int x)
{
    if (t == nullptr)
        return new node(x);
    if (findkey(t, x) == nullptr) {
        pair<pnode, pnode> v = split(t, x);
        t = merge(v.x, new node(x));
        t = merge(t, v.y);
        return t;
    }
    (findkey(t, x)->cnt)++;
    return t;
}

pnode del(pnode t, int x)
{
    if (findkey(t, x) == nullptr) {
        return t;
    }
    if (findkey(t, x)->cnt == 1) {
        pair<pnode, pnode> v = split(t, x);
        pair<pnode, pnode> u = split(v.x, x - 1);
        return merge(u.x, v.y);
    }
    (findkey(t, x)->cnt)--;
    return t;
}


// Хеш строки
typedef long long ll;
const int mod = 1000, basis = 1000;
struct Hash {
    ll n, _n;
    vector<int> _hash, _powBasis;
    ll _mod(ll x) {
        return ((x % mod) + mod) % mod;
    }
    ll _pow(ll x, ll k) {
        if (k == 0) {
            return 1;
        }
        else if (k == 1) {
            return x;
        }
        else {
            ll t = _pow(x, k / 2);
            return _mod(_mod(t * t) * (k % 2 == 0 ? 1 : x));
        }
    }
    void count_n() {
        _n = _pow(basis, (ll)(mod - 2) * (ll)n);
    }
public:
    Hash(string& s) {
        n = s.length();
        _hash.resize(n + 1);
        _powBasis.resize(n + 1);
        _powBasis[0] = 1;
        for (int i = 0; i < n; ++i) {
            _powBasis[i + 1] = _mod((ll)_powBasis[i] * (ll)basis);
        }
        _hash[0] = 0;
        for (int i = 0; i < n; ++i) {
            _hash[i + 1] = _mod(_hash[i] + _mod((ll)s[i] * (ll)_powBasis[i]));
        }
        count_n();
    }

    int getHashSubStr(int l, int r) { // [l, r)
        return _mod((ll)_mod((ll)_mod(_hash[r] - _hash[l]) * (ll)_powBasis[n - l]) * (ll)_n);
    }
};

// Система непересекающихся множеств(СНМ)

struct dsu {
    int* p;
    int* r;

    //добавляет новый элемент x, помещая его в новое множество, состоящее из одного него.
    dsu(int n) {
        p = new int[n];
        r = new int[n];
        for (int i = 0; i < n; i++) {
            p[i] = i;
            r[i] = 1;
        }
    }
    //возвращает, в каком множестве находится указанный элемент u. На самом деле при этом возвращается один из элементов множества \
        (называемый представителем или лидером. Этот представитель выбирается в каждом множестве самой структурой данны\
        х (и может меняться с течением времени, а именно, после вызовов {\rm union\_sets}()).
    int get(int u) {
        if (p[u] == u) return u;
        return p[u] = get(p[u]);
    }

    //объединяет два указанных множества(множество, в котором находится элемент u, и множество, в котором находится элемент v).
    bool join(int u, int v) {
        u = get(u);
        v = get(v);
        if (u != v) {
            if (r[u] > r[v]) {
                p[v] = u;
                r[u] = max(r[u], r[v] + 1);
            }
            else {
                p[u] = v;
                r[v] = max(r[v], r[u] + 1);
            }
            return true;
        }
        return false;
    }
};



//Фенвик или Дерево Фенвика
// отсчет массива с 1 !!!
typedef long long ll;
const ll maxn = 1e5;
ll t[maxn];

// возвращает сумму на префиксе [1, r]
ll sum(ll r) {
    ll res = 0;
    for (; r > 0; r -= r & -r) {
        res += t[r];
    }
    return res;
}

// sum [l, r]
ll getsum(ll l, ll r) {
    return sum(r) - sum(l - 1);
}

// обновляет нужные t\
    k - index, dt - разница между новым значением и старым, \
    или просто val, если добавляем новый
void add(ll n, ll k, ll dt) {
    for (; k <= n; k += k & -k) {
        t[k] += dt;
    }
}

// возвращает индекс, на котором сумма уже больше
ll lower_bound(ll n, ll s) {
    ll k = 0;
    for (ll l = log(n); l >= 0; l--) {
        if (k + (1 << l) <= n && t[k + (1 << l)] < s) {
            k += (1 << l);
            s -= t[k];
        }
    }
    return k;
}



//Спарс таблица (СПАРС)
//находит RMQ на отрезке за O(1), предпосчет и память за O(n*log(n))
const int maxn = 1e4;
const int logmaxn = 1e2; // log(maxn)
int a[maxn], lg[maxn], mn[maxn][logmaxn];

int rmq(int l, int r) { // полуинтервал [l; r) 
    int t = lg[r - l];
    return min(mn[l][t], mn[r - (1 << t)][t]);
}

int main()
{
    ifstream cin("input.txt");
    int n; cin >> n;
    for (int i = 0; i < n; i++)
    {
        cin >> a[i];
    }
    // строим log для таблицы
    for (int l = 0; l < log(n); l++) {
        for (int i = (1 << l); i < maxn; i++) {
            lg[i] = l;
        }
    }

    for (int i = n - 1; i >= 0; i--) {

        mn[i][0] = a[i];
        for (int l = 0; l < log(n) - 1; l++) {
            mn[i][l + 1] = min(mn[i][l], mn[i + (1 << l)][l]);
        }
    }
}


//Дерево отрезков(ДО)\
    позволяет выполнять ассоциативные операции(+, %, min, max) на отрезке за log, требует 4*n памяти\
    n - размер массива

typedef long long ll;
const int maxn = 1e4;
ll V[maxn], T[4 * maxn];

//build(1, 0, n - 1);  build Tree O(n*log(n))
void build(int v, int tl, int tr) {
    if (tl == tr) {
        T[v] = V[tl];
    }
    else {
        int tm = (tl + tr) / 2;
        build(v * 2, tl, tm);
        build(v * 2 + 1, tm + 1, tr);
        T[v] = T[v * 2] + T[v * 2 + 1];
    }
}

//getsum(1, 0, n - 1, l - 1, r - 1) -> получить сумму на отрезке [l, r]
ll getsum(int v, int tl, int tr, int l, int r) { // log n операций
    if (l > r) {
        return 0;
    }
    if (l == tl && r == tr) {
        return T[v];
    }
    int tm = (tl + tr) / 2;
    return getsum(v * 2, tl, tm, l, min(r, tm)) + getsum(v * 2 + 1, tm + 1, tr, max(l, tm + 1), r);
}

// обновляет значение элемента v на новое\
    update(1, 0, n - 1, pos, new_val);
void update(int v, int tl, int tr, int pos, int new_val) { // log n операций
    if (tl == tr) {
        T[v] = new_val;
    }
    else {
        int tm = (tl + tr) / 2;
        if (pos <= tm)
        {
            update(v * 2, tl, tm, pos, new_val);
        }
        else
        {
            update(v * 2 + 1, tm + 1, tr, pos, new_val);
        }
        T[v] = T[v * 2] + T[v * 2 + 1];
    }
}


// Персистентное ДО
typedef long long ll;
const int maxn = 1e4;
ll V[maxn];
struct pers;
typedef pers* ppers;
struct pers {
    ppers l, r;
    ll cnt;
};
pers Pt[4 * maxn];
//build(1, 0, n - 1);  build Tree O(n*log(n))
void build(int v, int tl, int tr) {
    if (tl == tr) {
        Pt[v].cnt = V[tl];
        Pt[v].l = nullptr; Pt[v].r = nullptr;
    }
    else {
        int tm = (tl + tr) / 2;
        build(v * 2, tl, tm);
        build(v * 2 + 1, tm + 1, tr);
        Pt[v].cnt = Pt[v * 2].cnt + Pt[v * 2 + 1].cnt;
        Pt[v].l = &Pt[v * 2];
        Pt[v].r = &Pt[v * 2 + 1];
    }
}






// Центройды
// нам лень явно удалять вершины: заведем массив used -- была ли вершина удалена 
const int maxn = 1e3;
bool used[maxn];
int s[maxn]; // размеры поддеревьев 
bitset<maxn>g[maxn];// матрица смежности 

void sizes(int v, int p) {
    s[v] = 1;
    for (int u = 0; u < g[v].size(); u++) {
        if (u != p && !used[u]) {
            sizes(u, v), s[v] += s[u];
        }
    }
}

int centroid(int v, int p, int n) {
    for (int u = 0; u < g[v].size(); u++) {
        if (u != p && !used[u] && s[u] > n / 2) {
            return centroid(u, v, n);
        }
    }
    return v;
}

void solve(int v) {
    used[v] = true;
    for (int u = 0; u < g[v].size(); u++) {
        if (!used[u]) {
            solve(centroid(u, v, s[u]));
        }
    }
}




//Хеш таблица
typedef long long ll;

//template for generic type 
template<typename K, typename V>

//Hashnode class 
class HashNode
{
public:
    V value;
    K key;

    //Constructor of hashnode  
    HashNode(K key, V value) {
        this->value = value;
        this->key = key;
    }
};

//Our own Hashmap class 
class HashMap
{
    //hash element array 
    HashNode<K, V>** arr;
    int capacity;
    //current size 
    int size;
    //dummy node 
    HashNode<K, V>* dummy;

public:
    HashMap()
    {
        //Initial capacity of hash array 
        capacity = 20;
        size = 0;
        arr = new HashNode<K, V> * [capacity];

        //Initialise all elements of array as NULL 
        for (int i = 0; i < capacity; i++)
            arr[i] = NULL;

        //dummy node with value and key -1 
        dummy = new HashNode<K, V>(-1, -1);
    }
    // This implements hash function to find index 
    // for a key 
    int hashCode(K key)
    {
        return key % capacity;
    }

    //Function to add key value pair 
    void insertNode(K key, V value)
    {
        HashNode<K, V>* temp = new HashNode<K, V>(key, value);

        // Apply hash function to find index for given key 
        int hashIndex = hashCode(key);

        //find next free space  
        while (arr[hashIndex] != NULL && arr[hashIndex]->key != key
            && arr[hashIndex]->key != -1)
        {
            hashIndex++;
            hashIndex %= capacity;
        }

        //if new node to be inserted increase the current size 
        if (arr[hashIndex] == NULL || arr[hashIndex]->key == -1)
            size++;
        arr[hashIndex] = temp;
    }

    //Function to delete a key value pair 
    V deleteNode(int key)
    {
        // Apply hash function to find index for given key 
        int hashIndex = hashCode(key);

        //finding the node with given key 
        while (arr[hashIndex] != NULL)
        {
            //if node found 
            if (arr[hashIndex]->key == key)
            {
                HashNode<K, V>* temp = arr[hashIndex];

                //Insert dummy node here for further use 
                arr[hashIndex] = dummy;

                // Reduce size 
                size--;
                return temp->value;
            }
            hashIndex++;
            hashIndex %= capacity;

        }

        //If not found return null 
        return NULL;
    }

    //Function to search the value for a given key 
    V get(int key)
    {
        // Apply hash function to find index for given key 
        int hashIndex = hashCode(key);
        int counter = 0;
        //finding the node with given key    
        while (arr[hashIndex] != NULL)
        {
            int counter = 0;
            if (counter++ > capacity)  //to avoid infinite loop 
                return NULL;
            //if node found return its value 
            if (arr[hashIndex]->key == key)
                return arr[hashIndex]->value;
            hashIndex++;
            hashIndex %= capacity;
        }

        //If not found return null 
        return NULL;
    }

    //Return current size  
    int sizeofMap()
    {
        return size;
    }

    //Return true if size is 0 
    bool isEmpty()
    {
        return size == 0;
    }

    //Function to display the stored key value pairs 
    void display()
    {
        for (int i = 0; i < capacity; i++)
        {
            if (arr[i] != NULL && arr[i]->key != -1)
                cout << "key = " << arr[i]->key
                << "  value = " << arr[i]->value << endl;
        }
    }
};

//Driver method to test map class 
int main()
{
    HashMap<int, int>* h = new HashMap<int, int>;
    h->insertNode(1, 1);
    h->insertNode(2, 2);
    h->insertNode(2, 3);
    h->display();
    cout << h->sizeofMap() << endl;
    cout << h->deleteNode(2) << endl;
    cout << h->sizeofMap() << endl;
    cout << h->isEmpty() << endl;
    cout << h->get(2);

    return 0;
}


const int MOD = 1e6; // Вариант 1
struct Hash {
    vector<pair<int, int>> f[MOD];
    int size;
    int& operator [](int val) {
        for (int i = 0; i < f[val % MOD].size(); i++) {
            if (f[val % MOD][i].first == val) {
                return f[val % MOD][i].second;
            }
        }
        f[val % MOD].push_back(make_pair(val, 0));
        return f[val % MOD][f[val % MOD].size() - 1].second;
    }
};


struct hm { // Вариант 2
    int n;
    int* t;
    ll* ind;
    void gethm(int _n) {// берем хеш по основанию n
        n = _n;
        t = new int[n];
        ind = new ll[n];
        for (int i = 0; i < n; ++i) {
            ind[i] = 0;
        }
    }
    int& operator[](ll x) {
        int i = ((x % n) + n) % n; // хеш функция
        while (ind[i] != x && ind[i] != 0) {
            i = (i + 1) % n;
        }
        if (ind[i] == 0) {
            ind[i] = x;
        }
        return t[i];
    }
};



// Карась. Магическим образом exists(s) понимает была ли такая строка добавлена в add
struct Node;
typedef Node* PNode;
struct Node {
    vector <PNode> ch = vector <PNode>(95, nullptr);
    PNode suf = nullptr, p = nullptr, super_suf = nullptr;
    int c = -1;
    bool is_term = false;
};

PNode Proot = new Node();
void add(const string& s) {
    PNode cur = Proot;
    for (int i = 0; i < s.size(); ++i) {
        int c = s[i] - 32;
        if (!cur->ch[c]) {
            cur->ch[c] = new Node();
            cur->ch[c]->c = c;
            cur->ch[c]->p = cur;
        }
        cur = cur->ch[c];
    }
    cur->is_term = true;
}

PNode go_suf(PNode n);

PNode go(PNode n, int c) {
    if (n->ch[c] != nullptr) {
        return n->ch[c];
    }
    return go(go_suf(n), c);
}

PNode go_suf(PNode n) {
    if (n->suf != nullptr) {
        return n->suf;
    }
    n->suf = go(go_suf(n->p), n->c);
    return n->suf;
}

PNode calc_super_suf(PNode t) {
    if (t == Proot) {
        return t;
    }
    if (t->is_term) {
        return t;
    }
    if (t->super_suf == nullptr) {
        t->super_suf = calc_super_suf(go_suf(t));
    }
    return t->super_suf;
}

bool exists(string& s) {
    PNode cur = Proot;
    for (int i = 0; i < s.size(); ++i) {
        cur = go(cur, s[i] - 32);
        if (cur->super_suf == nullptr) {
            cur->super_suf = calc_super_suf(go_suf(cur));
        }
        if (cur->is_term || cur->super_suf->is_term) {
            return true;
        }
    }
    return 0;
}

void makeFake() {
    Proot->suf = new Node();
    Proot->p = Proot->suf;
    Proot->super_suf = Proot->suf;
    for (int i = 0; i < Proot->suf->ch.size(); ++i) {
        Proot->suf->ch[i] = Proot;
    }
}

signed main() {
    ifstream cin("input.txt");
    int n;
    cin >> n;
    makeFake();
    string s;
    getline(cin, s);
    for (int i = 0; i < n; ++i) {
        getline(cin, s);
        add(s);
    }
    while (getline(cin, s)) {
        if (exists(s)) {
            cout << s << "\n";
        }
    }
    return 0;
}



//Бор — структура данных для хранения набора строк,\
    представляющая из себя подвешенное дерево с символами на рёбрах. \
    Строки получаются последовательной записью всех символов, \
    хранящихся на рёбрах между корнем бора и терминальной вершиной.\
    Размер бора линейно зависит от суммы длин всех строк, а поиск в бору занимает время,\
    пропорциональное длине образца.
struct myBor;
typedef myBor* pBor;

struct myBor {
    vector<pair<myBor, char>> next;
    void push(string s) {
        pBor v = this;
        for (int i = 0; i < s.size(); i++) {
            bool find = false;
            for (int j = 0; j < v->next.size(); j++) {
                if (v->next[j].second == s[i]) {
                    find = true;
                    v = &v->next[j].first;
                    break;
                }
            }
            if (!find) {
                myBor New;
                v->next.resize(v->next.size() + 1);
                int u = v->next.size() - 1;
                v->next[u].first = New;
                v->next[u].second = s[i];
                v = &v->next[u].first;
            }
        }
    }
};

// Корневая декомпозиция с операциями удаления, вставки в произвольное место и вычисления функции на отрезке.\
n - колво чисел, len = sqrt(n) - длина блоков, cnt = ⌈n/len⌉ с ОКРУГЛЕНИЕМ ВВЕРХ
vector<int> T, L, R, A, B;
int n, len, cnt;
#define pb(x) push_back(x)

void build() {
    for (int i = 0; i < n; i++) { // разделим массив A на блоки длины len 
        T[i] = i; // В массиве T хранится актуальный порядок блоков [0, cnt-1] 
        L[i] = i * len;//для каждого блока сохраним индекс самого левого и самого правого элемента в массивах L, R,
        R[i] = (i + 1) * len - 1;
    }
    for (int i = 0; i < n; i++) {//в каждом блоке заранее посчитаем необходимую операцию 
        B[i / len] = B[i / len] + A[i];// и запишем в B
    }
}

int createNewBlock(int l, int r) {//Вспомогательная функция для создания блока под split
    int result = 0;
    for (int i = l; i < r; i++) {
        result = result + A[i]; //перешитываем результат операции на отрезке
    }
    B.pb(result);
    L.pb(l);
    R.pb(r);
    cnt++;
    return cnt - 1;
}
int split(int x) {
    int ind = 0;
    for (int i = 0; i < T.size(); i++) {
        if (L[i] <= x && x <= R[i]) {
            ind = i;
        }
    }
    if (L[ind] == x || R[ind] == x || x < 0) { // x -- не подходит
        return 0;
    }

    int first = createNewBlock(L[ind], x),
        second = createNewBlock(x, R[ind]);
    T.erase(T.begin() + ind); // операций T.erase(x) удаляет элемент под номером x и сдвигает массив T.
    T.insert(T.begin() + ind, first); // операций T.insert(x, y) вставляет в массив T после индекса x значение y и сдвигает массив. Время работы  
    T.insert(T.begin() + ind + 1, second);
    return ind;
}
//Асимптотика: O(cnt) – поиск нужного блока, O(len) – подсчет функции.Итог : O(len + cnt).Также мы увеличили cnt на 2.

//rebuild  . После операции split мы увеличили cnt на 2, а время работы всех функция зависит от него.\
    Для того чтоб cnt не стало слишком большим будем полностью перестраивать структуру\
    изменяя cnt на базовое значение равное 
void rebuild() {
    vector<int> tempA; // временная актуальная копия массива  
    for (int i = 0; i < T.size(); i++) {
        for (int j = L[i]; j < R[i]; j++) {
            tempA.pb(A[j]);
        }
    }
    A = tempA;
    B.clear(), L.clear(), R.clear();// очистка
    build(); // перестройка
}//O(A.size() + cnt)

int get(int l, int r) {
    int result = 0;
    int indexL = split(l - 1), //Разделим наши блоки при помощи операции split
        indexR = split(r) - 1;
    for (int i = indexL; i < indexR; i++) {
        result = result + B[T[i]]; //Посчитаем операцию на целых блоках использую массив B
    }
    return result;
}//O(cnt + len)

void erase(int x) { // erase A[x]
    split(x - 1);
    int ind = split(x);
    T.erase(T.begin() + ind);
}//O(cnt + len)

void insert(int x, int y) {
    int ind = split(x);
    A.pb(y);
    int indexNewBlock = createNewBlock(A.size() - 1, A.size() - 1);
    T.insert(T.begin() + ind, indexNewBlock);
}//O(cnt + len)
#include <map>
#include <iostream>
#define ll unsigned long long

using namespace std;

const int N = 3;
// N = 5 时 25 条边 
// 0: unknown
// 1: null edge
// 2: real edge
// const int E = 25;
// const int U = 14400; // (5!)^2

int cnt = 0;
ll edge[N * N];
int adj[N * N * N * N];

struct ListNode {
    ll val;
    ListNode *next;
    ListNode() : val(-1), next(nullptr) {}  // 默认值为-1
    ListNode(int x) : val(x), next(nullptr) {}  // 按值 x 初始化
    ListNode(int x, ListNode *next) : val(x), next(next) {}  // 使用新的节点初始化
};

struct ListHash {
    ll graph;
    double val;
    ListHash *next;
    ListHash() : graph(-1), next(nullptr) {}  // 默认值为-1
    ListHash(int x) : graph(x), next(nullptr) {}  // 按值 x 初始化
    ListHash(int x, ListHash *next) : graph(x), next(next) {}  // 使用新的节点初始化
};

unordered_map<ll, ListHash*> dp;

void print_graph(ll graph)
{
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            int e = i * N + j;
            cout << (graph / edge[e]) % 3 << ' ';
        }
        cout << endl;
    }
}

namespace Isomorph {	
	ll value(int type, ll val) {
		if (type == 0)
			return 3 * val * val + 998244353;

		if (type == 2)
			return 6 * val * val * val + 19260817;

		return 9 * val - 84033703;
	}
    void _hash_points(ll graph, ll vi[], ll vj[])
    {   
        ll ui[N], uj[N];
		for (int iter = 1; iter <= 2 * N; iter ++) {
			for (int i = 0; i < N; i ++) {
				ui[i] = 9999;
				for (int j = 0; j < N; j ++) {
                    int e = i * N + j;
                    int a = (graph / edge[e]) % 3;
					ui[i] += value(a, vj[j]);
                }
			}
			for (int j = 0; j < N; j ++) {
				uj[j] = 471231;
				for (int i = 0; i < N; i ++) {
                    int e = i * N + j;
                    int a = (graph / edge[e]) % 3;
					uj[j] += value(a, vi[i]);
                }
			}
			for (int i = 0; i < N; i ++)
				vi[i] = ui[i];
			for (int j = 0; j < N; j ++)
				vj[j] = uj[j];
		}
    }
	ll hash(ll graph) {
        ll vi[N], vj[N];
        fill(vi, vi + N, 84986076);
		fill(vj, vj + N, 12345678);
        _hash_points(graph, vi, vj);
		ll mul = 1;
		for (int i = 0; i < N; i ++)
			mul *= vi[i];
		for (int j = 0; j < N; j ++)
			mul *= vj[j];
		for (int i = 0; i < N; i ++)
			mul += 2183 * vi[i] * vi[i];
		for (int j = 0; j < N; j ++)
			mul += 5189 * vj[j] + 6;
		return mul;
	}
    bool check_isomorph(ll graph_a, ll graph_b)
    {
        int num = 0;
        ll vi_a[N], vj_a[N], vi_b[N], vj_b[N];
        fill(vi_a, vi_a + N, 84986076);
        fill(vj_a, vj_a + N, 12345678);
        fill(vi_b, vi_b + N, 84986076);
		fill(vj_b, vj_b + N, 12345678);
        _hash_points(graph_a, vi_a, vj_a);
        _hash_points(graph_b, vi_b, vj_b);
        int L[N], R[N];
        for (int i = 0; i < N; ++i) L[i] = R[i] = i;
        do {
            bool flag = true;
            for (int i = 0; i < N; ++i)
                if (vj_a[R[i]] != vj_b[i]) { flag = false; break; } 
            if (!flag) continue;
            do {
                flag = true;
                for (int i = 0; i < N; ++i) if (vi_a[L[i]] != vi_b[i]) { flag = false; break; } 
                if (!flag) continue;
                for (int i = 0; i < N; ++i) {
                    for (int j = 0; j < N; ++j) {
                        int e_b = i * N + j;
                        int e_a = L[i] * N + R[j];
                        if ((graph_a / edge[e_a]) % 3 != (graph_b / edge[e_b]) % 3) { flag = false; break; } 
                    }
                }
                if (flag) return true;
            } while(next_permutation(L, L + N));
        } while(next_permutation(R, R + N));

        // print_graph(graph_a);
        // for (int i = 0; i < N; ++i) cout << vi_a[i] << ' ';
        // cout << endl;
        // for (int i = 0; i < N; ++i) cout << vj_a[i] << ' ';
        // cout << endl;
        // cout << endl;

        // print_graph(graph_b);
        // for (int i = 0; i < N; ++i) cout << vi_b[i] << ' ';
        // cout << endl;
        // for (int i = 0; i < N; ++i) cout << vj_b[i] << ' ';
        // cout << endl;
        // cout << endl;

        return false;
    }
}

int compute_list_length(ListNode* head)
{
    int num = 0;
    ListNode *p = head;
    while (p->val != -1) {
        num += 1;
        p = p->next;
    }
    return num;
}

bool valid_query(ll graph, int i)
{
    // 如果 edge_i 已经 query 过，则不能再问
    if ((graph / edge[i]) % 3 != 0) return false;
    for (int j = 0; j < 2 * (N - 1); ++j) {
        // 检查 edge_i 所有邻接边
        int e = adj[i * 2 * (N - 1) + j];
        // 如果邻边存在且被问到了
        if ((graph / edge[e]) % 3 == 2) return false;
    }
    return true;
}

bool check_valid(ll graph, ll cur)
{
    for (int i = 0; i < N * N; ++i) {
        if (cur % 3 != 0 && cur % 3 != graph % 3) return false;
        graph /= 3, cur /= 3;
    }
    return true;
}

ListNode* compute_valid_graph(ll cur, ListNode* pre_head)
{
    ListNode *head, *tail, *p;
    p = pre_head, tail = head = new ListNode;
    while (p->val != -1) {
        if (check_valid(p->val, cur)) {
            tail->val = p->val;
            ListNode *node = new ListNode;
            tail->next = node;
            tail = node;
        }
        p = p->next;
    }
    return head;
}

double dfs(ll graph, ListNode* pre_head)
{
    ll cur = Isomorph::hash(graph);
    if (dp.find(cur) != dp.end()) {
        ListHash* p = dp[cur];
        while (p->next != nullptr) {
            if (Isomorph::check_isomorph(graph, p->graph)) return p->val;
            p = p->next;
        }
    }
    // ++cnt;
    // if (cnt % 1000000 == 0) cout <<  cnt / 1000000 << endl;

    double universal_size = (double)1.0 * compute_list_length(pre_head);
    // 枚举当前选哪条边, sum 记录当前状态收益
    double sum = 0.;
    for (int i = 0; i < N * N; ++i) {
        // 检查是否可以 query edge_i
        if (!valid_query(graph, i)) continue;
        ll status_i_exists = graph + 2 * edge[i];
        ll status_i_null = graph + edge[i];

        ListNode *head_exists = compute_valid_graph(status_i_exists, pre_head);
        ListNode *head_null = compute_valid_graph(status_i_null, pre_head);
        double exists_size = (double)1.0 * compute_list_length(head_exists);
        double null_size = (double)1.0 * compute_list_length(head_null);
        double p_i_exists = exists_size / universal_size;
        double p_i_null = null_size / universal_size;
        // debug
        if (p_i_exists + p_i_null < 0.999) {
            cout << cur << endl;
            cout << graph << endl;
            cout << p_i_exists << endl;
            cout << p_i_null << endl;
            cout << "error!!!" << endl;
            return 0;
        }
        double cur_sum = 0;
        if (p_i_exists > 0) cur_sum += p_i_exists * (dfs(status_i_exists, head_exists) + 1);
        if (p_i_null > 0) cur_sum += p_i_null * dfs(status_i_null, head_null);
        sum = max(sum, cur_sum);
    }
    // if (sum != 0) dp[cur] = sum;
    if (dp.find(cur) != dp.end()) {
        ListHash* p = dp[cur];
        while (p->next != nullptr) p = p->next;
        ListHash* tail = new ListHash;
        tail->graph = graph;
        tail->val = sum;
        p->next = tail;
    } else {
        dp[cur] = new ListHash;
        dp[cur]->graph = graph;
        dp[cur]->val = sum;
    }
    return sum;
}

ListNode* preprocess()
{
    // 计算 edge_i 的位数
    ll cur = 1;
    for (int i = 0; i < N * N; ++i) {
        edge[i] = cur;
        cur *= 3;
    }
    // 计算 edge_i 的邻接边
    for (int i = 0; i < N * N; ++i) {
        for (int j = 0; j < (N - 1); ++j) {
            adj[i * 2 * (N - 1) + j] = (i % N + (j + 1)) % N + i / N * N;
            adj[i * 2 * (N - 1) + (N - 1) + j] = (i / N + (j + 1)) % N * N + i % N;
        }
    }
    // 保存所有的上三角矩阵
    ListNode *head = new ListNode;
    ListNode *tail = head;
    int L[N], R[N], cnt = 0;
    for (int i = 0; i < N; ++i) {
        L[i] = i;
        R[i] = i; 
    }
    do {
        do {
            cur = 0;
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    cur += ((L[i] + R[j]) >= N - 1)? 2 * edge[i * N + j]: edge[i * N + j];
                }
            }
            // 挂链表
            ListNode *node = new ListNode;
            tail->val = cur;
            tail->next = node;
            tail = node;
            // cnt += 1;
            // cout << cnt << ':' << cur << endl;
        } while (next_permutation(L, L + N));
    } while (next_permutation(R, R + N));

    // print_upper_triangle_graph(cnt);
    return head;
}

int main()
{   
    clock_t start,end;
    start = clock();
    ListNode *head = preprocess();
    double Ans = dfs(0, head);
    cout << Ans / N << endl;
    cout << dp.size() << endl;
    end = clock();
    cout<<(double)(end-start)/CLOCKS_PER_SEC;
    return 0;
}
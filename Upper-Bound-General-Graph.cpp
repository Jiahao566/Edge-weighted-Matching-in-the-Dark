#include <map>
#include <iostream>
#define ll unsigned long long

using namespace std;

const int N = 4;
// 0: unknown
// 1: null edge
// 2: real edge

ll edge[(2 * N - 1) * N];
int adj[(2 * N - 1) * N * 2 * (2 * N - 2)];

struct ListNode {
    ll val;
    ListNode *next;
    ListNode() : val(-1), next(nullptr) {}  // default: -1
    ListNode(int x) : val(x), next(nullptr) {}  
};

struct ListHash {
    ll graph;
    double val;
    ListHash *next;
    ListHash() : graph(-1), next(nullptr) {}  // default: -1
    ListHash(int x) : graph(x), next(nullptr) {}  
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

ll edge_label(int i, int j)
{
    if (i > j) swap(i, j);
    return (4 * N - 1 - i) * i / 2 + j - i - 1;
}

int left_vertex_label(ll e)
{
    for (int i = 0; i < 2 * N - 1; ++i)
    for (int j = i + 1; j < 2 * N; ++j)
        if (edge_label(i, j) == e) return i;
}

int right_vertex_label(ll e)
{
    for (int i = 0; i < 2 * N - 1; ++i)
    for (int j = i + 1; j < 2 * N; ++j)
        if (edge_label(i, j) == e) return j;
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
        ll ui[2 * N], uj[2 * N];
		for (int iter = 1; iter <= 2 * N; iter ++) {
			for (int i = 0; i < 2 * N; i ++) {
				ui[i] = 9999;
				for (int j = 0; j < 2 * N; j ++) {
                    if (i == j) continue;
                    int e = edge_label(i, j);
                    int a = (graph / edge[e]) % 3;
					ui[i] += value(a, vj[j]);
                }
			}
			for (int j = 0; j < 2 * N; j ++) {
				uj[j] = 471231;
				for (int i = 0; i < 2 * N; i ++) {
                    if (i == j) continue;
                    int e = edge_label(i, j);
                    int a = (graph / edge[e]) % 3;
					uj[j] += value(a, vi[i]);
                }
			}
			for (int i = 0; i < 2 * N; i ++)
				vi[i] = ui[i];
			for (int j = 0; j < 2 * N; j ++)
				vj[j] = uj[j];
		}
    }
	ll hash(ll graph) {
        ll vi[2 * N], vj[2 * N];
        fill(vi, vi + 2 * N, 84986076);
		fill(vj, vj + 2 * N, 12345678);
        _hash_points(graph, vi, vj);
		ll mul = 1;
		for (int i = 0; i < 2 * N; i ++)
			mul *= vi[i];
		for (int j = 0; j < 2 * N; j ++)
			mul *= vj[j];
		for (int i = 0; i < 2 * N; i ++)
			mul += 2183 * vi[i] * vi[i];
		for (int j = 0; j < 2 * N; j ++)
			mul += 5189 * vj[j] + 6;
		return mul;
	}
    bool check_isomorph(ll graph_a, ll graph_b)
    {
        int num = 0;
        ll vi_a[2 * N], vj_a[2 * N], vi_b[2 * N], vj_b[2 * N];
        fill(vi_a, vi_a + 2 * N, 84986076);
        fill(vj_a, vj_a + 2 * N, 12345678);
        fill(vi_b, vi_b + 2 * N, 84986076);
		fill(vj_b, vj_b + 2 * N, 12345678);
        _hash_points(graph_a, vi_a, vj_a);
        _hash_points(graph_b, vi_b, vj_b);
        int V[2 * N];
        for (int i = 0; i < 2 * N; ++i) V[i] = i;
        sort(vi_a, vi_a + 2 * N);
        sort(vi_b, vi_b + 2 * N);
        for (int i = 0; i < 2 * N; ++i) if (vi_a[i] != vi_b[i]) return false;
        do {
            bool flag = true;
            for (int i = 0; i < 2 * N; ++i)
                if (vj_a[V[i]] != vj_b[i]) { flag = false; break; } 
            if (!flag) continue;
            for (int i = 0; i < 2 * N; ++i) {
                for (int j = 0; j < 2 * N; ++j) {
                    if (i == j) continue;
                    int e_b = edge_label(i, j);
                    int e_a = edge_label(V[i], V[j]);
                    if ((graph_a / edge[e_a]) % 3 != (graph_b / edge[e_b]) % 3) { flag = false; break; } 
                }
            }
            if (flag) return true;
        } while(next_permutation(V, V + 2 * N));

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
    // If edge_i has been queried
    if ((graph / edge[i]) % 3 != 0) return false;
    for (int j = 0; j < 2 * (2 * N - 2); ++j) {
        // check all edges sharing common point with edge_i
        int e = adj[i * 2 * (2 * N - 2) + j];
        if ((graph / edge[e]) % 3 == 2) return false;
    }
    return true;
}

bool check_valid(ll graph, ll cur)
{
    for (int i = 0; i < (2 * N - 1) * N; ++i) {
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

    double universal_size = (double)1.0 * compute_list_length(pre_head);
    //Enumerate chosen edges, sum records the expected rewards
    double sum = 0.;
    for (int i = 0; i < (2 * N - 1) * N; ++i) {
        // Check whethr querying edge_i is legal
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
    // Label edges, use a ternary number to represent a graph
    ll cur = 1;
    for (int i = 0; i < (2 * N - 1) * N; ++i) {
        edge[i] = cur;
        cur *= 3;
    }
    // compute edges sharing a common point with edge_i
    for (int i = 0; i < (2 * N - 1) * N; ++i) {
        int u = left_vertex_label(i), v = right_vertex_label(i);
        // cout << i <<' ' << u << ' ' << v << endl;
        int cnt_u = 0, cnt_v = 0; 
        for (int j = 0; j < 2 * N - 1; ++j) {
            int new_v = (u + 1 + j) % (2 * N), new_u = (v + 1 + j) % (2 * N);
            if (new_v != u && new_v != v) adj[i * 2 * (2 * N - 2) + cnt_u++] = edge_label(u, new_v);
            if (new_u != u && new_u != v) adj[i * 2 * (2 * N - 2) + (2 * N - 2) + cnt_v++] = edge_label(new_u, v);
        }
    }

    // Compute all possible hard-instnce matrices
    ListNode *head = new ListNode;
    ListNode *tail = head;
    int degree[2 * N], cnt = 0;
    for (int i = 0; i < 2 * N; ++i) degree[i] = i;
    do {
        cur = 0;
        for (int i = 0; i < 2 * N - 1; ++i) {
            for (int j = i + 1; j < 2 * N; ++j) {
                ll e = edge_label(i, j);
                int degree_i = degree[i] < N ? degree[i] + 1: degree[i];
                int degree_j = degree[j] < N ? degree[j] + 1: degree[j];
                cur += ((degree_i + degree_j) >= 2 * N)? 2 * edge[e]: edge[e];
            }
        }
        // linked list
        ListNode *node = new ListNode;
        tail->val = cur;
        tail->next = node;
        tail = node;
    } while (next_permutation(degree, degree + 2 * N));

    return head;
}

int main()
{   
    clock_t start,end;
    start = clock();
    ListNode *head = preprocess();
    double Ans = dfs(0, head);
    cout << "Ratio: " << Ans / N << endl;
    cout << "Number of DP states: " << dp.size() << endl;
    end = clock();
    cout<< "Running time: " << (double)(end-start)/CLOCKS_PER_SEC;
    return 0;
}
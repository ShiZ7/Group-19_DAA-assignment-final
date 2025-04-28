#include <iostream>
#include <vector>
#include <algorithm>
#include <queue>
#include <cmath>
#include <map>
#include <set>
#include <fstream>
#include <unordered_map>
#include <sstream>
#include <chrono>

using namespace std;
using namespace chrono;

// Edge structure for max-flow
struct Edge {
    int to, rev;
    double cap;
    Edge(int t, int r, double c) : to(t), rev(r), cap(c) {}
};

// MaxFlow class using Dinic's algorithm
class MaxFlow {
public:
    MaxFlow(int n) : adj(n), level(n), ptr(n) {}
    void addEdge(int from, int to, double cap) {
        adj[from].emplace_back(to, adj[to].size(), cap);
        adj[to].emplace_back(from, adj[from].size() - 1, 0.0);
    }
    double dinic(int s, int t) {
        double flow = 0.0;
        while (bfs(s, t)) {
            fill(ptr.begin(), ptr.end(), 0);
            while (double f = dfs(s, t, 1e18))
                flow += f;
        }
        return flow;
    }
    vector<int> get_cut(int s) {
        vector<int> cut;
        vector<bool> visited(adj.size(), false);
        queue<int> q;
        q.push(s);
        visited[s] = true;
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            if (u >= 2 && u < 2 + adj.size()) cut.push_back(u - 2);
            for (auto &e : adj[u]) {
                if (e.cap > 1e-9 && !visited[e.to]) {
                    visited[e.to] = true;
                    q.push(e.to);
                }
            }
        }
        return cut;
    }
private:
    bool bfs(int s, int t) {
        fill(level.begin(), level.end(), -1);
        queue<int> q;
        level[s] = 0;
        q.push(s);
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            for (auto &e : adj[u]) {
                if (e.cap > 1e-9 && level[e.to] == -1) {
                    level[e.to] = level[u] + 1;
                    q.push(e.to);
                    if (e.to == t) return true;
                }
            }
        }
        return false;
    }
    double dfs(int u, int t, double flow) {
        if (u == t) return flow;
        for (; ptr[u] < adj[u].size(); ptr[u]++) {
            auto &e = adj[u][ptr[u]];
            if (e.cap > 1e-9 && level[e.to] == level[u] + 1) {
                double f = dfs(e.to, t, min(flow, e.cap));
                if (f > 1e-9) {
                    e.cap -= f;
                    adj[e.to][e.rev].cap += f;
                    return f;
                }
            }
        }
        return 0.0;
    }
public:
    vector<vector<Edge>> adj;
    vector<int> level, ptr;
};

// Global variables
vector<vector<int>> adj;
int n, m, h;
vector<int> clique_degree;
vector<vector<int>> h_minus_1_cliques;
vector<vector<int>> h_cliques;

// Check if a set of vertices forms a clique
bool is_clique(const vector<int>& vertices) {
    for (int i = 0; i < vertices.size(); i++) {
        for (int j = i + 1; j < vertices.size(); j++) {
            if (!binary_search(adj[vertices[i]].begin(), adj[vertices[i]].end(), vertices[j])) {
                return false;
            }
        }
    }
    return true;
}

// Enumerate k-cliques
void enumerate_k_cliques(int k, vector<vector<int>>& cliques) {
    cliques.clear();
    if (k == 1) {
        for (int v = 0; v < n; v++) cliques.push_back({v});
    } else if (k == 2) {
        for (int v1 = 0; v1 < n; v1++) {
            for (int v2 = v1 + 1; v2 < n; v2++) {
                if (is_clique({v1, v2})) cliques.push_back({v1, v2});
            }
        }
    } else {
        vector<vector<int>> smaller;
        enumerate_k_cliques(k - 1, smaller);
        for (auto &cl : smaller) {
            for (int v = cl.back() + 1; v < n; v++) {
                if (all_of(cl.begin(), cl.end(), [&](int u) { return binary_search(adj[u].begin(), adj[u].end(), v); })) {
                    auto new_cl = cl;
                    new_cl.push_back(v);
                    cliques.push_back(new_cl);
                }
            }
        }
    }
}

// Generate cliques
void generate_cliques() {
    enumerate_k_cliques(h - 1, h_minus_1_cliques);
    enumerate_k_cliques(h, h_cliques);
    clique_degree.assign(n, 0);
    for (const auto& clique : h_cliques) {
        for (int v : clique) clique_degree[v]++;
    }
}

// Core decomposition
vector<int> core_decomposition() {
    vector<int> deg(n), core(n), removed(n, 0);
    for (int v = 0; v < n; v++) deg[v] = adj[v].size();
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<>> pq;
    for (int v = 0; v < n; v++) pq.emplace(deg[v], v);
    int k = 0;
    while (!pq.empty()) {
        auto [d, v] = pq.top();
        pq.pop();
        if (removed[v]) continue;
        core[v] = k;
        removed[v] = 1;
        for (int u : adj[v]) {
            if (!removed[u]) {
                deg[u]--;
                pq.emplace(deg[u], u);
            }
        }
        k = max(k, d);
    }
    return core;
}

// Find connected components
vector<vector<int>> find_connected_components(const vector<int>& vertices) {
    vector<vector<int>> components;
    vector<bool> visited(n, false);
    for (int v : vertices) {
        if (!visited[v]) {
            vector<int> comp;
            queue<int> q;
            q.push(v);
            visited[v] = true;
            while (!q.empty()) {
                int u = q.front();
                q.pop();
                comp.push_back(u);
                for (int w : adj[u]) {
                    if (!visited[w] && find(vertices.begin(), vertices.end(), w) != vertices.end()) {
                        visited[w] = true;
                        q.push(w);
                    }
                }
            }
            components.push_back(comp);
        }
    }
    return components;
}

// Load dataset
void load_as733(const string& filename, vector<vector<int>>& adj_list, int& n, int& m) {
    ifstream infile(filename);
    if (!infile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        exit(1);
    }
    vector<pair<int, int>> edges;
    int u, v, max_vertex = 0;
    while (infile >> u >> v) {
        edges.emplace_back(u, v);
        max_vertex = max(max_vertex, max(u, v));
    }
    infile.close();
    n = max_vertex + 1;
    m = edges.size();
    adj_list.assign(n, {});
    for (const auto& edge : edges) {
        u = edge.first;
        v = edge.second;
        adj_list[u].push_back(v);
        adj_list[v].push_back(u);
    }
    for (auto& neighbors : adj_list) {
        sort(neighbors.begin(), neighbors.end());
        neighbors.erase(unique(neighbors.begin(), neighbors.end()), neighbors.end());
    }
}

int main() {
    string filename = "as733.txt";
    cout << "Loading As733 network...\n";
    load_as733(filename, adj, n, m);
    cout << "Vertices: " << n << " Edges: " << m << endl;

    ofstream density_out("density_vs_h.txt");

    for (h = 2; h <= 6; ++h) {
        cout << "\nProcessing h = " << h << "\n";

        auto start = high_resolution_clock::now();

        generate_cliques();
        cout << "Found " << h_cliques.size() << " cliques of size " << h << endl;
        cout << "Found " << h_minus_1_cliques.size() << " cliques of size " << (h - 1) << endl;

        vector<int> core = core_decomposition();
        int k_max = *max_element(core.begin(), core.end());

        double l = 0.0, u = k_max;
        vector<int> best_subgraph;
        double best_density = 0.0;

        vector<int> core_vertices;
        for (int v = 0; v < n; v++) {
            if (core[v] >= 1) core_vertices.push_back(v);
        }

        vector<vector<int>> components = find_connected_components(core_vertices);

        for (const auto& component : components) {
            if (component.size() < h) continue;

            double density = (double)h_cliques.size() / component.size();
            if (density > best_density) {
                best_density = density;
                best_subgraph = component;
            }
        }

        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<seconds>(stop - start);
        cout << "h = " << h << ", Best density = " << best_density << ", Time: " << duration.count() << " seconds\n";
        density_out << h << " " << best_density << endl;
    }

    density_out.close();
    return 0;
}

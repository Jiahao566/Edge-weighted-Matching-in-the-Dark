#include <vector>
#include <iostream>
#include <iomanip>
#include <functional>
#include "gurobi_c++.h"

using namespace std;

const int precision = 4;
const int max_iterations = 20;
const double mip_gap = 1e-4;
int n;
double eps;

// Generate the initial set of threshold pairs theta and beta_inv.
// Function theta have value n at all points in {0, 1, ..., n}.
// Function beta_inv is non-decreasing and from {0, 1, ..., n} to {0, 1, ..., n}.
vector<pair<vector<int>, vector<int>>> generate_initial_thresholds(int n)
{
    vector<pair<vector<int>, vector<int>>> res;
    vector<int> theta(n, n);
    vector<int> beta_inv(n);
    function<void(int)> dfs = [&](int i)
    {
        if (i == n)
            res.push_back({theta, beta_inv});
        else
            for (beta_inv[i] = (i > 0 ? beta_inv[i - 1] : 0); beta_inv[i] <= n; beta_inv[i]++)
                dfs(i + 1);
    };
    dfs(0);
    return res;
}

// Generate all possible pairs of threshold functions theta and beta_inv.
// Both functions are non-decreasing and from {0, 1, ..., n} to {0, 1, ..., n}.
// Function theta is pointwise greater than or equal to beta_inv.
// TODO: Halve the number of pairs by symmetry, ensuring sum(theta) >= sum(beta)
vector<pair<vector<int>, vector<int>>> generate_thresholds(int n)
{
    vector<pair<vector<int>, vector<int>>> res;
    vector<int> theta(n), beta_inv(n);
    function<void(int)> dfs = [&](int i)
    {
        if (i == n)
            res.push_back({theta, beta_inv});
        else
            for (theta[i] = (i > 0 ? theta[i - 1] : 0); theta[i] <= n; theta[i]++)
                for (beta_inv[i] = (i > 0 ? beta_inv[i - 1] : 0); beta_inv[i] <= theta[i]; beta_inv[i]++)
                    dfs(i + 1);
    };
    dfs(0);
    return res;
}

// Compute the inverse of non-decreasing function from {0, 1, ..., n} to {0, 1, ..., n}.
vector<int> inverse(vector<int> &arr)
{
    int n = arr.size();
    vector<int> res(n);
    int i = 0;

    for (int j = 0; j < n; j++)
        while (i < arr[j])
        {
            res[i] = j;
            i++;
        }

    while (i < n)
    {
        res[i] = n;
        i++;
    }

    return res;
}

// Solve the competitive ratio and the corresponding gain-sharing functions g and h,
// in the relaxed problem with competitiveness constraints for the given threshold pairs.
tuple<float, vector<double>, vector<double>> solve(vector<pair<vector<int>, vector<int>>> thresholds, tuple<double, vector<double>, vector<double>> init_x)
{
    double res_ratio;
    vector<double> res_g, res_h;

    try
    {
        // Create an environment
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "grid_function.log");
        env.start();

        // Create an empty model
        GRBModel model = GRBModel(env);

        model.set(GRB_DoubleParam_MIPGap, mip_gap);

        // Create variables
        GRBVar *g = model.addVars(n + 1, GRB_CONTINUOUS);
        GRBVar *h = model.addVars(n + 1, GRB_CONTINUOUS);
        for (int i = 0; i <= n; i++)
        {
            g[i].set(GRB_DoubleAttr_LB, 0.0);
            h[i].set(GRB_DoubleAttr_LB, 0.0);
        }
        GRBVar Gamma = model.addVar(0.0, 1.0, 1.0, GRB_CONTINUOUS, "Gamma");

        if (get<0>(init_x) > 0.0)
        {
            Gamma.set(GRB_DoubleAttr_Start, get<0>(init_x));
            for (int i = 0; i <= n; i++)
            {
                g[i].set(GRB_DoubleAttr_Start, get<1>(init_x)[i]);
                h[i].set(GRB_DoubleAttr_Start, get<2>(init_x)[i]);
            }
        }

        model.setObjective(1.0 * Gamma, GRB_MAXIMIZE);

        // Fixing g[0]
        // model.addConstr(g[0] == 0.82);

        // Gain sharing constraints
        for (int i = 0; i <= n; i++)
            for (int j = 0; j <= i; j++)
                model.addQConstr(g[i] * h[j] + h[i] * g[j] <= 1.0);

        // Monotonicity constraints
        for (int i = 0; i < n; i++)
        {
            model.addConstr(g[i] >= g[i + 1]);
            model.addConstr(h[i] <= h[i + 1]);
        }

        // Boundary constraint
        model.addQConstr(g[0] * g[0] + h[0] * h[0] == 1.0);
        model.addConstr(g[0] >= h[0]);
        model.addConstr(g[n] == 0.0);

        // Competitive ratio constraints
        for (auto &c : thresholds)
        {
            auto theta = c.first;
            auto beta_inv = c.second;
            auto beta = inverse(beta_inv);
            auto theta_inv = inverse(theta);

            GRBQuadExpr ratio = 0.0;
            for (int i = 0; i < n; i++)
                ratio += (theta[i] - beta_inv[i]);
            for (int i = 0; i < n; i++)
                ratio += (n - theta[i] + beta_inv[i]) * h[i] * g[theta[i]];
            for (int i = 0; i < n; i++)
                ratio += (n - beta[i] + theta_inv[i]) * h[i] * g[beta[i]];
            model.addQConstr(ratio >= n * n * Gamma);
        }

        // Optimize model
        model.optimize();

        // Print solution
        cout << "ratio = " << fixed << setprecision(precision) << Gamma.get(GRB_DoubleAttr_X) << endl;
        res_ratio = Gamma.get(GRB_DoubleAttr_X);

        cout << "g = [ ";
        for (int i = 0; i <= n; i++)
        {
            cout << fixed << setprecision(precision) << g[i].get(GRB_DoubleAttr_X) << " ";
            if (i < n)
                cout << ", ";
            res_g.push_back(g[i].get(GRB_DoubleAttr_X));
        }
        cout << "]" << endl;

        cout << "h = [ ";
        for (int i = 0; i <= n; i++)
        {
            cout << fixed << setprecision(precision) << h[i].get(GRB_DoubleAttr_X) << " ";
            if (i < n)
                cout << ", ";
            res_h.push_back(h[i].get(GRB_DoubleAttr_X));
        }
        cout << "]" << endl;
    }
    catch (GRBException e)
    {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }

    return {res_ratio, res_g, res_h};
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        cerr << "Usage: " << argv[0] << " <n>" << endl;
        return 1;
    }
    n = stoi(argv[1]);
    eps = 1.0 / n;

    auto thresholds = generate_initial_thresholds(n);
    auto test_thresholds = generate_thresholds(n);
    auto [ratio, g, h] = solve(thresholds, make_tuple(0.0, vector<double>(), vector<double>()));
    double margin = 0.1;

    for (int i = 0; i < max_iterations; i++)
    {
        cout << endl
             << "**** Iteration " << i + 1 << " ****"
             << endl
             << endl;

        double max_gap = 0.0;

        for (auto &c : test_thresholds)
        {
            auto theta = c.first;
            auto beta_inv = c.second;
            auto beta = inverse(beta_inv);
            auto theta_inv = inverse(theta);

            double test = 0.0;
            for (int i = 0; i < n; i++)
                test += (theta[i] - beta_inv[i]);
            for (int i = 0; i < n; i++)
                test += (n - theta[i] + beta_inv[i]) * h[i] * g[theta[i]];
            for (int i = 0; i < n; i++)
                test += (n - beta[i] + theta_inv[i]) * h[i] * g[beta[i]];
            double gap = ratio - test / (n * n);
            max_gap = max(max_gap, gap);
            if (gap > margin)
                thresholds.push_back(c);
        }

        cout << "Max gap = " << fixed << setprecision(precision + 2) << max_gap << endl;

        if (max_gap <= margin)
        {
            cout << "All constraints satisfied with margin " << fixed << setprecision(precision) << margin << endl;
            if (margin < 0.0002)
                break;
            else
                // margin /= 10;
                margin = max_gap / 2;
        }
        else
        {
            double init_ratio = ratio - max_gap - 1e-4;
            tie(ratio, g, h) = solve(thresholds, make_tuple(init_ratio, g, h));
        }
    }

    return 0;
}
#include <iostream>
#include <iomanip>
#include <vector>
#include <functional>
#include <algorithm>
#include <cmath>

#include <omp.h>

#include "progressbar.hpp"

using namespace std;

const int precision = 4;
const int chunk_size = 1;

const double c1_lin = 1.351;
const double c2_lin = 1 - c1_lin;
const double c1_arc = 0.890;
const double c2_arc = 1 - c1_arc;

inline double lin_ratio(double x)
{
    return c2_lin * x * x + c1_lin * x;
}

double arc_ratio(double x)
{
    return c2_arc * x * x + c1_arc * x;
}

// choice of n_lin and alpha for different n
// ------------------------------
//   n  | n_lin | alpha | ratio
// ------------------------------
//   5  |   4   | 0.741 | 0.6362
//   6  |   4   | 0.671 | 0.6436
//   7  |   4   | 0.626 | 0.6479
//   8  |   5   | 0.650 | 0.6506
//   9  |   5   | 0.615 | 0.6529
//  10  |   5   | 0.585 | 0.6549
//  10  |   6   | 0.640 | 0.6548
//  11  |   6   | 0.612 | 0.6565
//  11  |   7   | 0.657 | 0.6558
//  12  |   7   | 0.632 | 0.6575
//  13  |   7   | 0.609 | 0.6590
//  13  |   8   | 0.646 | 0.6580
//  14  |   ?
// ------------------------------
const int n = 14;
const int n_lin = 8;             // number of grid points in linear part
const int n_arc = n - n_lin + 1; // number of grid points in arc part
// ----------------
//  alpha | ratio
// ----------------
//  0.580 |
//  0.590 |
//  0.600 | 0.6557
//  0.605 | 0.6563
//  0.609 | 0.6564
//  0.610 | 0.6565
//  0.611 | 0.6565
//  0.612 | 0.6565
//  0.613 | 0.6565
//  0.614 | 0.6565
//  0.615 | 0.6563
//  0.620 | 0.6557
//  0.625 |
//  0.630 |
//  0.635 |
//  0.640 |
//  0.645 |
//  0.650 |
// ----------------
// const float alpha = 0.612;
// double alpha = 0.612;

const double eps = 1.0 / n;
// const int k = min(7, n); // cutoff index for parallelization
const int k = min(7, n - 6);

array<float, n + 1> g, h;
// float alpha;

void print(array<float, n + 1> &arr)
{
    for (auto &a : arr)
    {
        cout << fixed << setprecision(precision) << a << " ";
    }
    cout << endl;
}

vector<pair<array<int, k>, array<int, k>>> gen_combinations()
{
    vector<pair<array<int, k>, array<int, k>>> res;
    array<int, k> theta, beta_inv;
    function<void(int)> dfs = [&](int i)
    {
        if (i == k)
            res.push_back({theta, beta_inv});
        else
            for (theta[i] = theta[i - 1]; theta[i] <= n; theta[i]++)
                for (beta_inv[i] = beta_inv[i - 1]; beta_inv[i] <= theta[i]; beta_inv[i]++)
                    dfs(i + 1);
    };
    dfs(0);
    return res;
}

array<int, n> inverse(array<int, n> &arr)
{
    array<int, n> res;
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

float check(array<int, n> &beta_inv, array<int, n> &theta)
{
    double ratio = 0.0;

    array<int, n> beta = inverse(beta_inv);
    array<int, n> theta_inv = inverse(theta);

    for (int i = 0; i < n; i++)
        ratio += (theta[i] - beta_inv[i]) * eps * eps;

    for (int i = 0; i < n; i++)
        ratio += (n - theta[i] + beta_inv[i]) * h[i] * g[theta[i]] * eps * eps;

    for (int i = 0; i < n; i++)
        ratio += (n - beta[i] + theta_inv[i]) * h[i] * g[beta[i]] * eps * eps;

    return ratio;
}

float check_in_parallel(array<int, k> theta_start, array<int, k> beta_inv_start)
{
    float ratio = 1.0;

    array<int, n> theta, beta_inv;
    for (int i = 0; i < k; i++)
    {
        theta[i] = theta_start[i];
        beta_inv[i] = beta_inv_start[i];
    }

    function<void(int)> dfs = [&](int i)
    {
        if (i == n)
        {
            float res = check(beta_inv, theta);
            if (res < ratio)
            {
                ratio = res;
            }
        }
        else
            for (theta[i] = theta[i - 1]; theta[i] <= n; theta[i]++)
                for (beta_inv[i] = beta_inv[i - 1]; beta_inv[i] <= theta[i]; beta_inv[i]++)
                    dfs(i + 1);
    };
    dfs(k);
    return ratio;
}

int main()
{
    double alpha;
    double alpha_lb = 0.4;
    double alpha_ub = M_PI / 4;
    while (alpha_ub - alpha_lb > 1e-3)
    {
        alpha = (alpha_lb + alpha_ub) / 2;
        double h_sep = cos(alpha);

        double gamma = alpha + (M_PI / 2 - 2 * alpha) * arc_ratio((float)(n_arc - 2) / (n_arc - 1));
        double h_arc = sin(gamma);

        double g_lin = sin(alpha) * lin_ratio((float)(n_lin - 1) / n_lin);
        double h_lin = 1.0 / cos(alpha) - tan(alpha) * g_lin;

        if ((h_sep / h_arc) < (h_lin / h_sep))
            alpha_ub = alpha;
        else
            alpha_lb = alpha;
    }
    alpha = (alpha_lb + alpha_ub) / 2 - 3e-3;
    // alpha = 0.64;

    for (int i = 0; i < n_arc; i++)
    {
        double gamma = alpha + (M_PI / 2 - 2 * alpha) * arc_ratio((float)i / (n_arc - 1));
        g[i] = cos(gamma);
        h[i] = sin(gamma);
    }

    for (int i = n_arc; i <= n; i++)
    {
        g[i] = sin(alpha) * lin_ratio((float)(n - i) / n_lin);
        h[i] = 1.0 / cos(alpha) - tan(alpha) * g[i];
    }

    cout << "n = " << n << endl;
    cout << "alpha = " << fixed << setprecision(4) << alpha << endl;
    cout << "g: ";
    print(g);
    cout << "h: ";
    print(h);

    cout << "Gain matrix:" << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << fixed << setprecision(precision) << g[i] * h[j] + g[j] * h[i] << " ";
        }
        cout << endl;
    }

    cout << "Generating combinations..." << endl;
    auto combinations = gen_combinations();

    int num_procs = omp_get_num_procs();

    cout << "Computing ratio in parallel with " << num_procs << " processors..." << endl;
    float ratio = 1.0;

    progressbar bar(combinations.size());

#pragma omp parallel for reduction(min : ratio) schedule(dynamic, chunk_size)
    for (auto &c : combinations)
    {
        float res = check_in_parallel(c.first, c.second);
        if (res < ratio)
            ratio = res;

#pragma omp critical
        bar.update();
    }
    cout << endl
         << "RATIO: " << fixed << setprecision(precision) << ratio << endl;

    return 0;
}
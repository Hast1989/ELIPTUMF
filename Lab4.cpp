#include <iostream>
#include <fstream>
#include<cmath>
#include <iomanip>
#include<string>
double pi = 3.1415926535898;
double eps,h1,h2,tau;
int  n, m,indP;
double NormV2(double** x, double** y, int n,int m)
{
    double res;
    res = 0;
    for (int i = 0; i < n+1; i++)
    {
        for (int j = 0; j < m+1; j++)
        {
            res =res+ std::fabs(x[i][j] - y[i][j]);
        }
    }
    return res;
}
double KSI10(double x1)
{
    if (indP == 0)
    {
        return 0;
    }
    if (indP == 1)
    {
        return 1+x1;
    }
    if (indP == 2)
    {
        return -0;
    }
    if (indP == 3)
    {
        return 0;
    }
    if (indP == 4)
    {
        return x1;
    }
    return 0;
}
double KSI11(double x1)
{
    if (indP == 0)
    {
        return 0;
    }
    if (indP == 1)
    {
        return 1+x1;
    }
    if (indP == 2)
    {
        return 2;
    }
    if (indP == 3)
    {
        return 0;
    }
    if (indP == 4)
    {
        return pi*std::sin(x1)-x1;
    }
    return 0;
}
double KSI20(double x2)
{
    if (indP == 0)
    {
        return 0;
    }
    if (indP == 1)
    {
        return 1+x2;
    }
    if (indP == 2)
    {
        return x2*x2;
    }
    if (indP == 3)
    {
        return 0;
    }
    if (indP == 4)
    {
        return x2;
    }
    return 0;
}
double KSI21(double x2)
{
    if (indP == 0)
    {
        return 0;
    }
    if (indP == 1)
    {
        return 1+x2;
    }
    if (indP == 2)
    {
        return 1+x2*x2;
    }
    if (indP == 3)
    {
        return 0;
    }
    if (indP == 4)
    {
        return pi*std::sin(x2)-x2;
    }
    return 0;
}
double PSI10(double x1)
{
    if (indP == 0)
    {
        return 0;
    }
    if (indP == 1)
    {
        return 1;
    }
    if (indP == 2)
    {
        return 0;
    }
    if (indP == 3)
    {
        return -(1-x1)*x1;
    }
    if (indP == 4)
    {
        return std::cos(x1)+x1;
    }
    return 0;
}
double PSI11(double x1)
{
    if (indP == 0)
    {
        return 0;
    }
    if (indP == 1)
    {
        return 1;
    }
    if (indP == 2)
    {
        return 2;
    }
    if (indP == 3)
    {
        return (1-x1)*x1;
    }
    if (indP == 4)
    {
        return std::cos(x1) - x1;
    }
    return 0;
}
double PSI20(double x2)
{
    if (indP == 0)
    {
        return 0;
    }
    if (indP == 1)
    {
        return 0;
    }
    if (indP == 2)
    {
        return 0;
    }
    if (indP == 3)
    {
        return 0;
    }
    if (indP == 4)
    {
        return 0;
    }
    return 0;
}
double PSI21(double x2)
{
    if (indP == 0)
    {
        return 0;
    }
    if (indP == 1)
    {
        return 0;
    }
    if (indP == 2)
    {
        return 2;
    }
    if (indP == 3)
    {
        return 0;
    }
    if (indP == 4)
    {
        return 0;
    }
    return 0;
}
double U0(double x1, double x2)
{
    if (indP == 0)
    {
        return 8.;
    }
    if (indP == 1)
    {
        return 3.;
    }
    if (indP == 2)
    {
        return x1*x1+x2*x2;
    }
    if (indP == 3)
    {
        return 0;
    }
    return 0;
}
double F(double x1, double x2)
{
    if (indP == 0)
    {
        return 0;
    }
    if (indP == 1)
    {
        return 0;
    }
    if (indP == 2)
    {
        return -4.;
    }
    if (indP == 3)
    {
        return -(2 * x2 + x1 - (x2 * x2 + x1 * x1));
    }
    if (indP == 4)
    {
        return x2*std::cos(x1)+x1*std::sin(x2);
    }
    return 0;
}
void Progonka(int m, int n, double* L, double* D, double* U, double* rightb, double* x)  //n - количество не граничных узлов
{
    double* alpha;
    double* betta;
    double res;
    alpha = new double[n];
    betta = new double[n];
    alpha[m] = -U[m] / D[m];
    betta[m] = rightb[m] / D[m];
    alpha[n - 1] = 0;
    betta[n - 1] = 0;
    for (int i = m + 1; i < n - 1; i++)
    {
        res = L[i] * alpha[i - 1] + D[i];
        alpha[i] = -U[i] / res;
        betta[i] = (rightb[i] - L[i] * betta[i - 1]) / res;
    }
    x[n - 1] = (rightb[n - 1] - L[n - 1] * betta[n - 2]) / (L[n - 1] * alpha[n - 2] + D[n - 1]);
    for (int i = n - 2; i >= m; i--)
    {
        x[i] = alpha[i] * x[i + 1] + betta[i];
    }
    delete[] alpha;
    delete[] betta;
}
void MStept(double** res0, double** res1, int n, int m)
{
    double** res;
    double** resy;
    res = new double* [m + 1];
    for (int i = 0; i < m + 1; i++)
    {
        res[i] = new double[n + 1];
    }
    resy = new double* [n + 1];
    for (int i = 0; i < n + 1; i++)
    {
        resy[i] = new double[m + 1];
    }
    double* L1;
    L1 = new double[n + 1];
    double* D1;
    D1 = new double[n + 1];
    double* U1;
    U1 = new double[n + 1];
    double* rightb1;
    rightb1 = new double[n + 1];
    double* L2;
    L2 = new double[m + 1];
    double* D2;
    D2 = new double[m + 1];
    double* U2;
    U2 = new double[m + 1];
    double* rightb2;
    rightb2 = new double[m + 1];
    for (int j = 0; j < n + 1; j++)
    {
        L1[j] = 1;
        D1[j] = -2 * (1. + (h1 * h1) / tau);
        U1[j] = L1[j];
    }
    for (int j = 0; j < m + 1; j++)
    {
        L2[j] = 1.;
        D2[j] = -2 * (1. + (h2 * h2) / tau);
        U2[j] = L2[j];
    }
    L1[0] = 0;
    D1[0] = -((h1 * h1) / tau + 1);
    U1[0] = 1;
    L1[n] = 1;
    D1[n] = -((h1 * h1) / tau + 1);
    U1[n] = 0;
    for (int j = 0; j < m + 1; j++)
    {
        res0[j][0] = KSI20(j*h2);
        res0[j][m] = KSI21(j*h2);
    }
    for (int i = 0; i < n + 1; i++)
    {
        for (int j = 0; j < m + 1; j++)
        {
            res[j][i] = res0[i][j];
        }
    }
    for (int j = 1; j < m; j++)
    {
        for (int i = 1; i < n; i++)
        {
            rightb1[i] = -(2 * (h1 * h1) / tau) * res0[i][j] - ((h1 * h1) * (res0[i][j + 1] - 2 * res0[i][j] + res0[i][j - 1])) / (h2 * h2) - (h1 * h1) * F(i * h1, j * h2);
        }
        rightb1[0] = -((h1 * h1) / tau) * res0[0][j] + h1 * PSI10(j * h2) - ((h1 * h1) / (2 * h2 * h2)) * (res0[0][j + 1] - 2 * res0[0][j] + res0[0][j - 1]) - ((h1 * h1) / 4) * (F(0, j * h2) + F(h1 / 2, j * h2));
        rightb1[n] = -((h1 * h1) / tau) * res0[n][j] + h1 * PSI11(j * h2) - ((h1 * h1) / (2 * h2 * h2)) * (res0[n][j + 1] - 2 * res0[n][j] + res0[n][j - 1]) - ((h1 * h1) / 4) * (F(n * h1, j * h2) + F(h1 * (n - 1 / 2), j * h2));

        Progonka(0, n + 1, L1, D1, U1, rightb1, res[j]);
    }
    for (int i = 0; i < n + 1; i++)
    {

        for (int j = 0; j < m + 1; j++)
        {
            resy[i][j] = res[j][i];
        }
    }
    for (int i = 0; i < n + 1; i++)
    {

        for (int j = 0; j < m + 1; j++)
        {
            res1[i][j] = resy[i][j];
        }
    }
    for (int i = 1; i < n; i++)
    {
        for (int j = 1; j < m; j++)
        {
            rightb2[j] = -(2 * (h2 * h2) / tau) * resy[i][j] - ((h2 * h2) * (resy[i + 1][j] - 2 * resy[i][j] + resy[i - 1][j])) / (h1 * h1) - (h2 * h2) * F(i * h1, j * h2);
        }
        rightb2[1] = rightb2[1] - resy[i][0];
        rightb2[m - 1] = rightb2[m - 1] - resy[i][m];
        Progonka(1, m, L2, D2, U2, rightb2, res1[i]);
    }
    for (int j = 0; j < m + 1; j++)
    {
        L2[j] = 1.;
        D2[j] = -2 * (1. + (h2 * h2) / tau);
        U2[j] = L2[j];
    }
    for (int j = 1; j < m; j++)
    {
        rightb2[j] = -(2*(h2 * h2) / tau) * resy[n][j] - (2 * h2 * h2 / h1) * PSI11(j * h2) - (2 * h2 * h2 / (h1*h1)) * (resy[n][j] - resy[n-1][j]) - ((h1 * h1) / 4) * (F(0, j * h2) + F(h1 / 2, j * h2));
    }
    rightb2[1] = rightb2[1] - resy[n][0];
    rightb2[m - 1] = rightb2[m - 1] - resy[n][m];
    Progonka(1, m, L2, D2, U2, rightb2, res1[n]);
    for (int j = 1; j < m; j++)
    {
        rightb2[j] = -(2*(h2 * h2) / tau) * resy[0][j] + (2* h2*h2/h1)* PSI10(j * h2)- (2 * h2 * h2 / (h1*h1))* (resy[1][j] - resy[0][j]) - ((h2 * h2) / 4) * (F(n * h1, j * h2) + F(h1 * (n - 1 / 2), j * h2));
    }
    rightb2[1] = rightb2[1] - resy[0][0];
    rightb2[m - 1] = rightb2[m - 1] - resy[0][m];
    Progonka(1, m, L2, D2, U2, rightb2, res1[0]);
    for (int i = 0; i < m + 1; i++)
    {
        delete[] res[i];
    }
    for (int i = 0; i < n + 1; i++)
    {
        delete[] resy[i];
    }
    delete[] L1;
    delete[] D1;
    delete[] U1;
    delete[] rightb1;
    delete[] L2;
    delete[] D2;
    delete[] U2;
    delete[] rightb2;
    delete[] res;
    delete[] resy;
}
void MStept3(double** res0, double** res1, int n, int m)
{
    double** res;
    double** resy;
    res = new double* [m + 1];
    for (int i = 0; i < m + 1; i++)
    {
        res[i] = new double[n + 1];
    }
    resy = new double* [n + 1];
    for (int i = 0; i < n + 1; i++)
    {
        resy[i] = new double[m + 1];
    }
    double* L1;
    L1 = new double[n + 1];
    double* D1;
    D1 = new double[n + 1];
    double* U1;
    U1 = new double[n + 1];
    double* rightb1;
    rightb1 = new double[n + 1];
    double* L2;
    L2 = new double[m + 1];
    double* D2;
    D2 = new double[m + 1];
    double* U2;
    U2 = new double[m + 1];
    double* rightb2;
    rightb2 = new double[m + 1];
    for (int j = 0; j < n + 1; j++)
    {
        L1[j] = 1;
        D1[j] = -2 * (1. + (h1 * h1) / tau);
        U1[j] = L1[j];
    }
    for (int j = 0; j < m + 1; j++)
    {
        L2[j] = 1.;
        D2[j] = -2 * (1. + (h2 * h2) / tau);
        U2[j] = L2[j];
    }
    L1[0] = 0;
    D1[0] = -((h1 * h1) / tau + 1. );
    U1[0] = 1.;
    L1[n]= 1.;
    D1[n]= -((h1 * h1) / tau + 1.);
    U1[n] = 0;
    for (int i = 0; i < m+1 ; i++)
    {
        res0[i][0] = (i * h2) * (i * h2);
        res0[i][m] = 1.0+ (i * h2) * (i * h2);
    }
    for (int i = 0; i < n + 1; i++)
    {
        for (int j = 0; j < m + 1; j++)
        {
            res[j][i] = res0[i][j];
            resy[i][j] = res0[i][j];
        }

    }
    
    for (int j = 1; j < m; j++)
    {
        for (int i = 1; i < n; i++)
        {
            rightb1[i] = -(2 * (h1 * h1) / tau) * res0[i][j] - ((h1 * h1) * (res0[i][j + 1] - 2 * res0[i][j] + res0[i][j - 1])) / (h2 * h2) - (h1 * h1) * F(i * h1, j * h2);
        }
        rightb1[0] = -((h1*h1)/tau)*res0[0][j]-h1*PSI10(j*h2)-((h1*h1)/(2*h2*h2))*(res0[0][j+1]-2*res0[0][j]+res0[0][j-1])-((h1*h1)/4)*(F(0,j*h2)+F(h1/2,j*h2));
        rightb1[n] = -((h1 * h1) / tau) * res0[n][j] + h1 * PSI11(j * h2) - ((h1 * h1) / (2 * h2 * h2)) * (res0[n][j + 1] - 2 * res0[n][j] + res0[n][j - 1]) - ((h1 * h1) / 4) * (F(n*h1, j * h2) + F(h1*(n-1 / 2), j * h2));
        
        Progonka(0, n+1, L1, D1, U1, rightb1, res[j]);
    }
    
    for (int i = 0; i < n + 1; i++)
    {
        
        for (int j = 0; j < m + 1; j++)
        {
            resy[i][j] = res[j][i];
        }
    }
    for (int i = 0; i < n + 1; i++)
    {

        for (int j = 0; j < m +1; j++)
        {
            res1[i][j] = resy[i][j];
        }
    }
    for (int i = 1; i < n; i++)
    {
        for (int j = 1; j < m; j++)
        {
            rightb2[j] = -(2 * (h2 * h2) / tau) * resy[i][j] - ((h2 * h2) * (resy[i + 1][j] - 2 * resy[i][j] + resy[i + 1][j])) / (h1 * h1) - (h2 * h2) * F(i * h1, j * h2);
        }
        rightb2[1] = rightb2[1] - resy[i][0];
        rightb2[m - 1] = rightb2[m - 1] - resy[i][m];
        Progonka(1, m, L2, D2, U2, rightb2, res1[i]);
    }
    for (int j = 1; j < m; j++)
    {
        rightb2[j] = -(2 * (h2 * h2) / tau) * resy[n][j] - (h2 * h2) * (PSI11(j * h2)-(resy[n][j]-resy[n-1][j])/h1) / (h1)-(h2 * h2) * F(n * h1, j * h2);
    }
    rightb2[1] = rightb2[1] - resy[n][0];
    rightb2[m - 1] = rightb2[m - 1] - resy[n][m];
    Progonka(1, m, L2, D2, U2, rightb2, res1[n]);
    for (int j = 1; j < m; j++)
    {
        rightb2[j] = -(2 * (h2 * h2) / tau) * resy[0][j] -(h2 * h2) * ((resy[1][j] - resy[0][j]) / h1 - PSI10(j * h2)) / (h1)-(h2 * h2) * F(0 * h1, j * h2);
    }
    rightb2[1] = rightb2[1] - resy[0][0];
    rightb2[m - 1] = rightb2[m - 1] - resy[0][m];
    Progonka(1, m, L2, D2, U2, rightb2, res1[0]);
    for (int i = 0; i < m + 1; i++)
    {
        delete[] res[i];
    }
    for (int i = 0; i < n + 1; i++)
    {
        delete[] resy[i];
    }
    delete[] L1;
    delete[] D1;
    delete[] U1;
    delete[] rightb1;
    delete[] L2;
    delete[] D2;
    delete[] U2;
    delete[] rightb2;
    delete[] res;
    delete[] resy;
}
void MStept2(double** res0,double** res1,int n,int m)
{
    double** res;
    double** resy;
    res = new double*[m+1];
    for (int i = 0; i < m+1; i++)
    {
        res[i] = new double[n+1];
    }
    resy = new double* [n + 1];
    for (int i = 0; i < n + 1; i++)
    {
        resy[i] = new double[m + 1];
    }
    double* L1;
    L1 = new double[n+1];
    double* D1;
    D1 = new double[n+1];
    double* U1;
    U1 = new double[n+1];
    double* rightb1;
    rightb1 = new double[n+1];
    double* L2;
    L2 = new double[m+1];
    double* D2;
    D2 = new double[m+1];
    double* U2;
    U2 = new double[m+1];
    double* rightb2;
    rightb2= new double[m+1];
    for (int j = 0; j < n+1; j++)
    {
        L1[j] = 1;
        D1[j] = -2 * (1. + (h1 * h1) / tau);
        U1[j] = L1[j];
    }
    for (int j = 0; j < m + 1; j++)
    {
        L2[j] = 1.;
        D2[j] = -2 * (1. + (h2 * h2) / tau);
        U2[j] = L2[j];
    }
    for (int i = 1; i < n; i++)
    {
        res0[i][0] = 1.0 + i * h2;
        res0[i][m] = 1.0 + i * h2;
    }
    for (int j = 0; j < m + 1; j++)
    {
        res0[0][j] = res0[2][j] - 2 * h2;
        res0[n][j] = res0[m-2][j] + 2 * h2;
    }
    for (int i=0; i < n+1; i++)
    {
        for (int j = 0; j < m+1; j++)
        {
            res[j][i] = res0[i][j];
            res1[i][j] = res0[i][j];
            
        }
        
    }
    
    for (int j = 1; j < m; j++)
    {
        for (int i = 1; i < n; i++)
        {
            rightb1[i] = -(2* (h1 * h1) / tau) * res0[i][j] - ((h1 * h1) * (res0[i][j + 1] - 2 * res0[i][j] + res0[i][j - 1]))/ (h2 * h2) - (h1 * h1) *F(i * h1, j * h2);
        }
        rightb1[1] = rightb1[1] -   res0[0][j];
        rightb1[n-1] = rightb1[n-1] -  res0[n][j];
        Progonka(1, n, L1, D1, U1, rightb1, res[j]);
    }
    for (int i = 0; i < n + 1; i++)
    {
        for (int j = 0; j < m + 1; j++)
        {
            resy[i][j] = res[j][i];
        }
    }
    for (int i = 1; i < n; i++)
    {
        for (int j = 1; j < m; j++)
        {
            rightb2[j] = -(2* (h2 * h2) / tau) * resy[i][j] -((h2 * h2)* (resy[i+1][j] - 2 * resy[i][j] + resy[i+1][j])) / (h1 * h1) - (h2 * h2)*F(i * h1, j * h2);
        }
        rightb2[1] = rightb2[1] -  resy[i][0];
        rightb2[m - 1] = rightb2[m - 1] - resy[i][m];
        Progonka(1, m, L2, D2, U2, rightb2, res1[i]);
    }
    for (int i = 0; i < m + 1; i++)
    {
        delete[] res[i];
    }
    delete[] L1;
    delete[] D1;
    delete[] U1;
    delete[] rightb1;
    delete[] L2;
    delete[] D2;
    delete[] U2;
    delete[] rightb2;
    delete[] res;
}
void MStept1(double** res0, double** res1, int n, int m)
{
    double** res;
    double** resy;
    res = new double* [m + 1];
    for (int i = 0; i < m + 1; i++)
    {
        res[i] = new double[n + 1];
    }
    resy = new double* [n + 1];
    for (int i = 0; i < n + 1; i++)
    {
        resy[i] = new double[m + 1];
    }
    double* L1;
    L1 = new double[n + 1];
    double* D1;
    D1 = new double[n + 1];
    double* U1;
    U1 = new double[n + 1];
    double* rightb1;
    rightb1 = new double[n + 1];
    double* L2;
    L2 = new double[m + 1];
    double* D2;
    D2 = new double[m + 1];
    double* U2;
    U2 = new double[m + 1];
    double* rightb2;
    rightb2 = new double[m + 1];
    for (int j = 0; j < n + 1; j++)
    {
        L1[j] = 1;
        D1[j] = -2 * (1. + (h1 * h1) / tau);
        U1[j] = L1[j];
    }
    for (int j = 0; j < m + 1; j++)
    {
        L2[j] = 1.;
        D2[j] = -2 * (1. + (h2 * h2) / tau);
        U2[j] = L2[j];
    }
    for (int i = 0; i < n + 1; i++)
    {
        res0[i][0] = 1.0;
        res0[i][m] = 1.0;
    }
    for (int j = 0; j < m + 1; j++)
    {
        res0[0][j] = 1.0;
        res0[n][j] = 1.0;
    }
    for (int i = 0; i < n + 1; i++)
    {
        for (int j = 0; j < m + 1; j++)
        {
            res[j][i] = res0[i][j];
            res1[i][j] = res0[i][j];

        }

    }

    for (int j = 1; j < m; j++)
    {
        for (int i = 1; i < n; i++)
        {
            rightb1[i] = -(2 * (h1 * h1) / tau) * res0[i][j] - ((h1 * h1) * (res0[i][j + 1] - 2 * res0[i][j] + res0[i][j - 1])) / (h2 * h2) - (h1 * h1) * F(i * h1, j * h2);
        }
        rightb1[1] = rightb1[1] - res0[0][j];
        rightb1[n - 1] = rightb1[n - 1] - res0[n][j];
        Progonka(1, n, L1, D1, U1, rightb1, res[j]);
    }
    for (int i = 0; i < n + 1; i++)
    {
        for (int j = 0; j < m + 1; j++)
        {
            resy[i][j] = res[j][i];
        }
    }
    for (int i = 1; i < n; i++)
    {
        for (int j = 1; j < m; j++)
        {
            rightb2[j] = -(2 * (h2 * h2) / tau) * resy[i][j] - ((h2 * h2) * (resy[i + 1][j] - 2 * resy[i][j] + resy[i + 1][j])) / (h1 * h1) - (h2 * h2) * F(i * h1, j * h2);
        }
        rightb2[1] = rightb2[1] - resy[i][0];
        rightb2[m - 1] = rightb2[m - 1] - resy[i][m];
        Progonka(1, m, L2, D2, U2, rightb2, res1[i]);
    }
    for (int i = 0; i < m + 1; i++)
    {
        delete[] res[i];
    }
    delete[] L1;
    delete[] D1;
    delete[] U1;
    delete[] rightb1;
    delete[] L2;
    delete[] D2;
    delete[] U2;
    delete[] rightb2;
    delete[] res;
}
void test1()
{
    double** res0;
    double** res1;
    double L1, L2;
    L1 = 1;
    L2 = 1;
    h1 = 0.01;
    h2 = 0.01;
    tau = 0.0001;
    indP = 0;
    n = int(L1 / h1);
    m = int(L2 / h2);
    res0 = new double*[n+1];
    res1 = new double*[n+1];
    for (int i = 0; i < n+1; i++)
    {
        res0[i] = new double[m+1];
        res1[i] = new double[m+1];
    }
    
    for (int i = 0; i < n + 1; i++)
    {
        for (int j = 0; j < m + 1; j++)
        {
            res1[i][j] = U0(i * h1, i * h2)+0.1;
        }
    }
    
    while (NormV2(res0, res1, n, m) > eps/10)
    {
        for(int i = 0; i < n+1; i++)
        {
            for (int j = 0; j < m+1; j++)
            {
                res0[i][j] = res1[i][j];
            }
        }
        MStept1(res0, res1, n, m);
    }
    std::ofstream ans;
    ans.open("test1.txt");
    ans << h1 << ' ' << h2 << std::endl;
    for (int i = 0; i < n+1; i++)
    {
        for (int j = 0; j < m+1; j++)
        {
            ans << res1[i][j] << ' ';
        }
        ans << std::endl;
    }
    delete[] res0;
    delete[] res1;
}
void test2()
{
    double** res0;
    double** res1;
    double L1, L2;
    L1 = 1;
    L2 = 1;
    h1 = 0.1;
    h2 = 0.1;
    tau = 0.001;
    indP = 1;
    n = int(L1 / h1);
    m = int(L2 / h2);
    res0 = new double* [n + 1];
    res1 = new double* [n + 1];
    for (int i = 0; i < n + 1; i++)
    {
        res0[i] = new double[m + 1];
        res1[i] = new double[m + 1];
    }

    for (int i = 0; i < n + 1; i++)
    {
        for (int j = 0; j < m + 1; j++)
        {
            res1[i][j] = U0(i * h1, i * h2) + 0.1;
        }
    }

    while (NormV2(res0, res1, n, m) > eps / 10)
    {
        for (int i = 0; i < n + 1; i++)
        {
            for (int j = 0; j < m + 1; j++)
            {
                res0[i][j] = res1[i][j];
            }
        }
        MStept(res0, res1, n, m);
    }
    std::ofstream ans;
    ans.open("test2.txt");
    ans << h1 << ' ' << h2 << std::endl;
    for (int i = 0; i < n + 1; i++)
    {
        for (int j = 0; j < m + 1; j++)
        {
            ans << res1[i][j] << ' ';
        }
        ans << std::endl;
    }
    delete[] res0;
    delete[] res1;
}
void test3()
{
    double** res0;
    double** res1;
    double L1, L2;
    L1 = 1;
    L2 = 1;
    h1 = 0.05;
    h2 = 0.05;
    tau = 0.0001;
    indP = 2;
    n = int(L1 / h1);
    m = int(L2 / h2);
    res0 = new double* [n + 1];
    res1 = new double* [n + 1];
    for (int i = 0; i < n + 1; i++)
    {
        res0[i] = new double[m + 1];
        res1[i] = new double[m + 1];
    }
    
    for (int i = 0; i < n + 1; i++)
    {
        for (int j = 0; j < m + 1; j++)
        {
            res1[i][j] = U0(i * h1, i * h2);
        }
    }
    while (NormV2(res0, res1, n, m) > eps / 10)
    {
        for (int i = 0; i < n + 1; i++)
        {
            for (int j = 0; j < m + 1; j++)
            {
                res0[i][j] = res1[i][j];
            }
        }
        MStept(res0, res1, n, m);
    }
    std::ofstream ans;
    ans.open("test3.txt");
    ans << h1 << ' ' << h2 << std::endl;
    for (int i = 0; i < n + 1; i++)
    {
        for (int j = 0; j < m + 1; j++)
        {
            ans << res1[i][j] << ' ';
        }
        ans << std::endl;
    }
    delete[] res0;
    delete[] res1;
}
void testv16()
{
    double** res0;
    double** res1;
    double L1, L2;
    L1 = 1;
    L2 = 2;
    h1 = 0.05;
    h2 = 0.05;
    tau = 0.0001;
    indP = 3;
    n = int(L1 / h1);
    m = int(L2 / h2);
    res0 = new double* [n + 1];
    res1 = new double* [n + 1];
    for (int i = 0; i < n + 1; i++)
    {
        res0[i] = new double[m + 1];
        res1[i] = new double[m + 1];
    }

    for (int i = 0; i < n + 1; i++)
    {
        for (int j = 0; j < m + 1; j++)
        {
            res1[i][j] = U0(i * h1, i * h2) + 0.1;
        }
    }

    while (NormV2(res0, res1, n, m) > eps / 10)
    {
        for (int i = 0; i < n + 1; i++)
        {
            for (int j = 0; j < m + 1; j++)
            {
                res0[i][j] = res1[i][j];
            }
        }

        MStept(res0, res1, n, m);

    }
    std::ofstream ans;
    ans.open("testv16.txt");
    ans << h1 << ' ' << h2 << std::endl;
    for (int i = 0; i < n + 1; i++)
    {
        for (int j = 0; j < m + 1; j++)
        {
            ans << res1[i][j] << ' ';
        }
        ans << std::endl;
    }
    delete[] res0;
    delete[] res1;
}
void testv5()
{
    double** res0;
    double** res1;
    double L1, L2;
    L1 =pi;
    L2 = pi;
    h1 = 0.05;
    h2 = 0.05;
    tau = 0.0001;
    indP = 4;
    n = int(L1 / h1);
    m = int(L2 / h2);
    res0 = new double* [n + 1];
    res1 = new double* [n + 1];
    for (int i = 0; i < n + 1; i++)
    {
        res0[i] = new double[m + 1];
        res1[i] = new double[m + 1];
    }

    for (int i = 0; i < n + 1; i++)
    {
        for (int j = 0; j < m + 1; j++)
        {
            res1[i][j] = U0(i * h1, i * h2) + 0.1;
        }
    }

    while (NormV2(res0, res1, n, m) > eps / 10)
    {
        for (int i = 0; i < n + 1; i++)
        {
            for (int j = 0; j < m + 1; j++)
            {
                res0[i][j] = res1[i][j];
            }
        }

        MStept(res0, res1, n, m);

    }
    std::ofstream ans;
    ans.open("testv5.txt");
    ans << h1 << ' ' << h2 << std::endl;
    for (int i = 0; i < n + 1; i++)
    {
        for (int j = 0; j < m + 1; j++)
        {
            ans << res1[i][j] << ' ';
        }
        ans << std::endl;
    }
    delete[] res0;
    delete[] res1;
}
int main()
{
    eps = 0.001;
    //test1();
    //test2();
    //test3();
    testv16();
    //testv5();
    
    std::cout << "Hello World!\n";
}


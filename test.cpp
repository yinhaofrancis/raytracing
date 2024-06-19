#include <iostream>
#include <math.h>
#include <random>

double random_double(){
    static std::uniform_real_distribution<double> distribution(0, 1);
    static std::mt19937 generator;
    return distribution(generator);
}
double f(double d) {
    // return 8.0 * pow(d, 1.0/3.0);
    return 2.0 * d;
}

double pdf(double x) {
    // return (3.0/8.0) * x*x;
    return 0.5;
}

int main() {
    int N = 100000;
    auto sum = 0.0;
    for (int i = 0; i < N; i++) {
        auto x = f(random_double());
        sum += x*x / pdf(x);
    }
    std::cout << "I = " << sum / N << '\n';
}
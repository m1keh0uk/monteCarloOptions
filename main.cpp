#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>

// Function to generate a standard normal random variable
double generateNormalRandom() {
    double u1 = static_cast<double>(std::rand()) / RAND_MAX;
    double u2 = static_cast<double>(std::rand()) / RAND_MAX;
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
}

double sumPathResults(const std::vector<double>& path, double trial_total, double S){
    return trial_total += (path.back() - S) / S;
}

std::vector<double> simGBM(double S0, double mu, double sigma, double T, int numSteps) {
    std::vector<double> path(numSteps + 1);
    double dt = T / numSteps; // Time each step

    path[0] = S0; // Initial  price

    for (int i = 1; i <= numSteps; ++i) {
        double Z = generateNormalRandom();
        path[i] = path[i - 1] * std::exp((mu - 0.5 * sigma * sigma) * dt + sigma * std::sqrt(dt) * Z);
    }

    return path;
}

// Function to print the GBM path
void printPath(const std::vector<double>& path, int i) {
        std::cout << "Trial " << i << ": " << path.back() << std::endl;
}

int main() {
    std::srand(static_cast<unsigned int>(std::time(0))); // Seed the random number generator

    double S0 = 100.0;   // price
    double mu = 0.05;    // mu
    double sigma = 0.2;  // Vol
    double T = 1.0;      // Time
    int numSteps = 365;  // Num steps

    int trials = 100; // num paths created
    double trial_total = 0;

    // Simulate the GBM path
    for (int i = 1; i <= trials; ++i) {
        std::vector<double> gbmPath = simGBM(S0, mu, sigma, T, numSteps);
        printPath(gbmPath, i); // Print the GBM path
        trial_total = sumPathResults(gbmPath, trial_total, S0);
    }

    double avgReturn = trial_total / trials;
    std::cout << "Average Anuallized Percentage Return: " << avgReturn <<"%" << std::endl;
    

    return 0;
}

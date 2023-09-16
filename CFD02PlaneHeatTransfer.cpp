// CFD02PlaneHeatTransfer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <HeatEquationSolver.hpp>
#include <chrono>
#include <thread>

int main()
{
    Range2D domain;
    Index2D counts;
    std::cout << "Enter lower value of x: ";
    std::cin >> domain.xrange.first;
    std::cout << "\n";
    std::cout << "Enter upper value of x: ";
    std::cin >> domain.xrange.second;
    std::cout << "\n";
    std::cout << "Enter lower value of y: ";
    std::cin >> domain.yrange.first;
    std::cout << "\n";
    std::cout << "Enter upper value of y: ";
    std::cin >> domain.yrange.second;
    std::cout << "\n";
    
    std::cout << "Enter number of nodes in x direction: ";
    std::cin >> counts.xcount;
    std::cout << "\n";
    std::cout << "Enter number of nodes in y direction: ";
    std::cin >> counts.ycount;
    std::cout << "\n";
    ArshiaCFD::HeatEquationSolver P(domain, counts, 1.9);
    auto t1 = std::chrono::high_resolution_clock::now();
    P.run();
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> ms_double = t2 - t1;
    std::cout << "Total time = " << ms_double.count() << "ms\n";
    //Index2D counts1;
    //domain.xrange.first = 0;
    //domain.xrange.second = 1;
    //domain.yrange.first = 0;
    //domain.yrange.second = 1;
    //counts1.xcount = 100;
    //counts1.ycount = 100;
    
    //ArshiaCFD::HeatEquationSolver P1(domain, counts1, 1.65);
    //ArshiaCFD::HeatEquationSolver P2(domain, counts1, 1.87);
    //ArshiaCFD::HeatEquationSolver P3(domain, counts1, 1);
    //ArshiaCFD::HeatEquationSolver P4(domain, counts1, 1.70);
    //ArshiaCFD::HeatEquationSolver P5(domain, counts1, 1.80);
    //ArshiaCFD::HeatEquationSolver P6(domain, counts1, 1.85);
    //ArshiaCFD::HeatEquationSolver P7(domain, counts1, 1.90);
    //ArshiaCFD::HeatEquationSolver P8(domain, counts1, 1.95);
    //
    //auto t1 = std::chrono::high_resolution_clock::now();
    //std::thread Work1(&ArshiaCFD::HeatEquationSolver::run, P1);
    //std::thread Work2(&ArshiaCFD::HeatEquationSolver::run, P2);
    //std::thread Work3(&ArshiaCFD::HeatEquationSolver::run, P3);
    //std::thread Work4(&ArshiaCFD::HeatEquationSolver::run, P4);
    //std::thread Work5(&ArshiaCFD::HeatEquationSolver::run, P5);
    //std::thread Work6(&ArshiaCFD::HeatEquationSolver::run, P6);
    //std::thread Work7(&ArshiaCFD::HeatEquationSolver::run, P7);
    //std::thread Work8(&ArshiaCFD::HeatEquationSolver::run, P8);
    //Work1.join();
    //Work2.join();
    //Work3.join();
    //Work4.join();
    //Work5.join();
    //Work6.join();
    //Work7.join();
    //Work8.join();
    //auto t2 = std::chrono::high_resolution_clock::now();
    //std::chrono::duration<double, std::milli> ms_double = t2 - t1;
    //std::cout << "Total time = " << ms_double.count() << "ms\n";
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file

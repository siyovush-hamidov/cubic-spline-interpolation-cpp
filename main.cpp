#include "cubicSpline.h"
#include <iostream>
#include <fstream>
#include <cmath>

float f(float x) {
    return log(x) + exp(1) / x;
}

int main()
{
    float x1[] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2 };
    float y1[20];

    for (int i = 0; i < 20; i++) {
        y1[i] = f(x1[i]);
    }

    float y2[39];
    float c0[19], c1[19], c2[19], c3[19];

    cubicSpline(c0, c1, c2, c3, x1, y1, 20);

    float x2[39];
    for (int i = 0; i < 39; i++) {
        x2[i] = 0.1 + i * 0.05;
    }

    evaluate_polynomial(y2, x2, 39, x1, 20, c0, c1, c2, c3);

    std::ofstream input_file("input.txt");
    
    for (int i = 0; i < 20; i++) {
        input_file << x1[i] << " " << y1[i] << "\n";
    }

    input_file.close();

    std::ofstream output_file("output.txt");

    for (int i = 0; i < 39; i++) {
        output_file << x2[i] << " " << y2[i] << "\n";
    }
    output_file.close();

    std::cout << "Interpolation results have been written to input.txt and output.txt\n";

    return 0;
}

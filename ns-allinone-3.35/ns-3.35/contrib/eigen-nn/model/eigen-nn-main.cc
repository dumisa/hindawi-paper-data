/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */

#include "eigen-nn.h"

namespace ns3 {

// main.cpp
 
// don't forget to include out neural network
#include "NeuralNetwork.hpp"
 
void ReadCSV(std::string filename, std::vector<RowVector*>& data)
{
    data.clear();
    std::ifstream file(filename);
    std::string line, word;
    // determine number of columns in file
    getline(file, line, '\n');
    std::stringstream ss(line);
    std::vector<Scalar> parsed_vec;
    while (getline(ss, word, ', ')) {
        parsed_vec.push_back(Scalar(std::stof(&word[0])));
    }
    uint cols = parsed_vec.size();
    data.push_back(new RowVector(cols));
    for (uint i = 0; i < cols; i++) {
        data.back()->coeffRef(1, i) = parsed_vec[i];
    }
 
    // read the file
    if (file.is_open()) {
        while (getline(file, line, '\n')) {
            std::stringstream ss(line);
            data.push_back(new RowVector(1, cols));
            uint i = 0;
            while (getline(ss, word, ', ')) {
                data.back()->coeffRef(i) = Scalar(std::stof(&word[0]));
                i++;
            }
        }
    }
}

void genData(std::string filename)
{
    std::ofstream file1(filename + "-in");
    std::ofstream file2(filename + "-out");
    for (uint r = 0; r < 1000; r++) {
        Scalar x = rand() / Scalar(RAND_MAX);
        Scalar y = rand() / Scalar(RAND_MAX);
        file1 << x << ", " << y << std::endl;
        file2 << 2 * x + 10 + y << std::endl;
    }
    file1.close();
    file2.close();
}
 
typedef std::vector<RowVector*> data;
int main()
{
    NeuralNetwork n({ 2, 3, 1 });
    data in_dat, out_dat;
    genData("test");
    ReadCSV("test-in", in_dat);
    ReadCSV("test-out", out_dat);
    n.train(in_dat, out_dat);
    return 0;
}


}


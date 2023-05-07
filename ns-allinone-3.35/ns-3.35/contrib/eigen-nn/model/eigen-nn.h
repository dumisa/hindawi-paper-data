/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
#ifndef EIGEN_NN_H
#define EIGEN_NN_H

// eigen-nn.h
#include <eigen3/Eigen/Eigen>
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>

// use typedefs for future ease for changing data types like : float to double
typedef float Scalar;
typedef Eigen::MatrixXf Matrix;
typedef Eigen::RowVectorXf RowVector;
typedef Eigen::VectorXf ColVector;


namespace ns3 {

 
// neural network implementation class!
class NeuralNetwork 
{
public:
    // constructor
    NeuralNetwork(std::vector<uint> topology, Scalar learningRate = Scalar(0.005));
 
    // function for forward propagation of data
    void propagateForward(RowVector& input);
 
    // function for backward propagation of errors made by neurons
    void propagateBackward(RowVector& output);
 
    // function to calculate errors made by neurons in each layer
    void calcErrors(RowVector& output);
 
    // function to update the weights of connections
    void updateWeights();
 
    // function to train the neural network give an array of data points
    void train(std::vector<RowVector*> input_data, std::vector<RowVector*> output_data);
 
    // storage objects for working of neural network
    /*
          use pointers when using std::vector<Class> as std::vector<Class> calls destructor of
          Class as soon as it is pushed back! when we use pointers it can't do that, besides
          it also makes our neural network class less heavy!! It would be nice if you can use
          smart pointers instead of usual ones like this
        */
    std::vector<RowVector*> neuronLayers; // stores the different layers of out network
    std::vector<RowVector*> cacheLayers; // stores the unactivated (activation fn not yet applied) values of layers
    std::vector<RowVector*> deltas; // stores the error contribution of each neurons
    std::vector<Matrix*> weights; // the connection weights itself
    Scalar learningRate;

    std::vector<uint> topology;
};




}

#endif /* EIGEN_NN_H */


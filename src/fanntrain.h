/*
 *
 *  Fast Artificial Neural Network (fann) C++ Wrapper Sample
 *
 *  C++ wrapper XOR sample with functionality similar to xor_train.c
 *
 *  Copyright (C) 2004-2006 created by freegoldbar (at) yahoo dot com
 *
 *  This wrapper is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This wrapper is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include "floatfann.h"
#include "fann_cpp.h"

#include <ios>
#include <iostream>
#include <iomanip>
using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::left;
using std::right;
using std::showpos;
using std::noshowpos;

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <fstream>
#include <ctime>

void write_net(FANN::neural_net &net, std::string filename, std::string headerline = "");
void read_net(FANN::neural_net &net, std::string filename, std::string &headerline);
void read_net(FANN::neural_net &net, std::string filename);
void simpler_train(FANN::neural_net &net, FANN::training_data &data,
		  std::string network_file, std::string parameter_file,
                  std::string continue_train_file = "", std::string custom_header = "",
		  int iterations_between_save = 10);
void simpler_test(FANN::neural_net &net, std::vector<std::vector<float> > &input,
                 std::vector<std::vector<float> > &output, std::string parameter_file = "",
                 unsigned int num_output = 2);

// Callback function that simply prints the information to cout
int print_callback(FANN::neural_net &net, FANN::training_data &train,
    unsigned int max_epochs, unsigned int epochs_between_reports,
    float desired_error, unsigned int epochs, void *user_data)
{
      cout << "EPOCHS     " << setw(8) << epochs << ". "
         << "CURRENT ERROR: " << left << net.get_MSE() << right << endl;
    return 0;
}

// Test function that demonstrates usage of the fann C++ wrapper
void simple_train(FANN::neural_net &net, FANN::training_data &data,
		  std::string network_file, std::string parameter_file)
{
    // Load parameters 
    std::ifstream parameter(parameter_file);
    float learning_rate = 0.7f;
    unsigned int num_layers = 3;
    unsigned int num_input = 4;
    unsigned int num_hidden[3] = {3, 0, 0};
    unsigned int num_output = 1;
    float desired_error = 0.001f;
    unsigned int max_iterations = 30000;
    unsigned int iterations_between_reports = 1000;
    std::string line_description; 
    parameter >> line_description >> learning_rate;
    parameter >> line_description >> num_layers;
    parameter >> line_description >> num_input;
    for(int i = 0; i < num_layers - 2; ++i)
      parameter >> line_description >> num_hidden[i];
    parameter >> line_description >> num_output;
    parameter >> line_description >> desired_error;
    parameter >> line_description >> max_iterations;
    parameter >> line_description >> iterations_between_reports;

    // Create network
    cout << endl << "CREATING NETWORK" << endl;
    switch(num_layers) {
    case 3:
      net.create_standard(num_layers, num_input, num_hidden[0], num_output);
      break;
    case 4:
      net.create_standard(num_layers, num_input, num_hidden[0], num_hidden[1], num_output);
      break;
    case 5:
      net.create_standard(num_layers, num_input, num_hidden[0], num_hidden[1], num_hidden[2], num_output);
      break;
    default:
      break;
    }
    net.set_learning_rate(learning_rate);
    net.set_activation_steepness_hidden(1.0);
    net.set_activation_steepness_output(1.0);
    net.set_activation_function_hidden(FANN::SIGMOID_SYMMETRIC_STEPWISE);
    net.set_activation_function_output(FANN::SIGMOID_SYMMETRIC_STEPWISE);
    // Set additional properties such as the training algorithm
    //net.set_training_algorithm(FANN::TRAIN_QUICKPROP);
    
    // Output network type and parameters
    cout << endl << "Network Type                         :  ";
    switch (net.get_network_type()) {
	    case FANN::LAYER:
		cout << "LAYER" << endl;
		break;
	    case FANN::SHORTCUT:
		cout << "SHORTCUT" << endl;
		break;
	    default:
		cout << "UNKNOWN" << endl;
		break;
    }
    net.print_parameters();
    
    // Training network
    cout << endl << "TRAINING NETWORK" << endl;

    // Initialize and train the network with the data
    net.init_weights(data);

    cout << "MAX EPOCHS " << setw(8) << max_iterations << ". "
            << "DESIRED ERROR: " << left << desired_error << right << endl;
    net.set_callback(print_callback, NULL);
    net.train_on_data(data, max_iterations,
        iterations_between_reports, desired_error);
	  cout << endl << "SAVING NETWORK." << endl;
        // Save the network in floating point and fixed point
    net.save(network_file);
}

void write_net(FANN::neural_net &net, std::string filename, std::string headerline) {
  std::string tmpnetfile = ".";
  srand(time(NULL));
  for(int i = 0; i < 12; ++i) tmpnetfile += (char)('a'+rand()%26);
  net.save(tmpnetfile);
  std::ofstream out(filename, std::ios::trunc);
  if(headerline != "") out << headerline << endl;
  std::ifstream tmpnet(tmpnetfile);
  out << tmpnet.rdbuf();
  out.close();
  tmpnet.close();
  remove(tmpnetfile.c_str());
}

void read_net(FANN::neural_net &net, std::string filename, std::string &headerline) {
  std::string tmpnetfile = ".";
  srand(time(NULL));
  for(int i = 0; i < 12; ++i) tmpnetfile += (char)('a'+rand()%26);
  std::ifstream network_file(filename);
  std::ofstream tmpnet(tmpnetfile, std::ios::trunc);
  std::string buf;
  headerline = "";
  getline(network_file, buf);
  if(buf.substr(0, 4) != "FANN") headerline = buf;
  else tmpnet << buf << endl;
  while(getline(network_file, buf)) tmpnet << buf << endl;
  tmpnet.close();
  network_file.close();
  net.create_from_file(tmpnetfile);
  remove(tmpnetfile.c_str());
}

void read_net(FANN::neural_net &net, std::string filename) {
  std::string tempheader;
  read_net(net, filename, tempheader);
}

void simpler_train(FANN::neural_net &net, FANN::training_data &data, std::string network_file,
                  std::string parameter_file,
                  std::string continue_train_file, std::string custom_header, int iterations_between_save)
{
  // Load parameters 
  std::ifstream parameter(parameter_file);
  float learning_rate = 0.7f;
  unsigned int num_layers = 3;
  unsigned int num_input = 4;
  unsigned int num_hidden[3] = {3, 0, 0};
  unsigned int num_output = 1;
  float desired_error = 0.001f;
  unsigned int max_iterations = 30000;
  unsigned int iterations_between_reports = 1000;
  std::string line_description; 
  parameter >> line_description >> learning_rate;
  parameter >> line_description >> num_layers;
  parameter >> line_description >> num_input;
  for(int i = 0; i < num_layers - 2; ++i)
    parameter >> line_description >> num_hidden[i];
  parameter >> line_description >> num_output;
  parameter >> line_description >> desired_error;
  parameter >> line_description >> max_iterations;
  parameter >> line_description >> iterations_between_reports;
  
  // Create network
  if(continue_train_file == "") {
    cout << endl << "CREATING NETWORK" << endl;
    switch(num_layers) {
    case 3:
      net.create_standard(num_layers, num_input, num_hidden[0], num_output);
      break;
    case 4:
      net.create_standard(num_layers, num_input, num_hidden[0], num_hidden[1], num_output);
      break;
    case 5:
      net.create_standard(num_layers, num_input, num_hidden[0], num_hidden[1], num_hidden[2], num_output);
      break;
    default:
      break;
    }
    net.set_learning_rate(learning_rate);
    net.set_activation_steepness_hidden(1.0);
    net.set_activation_steepness_output(1.0);
    net.set_activation_function_hidden(FANN::SIGMOID_SYMMETRIC_STEPWISE);
    net.set_activation_function_output(FANN::SIGMOID_SYMMETRIC_STEPWISE);
  // Set additional properties such as the training algorithm
    net.set_training_algorithm(FANN::TRAIN_RPROP);
  } else {
    cout << endl << "CONTINUING FROM NETWORK FILE" << endl;
    read_net(net, continue_train_file, custom_header);
    net.set_training_algorithm(FANN::TRAIN_BATCH);
  }
  
  // Output network type and parameters
  cout << endl << "Network Type                         :  ";
  switch (net.get_network_type()) {
  case FANN::LAYER:
    cout << "LAYER" << endl;
    break;
  case FANN::SHORTCUT:
    cout << "SHORTCUT" << endl;
    break;
  default:
    cout << "UNKNOWN" << endl;
  break;
  }
  net.print_parameters();
  
  // Training network
  cout << endl << "TRAINING NETWORK" << endl;
  
  // Initialize and train the network with the data
  if(continue_train_file == "") net.init_weights(data);
  
  cout << "MAX EPOCHS " << setw(8) << max_iterations << ". "
       << "DESIRED ERROR: " << left << desired_error << right << endl;
  net.set_callback(print_callback, NULL);
  if(iterations_between_save < iterations_between_reports)
    iterations_between_save = iterations_between_reports;
  float min_MSE = 1;
  for(int current_epoch = 1; current_epoch <= max_iterations; ++current_epoch) {
    float cur_MSE = net.train_epoch(data);
    if(current_epoch == 1 || current_epoch % iterations_between_reports == 0)
      cout << "EPOCHS     " << setw(8) << current_epoch << ". "
         << "CURRENT ERROR: " << left << net.get_MSE() << right << endl;
    if(current_epoch % iterations_between_save == 0 && cur_MSE < min_MSE) {
      write_net(net, network_file, custom_header);
      min_MSE = cur_MSE;
    }
    if(cur_MSE <= desired_error) break;
  }
}


void simpler_test(FANN::neural_net &net, std::vector<std::vector<float> > &input,
                 std::vector<std::vector<float> > &output, std::string parameter_file,
                 unsigned int num_output)
{
  // Load parameters 
  if(parameter_file != "") {
    std::ifstream parameter(parameter_file);
    float learning_rate = 0.7f;
    unsigned int num_layers = 3;
    unsigned int num_input = 4;
    unsigned int num_hidden[3] = {3, 0, 0};
    float desired_error = 0.001f;
    unsigned int max_iterations = 30000;
    unsigned int iterations_between_reports = 1000;
    std::string line_description; 
    parameter >> line_description >> learning_rate;
    parameter >> line_description >> num_layers;
    parameter >> line_description >> num_input;
    for(int i = 0; i < num_layers - 2; ++i)
      parameter >> line_description >> num_hidden[i];
    parameter >> line_description >> num_output;
    parameter >> line_description >> desired_error;
    parameter >> line_description >> max_iterations;
    parameter >> line_description >> iterations_between_reports;
    parameter.close();
  }
  output.resize(input.size());
  // Run the network on the test data
  for (int i = 0; i < input.size(); ++i)
  {
    if(i%1000==0) {
      cout<<"NETWORK ON "<<i/1000<<"K / "<<input.size()/1000<<"K\r";
      cout.flush();
    }
    fann_type *calc_out = net.run(&(input[i][0]));
    for(int j = 0; j < num_output; ++j) output[i].push_back(*(calc_out + j)) ;
  }
}



/******************************************************************************/

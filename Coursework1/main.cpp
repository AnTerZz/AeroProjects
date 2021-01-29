#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include "Vector.h"

using namespace std;

Vector f(Vector y, double params[]); //Prototype for the main function


int main() {

	double parameters[8] = { 0 }; //Initialise a vector which will hold all the txt-read parameters

	ifstream cParams("parameters.txt");
	if (!cParams) { //check if the parameters file exists

		cout << "The parameters file was not found in the directory and has been created with default values" << endl;
		ofstream wParams;
		wParams.open("parameters.txt"); //creates a blank parameters file and writes default values to is as follows
		wParams << "1.5 "; //m1
		wParams << "2.0 "; //m2
		wParams << "1.0 "; //l1
		wParams << "1.0 "; //l2
		wParams << "1.571 "; //theta1_0
		wParams << "2.356 "; //theta2_0
		wParams << "0.01 "; //timestep h
		wParams << "10.0 "; //total duration t
		wParams.close();
	}

	ifstream rParams("parameters.txt"); //reads the parameters file
	if (rParams.is_open()) {
		string value[] = { "mass 1","mass 2","length 1","length 2","theta 1 start", "theta 2 start","timestep","duration" };
		for (int i = 0; i < 8; i++) {
			rParams >> parameters[i]; //store values from the read file into an array of 8 doubles 
			cout << "Read " << value[i] << " = " << parameters[i] << " from parameters.txt" << endl;
		}
		rParams.close();
	}

	else {
		cout << "Unable to open file";
	}

	ofstream resultsFile("output.txt"); //overwrites results in output.txt file with new solution
	resultsFile << "t\tx1\ty1\tx2\ty2\n"; //file header
	resultsFile << "---------------------------------------\n";


	Vector y(4);
	y[0] = parameters[4];
	y[1] = 0;
	y[2] = parameters[5];
	y[3] = 0;

	double h = parameters[6];
	double N = parameters[7] / h;

	for (int i = 0; i <= N; i++) {
		Vector k1 = f(y, parameters) * h; //As per operator overload, right side scalar multiplication
		Vector k2 = f(y + k1 / 2, parameters) * h;
		Vector k3 = f(y + k2 / 2, parameters) * h;
		Vector k4 = f(y + k3, parameters) * h;

		y = y + ((k1 + k4) / 6 + (k2 + k3) / 3);

		double l1 = parameters[2];
		double l2 = parameters[3];

		double x1, y1, x2, y2; //declaration and calculation of the position variables using the result from RK4
		x1 = l1 * sin(y[0]);
		y1 = -l1 * cos(y[0]);
		x2 = x1 + l2 * sin(y[2]);
		y2 = y1 - l2 * cos(y[2]);

		//Writes results to the output file, rounded to 4 decimal places
		resultsFile << i << "\t" << round(x1 * 10000) / 10000 << "\t" << round(y1 * 10000) / 10000 << "\t" << round(x2 * 10000) / 10000 << "\t" << round(y2 * 10000) / 10000 << endl;
	}
	resultsFile.close(); //closes output file once the loop is finished and no more data needs to be recorded
	cout << "The solution has been computed and can be found in the output.txt file" << endl; //Success message to user
	return 0;
}


//Main function which englobes the 4 equations in the first order system and calculates them for a vector of size 4 using the parameters array passed to it and previous step conditions
Vector f(Vector y, double params[]) {
	const double g = 9.81;
	double theta1 = y[0];
	double theta11 = y[1];
	double theta2 = y[2];
	double theta22 = y[3];

	double omp1 = ((-g * (2.0 * params[0] + params[1]) * sin(theta1) - params[1] * g * sin(theta1 - 2.0 * theta2) - 2.0 * sin(theta1 - theta2) * params[1] * (pow(theta22, 2) * params[3] + pow(theta11, 2) * params[2] * cos(theta1 - theta2))) / (params[2] * (2.0 * params[0] + params[1] - params[1] * cos(2.0 * theta1 - 2.0 * theta2))));
	double omp2 = ((2.0 * sin(theta1 - theta2) * (pow(theta11, 2) * params[2] * (params[0] + params[1]) + g * (params[0] + params[1]) * cos(theta1) + pow(theta22, 2) * params[3] * params[1] * cos(theta1 - theta2))) / (params[3] * (2.0 * params[0] + params[1] - params[1] * cos(2.0 * theta1 - 2.0 * theta2))));

	y[0] = y[1];
	y[1] = omp1;
	y[2] = y[3];
	y[3] = omp2;

	return y;
}
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include "Vector.h"

using namespace std;

Vector f(Vector y, double params[]);


int main() {

	double parameters[8] = { 0 };

	ifstream cParams("parameters.txt");
	if (!cParams) { //check if the parameters file exists

		cout << "The parameters file was not found in the directory and has been created with default values" << endl;
		ofstream wParams;
		wParams.open("parameters.txt"); //creates a blank parameters file and writes default values to is as follows
		wParams << "1.0 "; //m1
		wParams << "1.0 "; //m2
		wParams << "1.0 "; //l1
		wParams << "1.0 "; //l2
		wParams << "3.14159 "; //theta1_0
		wParams << "1.571 "; //theta2_0
		wParams << "0.001 "; //timestep h
		wParams << "10.0 "; //total duration t
		wParams.close();
	}

	ifstream fParams("parameters.txt");
	if (fParams.is_open()) {
		string value[] = { "mass 1","mass 2","length 1","length 2","theta 1 start", "theta 2 start","timestep","duration" };
		for (int i = 0; i < 8; i++) {
			fParams >> parameters[i]; //store values from the file into an array of 8 strings
			cout << "Read " << value[i] << " = " << parameters[i] << endl;
		}
		fParams.close();
	}

	else {
		cout << "Unable to open file";
	}

	ofstream resultsFile("output.txt"); //write results in overwritten output.txt file
	resultsFile << "t\tx1\ty1\tx2\ty2\n";
	resultsFile << "---------------------------------------\n";


	Vector y(4);
	y[0] = 1.571;
	y[1] = 0;
	y[2] = 2.356;
	y[3] = 0;

	double h = parameters[6];
	double n = parameters[7] / h;

	for (int i = 0; i < n; i++) {
		Vector k1 = f(y, parameters) * h;
		Vector k2 = f(y + k1 / 2, parameters) * h;
		Vector k3 = f(y + k2 / 2, parameters) * h;
		Vector k4 = f(y + k3, parameters) * h;

		y = y + ((k1 + k4) / 6 + (k2 + k3) / 3);
		double x1, y1, x2, y2;

		double l1 = parameters[2];
		double l2 = parameters[3];

		x1 = l1 * sin(y[0]);
		y1 = -l1 * cos(y[0]);
		x2 = x1 + l2 * sin(y[2]);
		y2 = y1 - l2 * cos(y[2]);

		cout << i << " " << x2 << " " << y2 << endl;
		resultsFile << i << "\t" << round(x1 * 10000) / 10000 << "\t" << round(y1 * 10000) / 10000 << "\t" << round(x2 * 10000) / 10000 << "\t" << round(y2 * 10000) / 10000 << endl;
	}
	resultsFile.close();
}

Vector f(Vector y, double params[]) {
	const double g = 9.81;
	double theta1 = y[0];
	double theta11 = y[1];
	double theta2 = y[2];
	double theta22 = y[3];

	double calc1 = ((-g * (2.0 * params[0] + params[1]) * sin(theta1) - params[1] * g * sin(theta1 - 2.0 * theta2) - 2.0 * sin(theta1 - theta2) * params[1] * (pow(theta22, 2) * params[3] + pow(theta11, 2) * params[2] * cos(theta1 - theta2))) / (params[2] * (2.0 * params[0] + params[1] - params[1] * cos(2.0 * theta1 - 2.0 * theta2))));
	double calc2 = ((2.0 * sin(theta1 - theta2) * (pow(theta11, 2) * params[2] * (params[0] + params[1]) + g * (params[0] + params[1]) * cos(theta1) + pow(theta22, 2) * params[3] * params[1] * cos(theta1 - theta2))) / (params[3] * (2.0 * params[0] + params[1] - params[1] * cos(2.0 * theta1 - 2.0 * theta2))));

	y[0] = y[1];
	y[1] = calc1;
	y[2] = y[3];
	y[3] = calc2;

	return y;
}
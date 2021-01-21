#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
using namespace std;

const double g_g = 9.81; //Global variable for gravitational constant

double f1(double t, double theta1, double theta2, double theta11, double theta22, double params[]);
double f2(double t, double theta1, double theta2, double theta11, double theta22, double params[]);
double f3(double t, double theta1, double theta2, double theta11, double theta22, double params[]);
double f4(double t, double theta1, double theta2, double theta11, double theta22, double params[]);


int main() {


    double parameters[8] = {0};

    ifstream cParams("parameters.txt");
    if (!cParams) { //check if the parameters file exists

        cout << "The parameters file was not found in the directory and has been created with default values" << endl;
        ofstream wParams;
        wParams.open("parameters.txt"); //creates a blank parameters file and writes default values to is as follows
        wParams << "1.5 "; //m1
        wParams << "2 "; //m2
        wParams << "1.0 "; //l1
        wParams << "1.0 "; //l2
        wParams << "1.571 "; //theta1_0
        wParams << "2.356 "; //theta2_0
        wParams << "0.01 "; //timestep h
        wParams << "10.0 "; //total duration t
        wParams.close();
    }

    ifstream fParams("parameters.txt");
    if (fParams.is_open()) {
        string value[] = {"mass 1","mass 2","length 1","length 2","theta 1 start", "theta 2 start","timestep","duration"};
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



    //Obtaining values from parameters
    double l1 = parameters[2];
    double l2 = parameters[3];
    double theta1 = parameters[4]; 
    double theta2 = parameters[5]; 
    double dTheta1 = 0;
    double dTheta2 = 0;
    double h = parameters[6];
    double n = parameters[7] / h;

    double k1_1, k1_2, k1_3, k1_4, k2_1, k2_2, k2_3, k2_4, k3_1, k3_2, k3_3, k3_4, k4_1, k4_2, k4_3, k4_4, x1, x2, y1, y2;

    double t = 0;

    //Runge Kutta 4 Method for a system of 4 Odes
    while (t < parameters[7]) {
        k1_1 = f2(t, theta1, theta2, dTheta1, dTheta2, parameters);
        k1_2 = f1(t, theta1, theta2, dTheta1, dTheta2, parameters);
        k1_3 = f4(t, theta1, theta2, dTheta1, dTheta2, parameters);
        k1_4 = f3(t, theta1, theta2, dTheta1, dTheta2, parameters);


        k2_1 = f2((t + h / 2), (theta1 + k1_2 * h / 2), (theta2 + k1_4 * h / 2), (dTheta1 + k1_1 * h / 2), (dTheta2 + k1_3 * h / 2), parameters);
        k2_2 = f1((t + h / 2), (theta1 + k1_2 * h / 2), (theta2 + k1_4 * h / 2), (dTheta1 + k1_1 * h / 2), (dTheta2 + k1_3 * h / 2), parameters);
        k2_3 = f4((t + h / 2), (theta1 + k1_2 * h / 2), (theta2 + k1_4 * h / 2), (dTheta1 + k1_1 * h / 2), (dTheta2 + k1_3 * h / 2), parameters);
        k2_4 = f3((t + h / 2), (theta1 + k1_2 * h / 2), (theta2 + k1_4 * h / 2), (dTheta1 + k1_1 * h / 2), (dTheta2 + k1_3 * h / 2), parameters);


        k3_1 = f2((t + h / 2), (theta1 + k2_2 * h / 2), (theta2 + k2_4 * h / 2), (dTheta1 + k2_1 * h / 2), (dTheta2 + k2_3 * h / 2), parameters);
        k3_2 = f1((t + h / 2), (theta1 + k2_2 * h / 2), (theta2 + k2_4 * h / 2), (dTheta1 + k2_1 * h / 2), (dTheta2 + k2_3 * h / 2), parameters);
        k3_3 = f4((t + h / 2), (theta1 + k2_2 * h / 2), (theta2 + k2_4 * h / 2), (dTheta1 + k2_1 * h / 2), (dTheta2 + k2_3 * h / 2), parameters);
        k3_4 = f3((t + h / 2), (theta1 + k2_2 * h / 2), (theta2 + k2_4 * h / 2), (dTheta1 + k2_1 * h / 2), (dTheta2 + k2_3 * h / 2), parameters);


        k4_1 = f2((t + h), (theta1 + k3_2 * h), (theta2 + k3_4 * h), (dTheta1 + k3_1 * h), (dTheta2 + k3_3 * h), parameters);
        k4_2 = f1((t + h), (theta1 + k3_2 * h), (theta2 + k3_4 * h), (dTheta1 + k3_1 * h), (dTheta2 + k3_3 * h), parameters);
        k4_3 = f4((t + h), (theta1 + k3_2 * h), (theta2 + k3_4 * h), (dTheta1 + k3_1 * h), (dTheta2 + k3_3 * h), parameters);
        k4_4 = f3((t + h), (theta1 + k3_2 * h), (theta2 + k3_4 * h), (dTheta1 + k3_1 * h), (dTheta2 + k3_3 * h), parameters);


        dTheta1 = dTheta1 + (k1_1 + k2_1 * 2 + k3_1 * 2 + k4_1) / 6.0 * h;
        theta1 = theta1 + (k1_2 + k2_2 * 2 + k3_2 * 2 + k4_2) / 6.0 * h;
        dTheta2 = dTheta2 + (k1_3 + k2_3 * 2 + k3_3 * 2 + k4_3) / 6.0 * h;
        theta2 = theta2 + (k1_4 + k2_4 * 2 + k3_4 * 2 + k4_4) / 6.0 * h;

        x1 = l1 * sin(theta1);
        y1 = -l1 * cos(theta1);
        x2 = x1 + l2 * sin(theta2);
        y2 = y1 - l2 * cos(theta2);


        
        cout << t << " " << x2 << " " << y2 << endl;
        resultsFile << t << "\t" << round(x1 * 10000) / 10000 << "\t" << round(y1 * 10000) / 10000 << "\t" << round(x2 * 10000) / 10000 << "\t" << round(y2 * 10000) / 10000 << endl; //Print results to file for all timesteps
        t += h;
    }


    resultsFile.close();

    return 0;
}


double f1(double t, double theta1, double theta2, double theta11, double theta22,double params[]) {
    return theta11;
}


double f2(double t, double theta1, double theta2, double theta11, double theta22, double params[]) {
    return ((-g_g * (2.0 * params[0] + params[1]) * sin(theta1) - params[1] * g_g * sin(theta1 - 2.0 * theta2) - 2.0 * sin(theta1 - theta2) * params[1] * (pow(theta22, 2) * params[3] + pow(theta11, 2) * params[2] * cos(theta1 - theta2))) / (params[2] * (2.0 * params[0] + params[1] - params[1] * cos(2.0 * theta1 - 2.0 * theta2))));
}


double f3(double t, double theta1, double theta2, double theta11, double theta22, double params[]) {
    return theta22;
}


double f4(double t, double theta1, double theta2, double theta11, double theta22, double params[]) {
    return((2.0 * sin(theta1 - theta2) * (pow(theta11, 2) * params[2] * (params[0] + params[1]) + g_g * (params[0] + params[1]) * cos(theta1) + pow(theta22, 2) * params[3] * params[1] * cos(theta1 - theta2))) / (params[3] * (2.0 * params[0] + params[1] - params[1] * cos(2.0 * theta1 - 2.0 * theta2))));
}

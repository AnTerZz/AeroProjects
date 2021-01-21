#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
using namespace std;


/* https://www.youtube.com/watch?v=JECaqGiKUXg
https://math.stackexchange.com/questions/1594617/difficulty-using-excel-to-solve-the-double-pendulum-problem-using-rk4-to-solve-f 
https://www.math24.net/double-pendulum/#:~:text=Substituting%20this%20in%20the%20original,l2%CE%B122.&text=dh%E2%88%82L%E2%88%82%CB%99%CE%B11%E2%88%92%E2%88%82,L%E2%88%82%CE%B12%3D0. */

const int g_g = 9.81;

double f1(double t, double theta1, double theta2, double theta11, double theta22, double params[]);
double f2(double t, double theta1, double theta2, double theta11, double theta22, double params[]);
double f3(double t, double theta1, double theta2, double theta11, double theta22, double params[]);
double f4(double t, double theta1, double theta2, double theta11, double theta22, double params[]);



int main() {


    double parameters[8]; //consider using a struct for code readability
    ifstream cParams("parameters.txt");

    if (!cParams) { //check if the parameters file exists

        cout << "The parameters file was not found in the directory and has been created with default values" << endl;
        ofstream wParams;
        wParams.open("parameters.txt"); //creates a blank parameters file and writes default values to is as follows
        wParams << "1.5 "; //m1
        wParams << "2 "; //m2
        wParams << "1.0 "; //l1
        wParams << "1.0 "; //l2
        wParams << "1.571 "; //theta10
        wParams << "2.356 "; //theta20
        wParams << "0.01 "; //Dt
        wParams << "10.0 "; //T
        wParams.close();
    }
    ifstream fParams("parameters.txt");

    if (fParams.is_open()) {

        for (int i = 0; i < 8; i++) {
            fParams >> parameters[i]; //store values from the file into an array of 8 strings
            cout << "Read parameter " << i << " : " << parameters[i] << endl;
        }
        fParams.close();
    }

    else {
        cout << "Unable to open file";
    }

    ofstream resultsFile("output.txt"); //write results in overwritten output.txt file
    resultsFile << "t\tx\ty\n";
    resultsFile << "-----------------------\n";



    //Runge Kutta 4 Method
    // x = Theta1 x1, u= Theta2 x2, y= dTheta1/h y1, v=dTheta2/h y2
    double m1 = parameters[0];
    double m2 = parameters[1];
    double l1 = parameters[2];
    double l2 = parameters[3];
    double theta1 = parameters[4];  //x0
    double theta2 = parameters[5];  //u0
    double dTheta1 = 0;
    double dTheta2 = 0;
    double h = parameters[6];
    double n = parameters[7] / h;

    double k11, k12, k13, k14, k21, k22, k23, k24, k31, k32, k33, k34, k41, k42, k43, k44, x2, y2;

    double t = 0;


    while (t <= 10) {
        k11 = f2(t, theta1, theta2, dTheta1, dTheta2, parameters);
        k12 = f1(t, theta1, theta2, dTheta1, dTheta2, parameters);
        k13 = f4(t, theta1, theta2, dTheta1, dTheta2, parameters);
        k14 = f3(t, theta1, theta2, dTheta1, dTheta2, parameters);


        k21 = f2((t + h / 2), (theta1 + k12 * h / 2), (theta2 + k14 * h / 2), (dTheta1 + k11 * h / 2), (dTheta2 + k13 * h / 2), parameters);
        k22 = f1((t + h / 2), (theta1 + k12 * h / 2), (theta2 + k14 * h / 2), (dTheta1 + k11 * h / 2), (dTheta2 + k13 * h / 2), parameters);
        k23 = f4((t + h / 2), (theta1 + k12 * h / 2), (theta2 + k14 * h / 2), (dTheta1 + k11 * h / 2), (dTheta2 + k13 * h / 2), parameters);
        k24 = f3((t + h / 2), (theta1 + k12 * h / 2), (theta2 + k14 * h / 2), (dTheta1 + k11 * h / 2), (dTheta2 + k13 * h / 2), parameters);


        k31 = f2((t + h / 2), (theta1 + k22 * h / 2), (theta2 + k24 * h / 2), (dTheta1 + k21 * h / 2), (dTheta2 + k23 * h / 2), parameters);
        k32 = f1((t + h / 2), (theta1 + k22 * h / 2), (theta2 + k24 * h / 2), (dTheta1 + k21 * h / 2), (dTheta2 + k23 * h / 2), parameters);
        k33 = f4((t + h / 2), (theta1 + k22 * h / 2), (theta2 + k24 * h / 2), (dTheta1 + k21 * h / 2), (dTheta2 + k23 * h / 2), parameters);
        k34 = f3((t + h / 2), (theta1 + k22 * h / 2), (theta2 + k24 * h / 2), (dTheta1 + k21 * h / 2), (dTheta2 + k23 * h / 2), parameters);


        k41 = f2((t + h), (theta1 + k32 * h), (theta2 + k34 * h), (dTheta1 + k31 * h), (dTheta2 + k33 * h), parameters);
        k42 = f1((t + h), (theta1 + k32 * h), (theta2 + k34 * h), (dTheta1 + k31 * h), (dTheta2 + k33 * h), parameters);
        k43 = f4((t + h), (theta1 + k32 * h), (theta2 + k34 * h), (dTheta1 + k31 * h), (dTheta2 + k33 * h), parameters);
        k44 = f3((t + h), (theta1 + k32 * h), (theta2 + k34 * h), (dTheta1 + k31 * h), (dTheta2 + k33 * h), parameters);


        dTheta1 = dTheta1 + (k11 + k21 * 2 + k31 * 2 + k41) / 6.0 * h;
        theta1 = theta1 + (k12 + k22 * 2 + k32 * 2 + k42) / 6.0 * h;
        dTheta2 = dTheta2 + (k13 + k23 * 2 + k33 * 2 + k43) / 6.0 * h;
        theta2 = theta2 + (k14 + k24 * 2 + k34 * 2 + k44) / 6.0 * h;


        x2 = l1*sin(theta1) + l2 * sin(theta2);
        y2 = -l1*cos(theta1) - l2 * cos(theta2);


        t += h;
        cout << t << " " << x2 << " " << y2 << endl;
        resultsFile << t << "\t" << x2 << "\t" << y2 << endl;
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

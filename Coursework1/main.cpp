#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
using namespace std;


/* https://www.youtube.com/watch?v=JECaqGiKUXg
https://math.stackexchange.com/questions/1594617/difficulty-using-excel-to-solve-the-double-pendulum-problem-using-rk4-to-solve-f 
https://www.math24.net/double-pendulum/#:~:text=Substituting%20this%20in%20the%20original,l2%CE%B122.&text=ddt%E2%88%82L%E2%88%82%CB%99%CE%B11%E2%88%92%E2%88%82,L%E2%88%82%CE%B12%3D0.*/

double f1(double T, double x, double u, double y, double v);
double f2(double T, double x, double u, double y, double v);
double f3(double T, double x, double u, double y, double v, string parameters[]);
double f4(double T, double x, double u, double y, double v, string parameters[]);

int main() {

    const double g = 9.81;



    string parameters[8]; //consider using a struct for code readability
    ifstream cParams("parameters.txt"); 

    if (!cParams) { //check if the parameters file exists

        cout << "The parameters file was not found in the directory and has been created with default values" << endl;
        ofstream wParams;
        wParams.open("parameters.txt"); //creates a bland parameters file and writes default values to is as follows
        wParams << "1.5 "; //m1
        wParams << "2.0 "; //m2
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
    // x = Theta1, u= Theta2, y= dTheta1/dt v=dTheta2/dt
    double ma1 = stof(parameters[0]);
    double ma2 = stof(parameters[1]);
    double l1 = stof(parameters[2]);
    double l2 = stof(parameters[3]);
    double x = stof(parameters[4]);  //x0
    double u = stof(parameters[5]);  //u0
    
    double T = 0;
    double v = 0;
    double y = 0;

/*#define f1(T,x,u,y,v) y
#define f2(T,x,u,y,v) v
#define f3(T,x,u,y,v) (-g*(2*ma1+ma2)*sin(x)-ma2*g*sin(x-2u)-2*sin(x-u)*ma2*(pow(v,2)*l2+pow(y,2)*l1*cos(x-u)))/(l1*(2*ma1+ma2-ma2*cos(2*x-2*u)))
#define f4(T,x,u,y,v) (2*sin(x-u)*(pow(y,2)*l1*(ma1+ma2)+g*(ma1+ma2)*cos(x)+pow(v,2)*l2*ma2*cos(x-u))) / (l2*(2*ma1+ma2-ma2*cos(2*x-2*u)))*/

    double h = stof(parameters[6]);
    double n = stof(parameters[7]) / h;


    resultsFile << T << "\t" << x << "\t" << u << endl;


    for (int i = 1; i < n; i++) {

        //cout << f1(T, x, u, y, v) << " " << f2(T, x, u, y, v) << " " << f3(T, x, u, y, v,parameters) << " " << f4(T, x, u, y, v,parameters) << endl;


        double k1 = h * (f1(T,x,u,y,v));
        double l1 = h * (f2(T, x, u, y, v));
        double m1 = h * (f3(T, x, u, y, v,parameters));
        double o1 = h * (f4(T, x, u, y, v,parameters));

        //cout << f1(T, x, u, y, v) << " " << f2(T, x, u, y, v) << " " << f3(T, x, u, y, v,parameters) << " " << f4(T, x, u, y, v,parameters) << endl;


        double k2 = h * (f1(T + 0.5 * h, x + 0.5 * h * k1, u + 0.5 * h * l1, y + 0.5 * h * m1, v + 0.5 * h * o1));
        double l2 = h * (f2(T + 0.5 * h, x + 0.5 * h * k1, u + 0.5 * h * l1, y + 0.5 * h * m1, v + 0.5 * h * o1));
        double m2 = h * (f3(T + 0.5 * h, x + 0.5 * h * k1, u + 0.5 * h * l1, y + 0.5 * h * m1, v + 0.5 * h * o1,parameters));
        double o2 = h * (f4(T + 0.5 * h, x + 0.5 * h * k1, u + 0.5 * h * l1, y + 0.5 * h * m1, v + 0.5 * h * o1,parameters));

        double k3 = h * (f1(T + 0.5 * h, x + 0.5 * h * k2, u + 0.5 * h * l2, y + 0.5 * h * m2, v + 0.5 * h * o2));
        double l3 = h * (f2(T + 0.5 * h, x + 0.5 * h * k2, u + 0.5 * h * l2, y + 0.5 * h * m2, v + 0.5 * h * o2));
        double m3 = h * (f3(T + 0.5 * h, x + 0.5 * h * k2, u + 0.5 * h * l2, y + 0.5 * h * m2, v + 0.5 * h * o2,parameters));
        double o3 = h * (f4(T + 0.5 * h, x + 0.5 * h * k2, u + 0.5 * h * l2, y + 0.5 * h * m2, v + 0.5 * h * o2,parameters));

        double k4 = h * (f1(T + h, x + h * k3, u + h * l3, y + h * m3, v + h * o3));
        double l4 = h * (f2(T + h, x + h * k3, u + h * l3, y + h * m3, v + h * o3));
        double m4 = h * (f3(T + h, x + h * k3, u + h * l3, y + h * m3, v + h * o3,parameters));
        double o4 = h * (f4(T + h, x + h * k3, u + h * l3, y + h * m3, v + h * o3,parameters));
        
        u = u + (k1 + k4) / 6 + (k2 + k3) / 3;
        x = x + (l1 + l4) / 6 + (l2 + l3) / 3;
        y = y + (m1 + m4) / 6 + (m2 + m3) / 3;
        v = v + (o1 + o4) / 6 + (o2 + o3) / 3;

        //cout << k1 << " " << k2 << " " << k3 << " " << k4 << endl; //returns ind after l2 - probably illegal opperation on a float
        resultsFile << T << "\t" << x << "\t" << u << endl;

        T = T + h;
    }


    resultsFile.close();

    return 0;
}

 double f1(double T,double x,double u,double y,double v) {
    double res = y;
    return res;
}

 double f2(double T, double x, double u, double y, double v) {
     double res = v;
     return res;
 }

 double f3(double T, double x, double u, double y, double v, string parameters[]) {
     const double g = 9.81;
     double ma1 = stof(parameters[0]);
     double ma2 = stof(parameters[1]);
     double l1 = stof(parameters[2]);
     double l2 = stof(parameters[3]);
     double res = (-g * (2 * ma1 + ma2) * sin(x) - ma2 * g * sin(x - 2u) - 2 * sin(x - u) * ma2 * (pow(v, 2) * l2 + pow(y, 2) * l1 * cos(x - u))) / (l1 * (2 * ma1 + ma2 - ma2 * cos(2 * x - 2 * u)));
     return res;
 }

 double f4(double T, double x, double u, double y, double v, string parameters[]) {
     const double g = 9.81;
     double ma1 = stof(parameters[0]);
     double ma2 = stof(parameters[1]);
     double l1 = stof(parameters[2]);
     double l2 = stof(parameters[3]);
     double res = (2 * sin(x - u) * (pow(y, 2) * l1 * (ma1 + ma2) + g * (ma1 + ma2) * cos(x) + pow(v, 2) * l2 * ma2 * cos(x - u))) / (l2 * (2 * ma1 + ma2 - ma2 * cos(2 * x - 2 * u)));
     return res;
 }
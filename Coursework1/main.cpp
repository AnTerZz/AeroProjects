#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int main() {


    string parameters[8]; //consider using a struct for code readability
    ifstream cParams("parameters.txt"); 

    if (!cParams) { //check if the parameters file exists

        cout << "The parameters file does not exist and has been instanciated with default values" << endl;
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
    resultsFile.close();

    return 0;
}



double  m1 = 1.0, m2 = 1.0, g = 9.81, L1 = 1.0, L2 = 1.0;



double fx1(double t, double theta1, double theta2, double theta11, double theta22) {

    return theta11;
}



double fv1(double t, double theta1, double theta2, double theta11, double theta22) {

    return ((-g * (2.0 * m1 + m2) * sin(theta1) - m2 * g * sin(theta1 - 2.0 * theta2) - 2.0 * sin(theta1 - theta2) * m2 * (pow(theta22, 2) * L2 + pow(theta11, 2) * L1 * cos(theta1 - theta2))) / (L1 * (2.0 * m1 + m2 - m2 * cos(2.0 * theta1 - 2.0 * theta2))));
}



double fx2(double t, double theta1, double theta2, double theta11, double theta22) {

    return theta22;
}



double fv2(double t, double theta1, double theta2, double theta11, double theta22) {

    return((2.0 * sin(theta1 - theta2) * (pow(theta11, 2) * L1 * (m1 + m2) + g * (m1 + m2) * cos(theta1) + pow(theta22, 2) * L2 * m2 * cos(theta1 - theta2))) / (L2 * (2.0 * m1 + m2 - m2 * cos(2.0 * theta1 - 2.0 * theta2))));
}



int main()

{


    ofstream resultsFile("output2.txt"); //write results in overwritten output2.txt file
    resultsFile << "t\tx\ty\n";
    resultsFile << "-----------------------\n";
    double t, dt, L1, L2, t1, t2, td1, td2, k1x1, k1x2, k1v1, k1v2, k2x1, k2x2, k2v1, k2v2, k3x1, k3x2, k3v1, k3v2, k4x1, k4x2, k4v1, k4v2, x1, x2, y1, y2;



    t = 0;

    dt = 0.001;

    L1 = 1.0;

    L2 = 1.0;

    t1 = 3.14159;

    t2 = 3.14159 / 2;

    td1 = 0;

    td2 = 0;



    for (double t = 0; t <= 10; t += dt)

    {

        k1x1 = fv1(t, t1, t2, td1, td2);

        k1x2 = fx1(t, t1, t2, td1, td2);

        k1v1 = fv2(t, t1, t2, td1, td2);

        k1v2 = fx2(t, t1, t2, td1, td2);



        k2x1 = fv1((t + dt / 2), (t1 + k1x2 * dt / 2), (t2 + k1v2 * dt / 2), (td1 + k1x1 * dt / 2), (td2 + k1v1 * dt / 2));

        k2x2 = fx1((t + dt / 2), (t1 + k1x2 * dt / 2), (t2 + k1v2 * dt / 2), (td1 + k1x1 * dt / 2), (td2 + k1v1 * dt / 2));

        k2v1 = fv2((t + dt / 2), (t1 + k1x2 * dt / 2), (t2 + k1v2 * dt / 2), (td1 + k1x1 * dt / 2), (td2 + k1v1 * dt / 2));

        k2v2 = fx2((t + dt / 2), (t1 + k1x2 * dt / 2), (t2 + k1v2 * dt / 2), (td1 + k1x1 * dt / 2), (td2 + k1v1 * dt / 2));



        k3x1 = fv1((t + dt / 2), (t1 + k2x2 * dt / 2), (t2 + k2v2 * dt / 2), (td1 + k2x1 * dt / 2), (td2 + k2v1 * dt / 2));

        k3x2 = fx1((t + dt / 2), (t1 + k2x2 * dt / 2), (t2 + k2v2 * dt / 2), (td1 + k2x1 * dt / 2), (td2 + k2v1 * dt / 2));

        k3v1 = fv2((t + dt / 2), (t1 + k2x2 * dt / 2), (t2 + k2v2 * dt / 2), (td1 + k2x1 * dt / 2), (td2 + k2v1 * dt / 2));

        k3v2 = fx2((t + dt / 2), (t1 + k2x2 * dt / 2), (t2 + k2v2 * dt / 2), (td1 + k2x1 * dt / 2), (td2 + k2v1 * dt / 2));



        k4x1 = fv1((t + dt), (t1 + k3x2 * dt), (t2 + k3v2 * dt), (td1 + k3x1 * dt), (td2 + k3v1 * dt));

        k4x2 = fx1((t + dt), (t1 + k3x2 * dt), (t2 + k3v2 * dt), (td1 + k3x1 * dt), (td2 + k3v1 * dt));

        k4v1 = fv2((t + dt), (t1 + k3x2 * dt), (t2 + k3v2 * dt), (td1 + k3x1 * dt), (td2 + k3v1 * dt));

        k4v2 = fx2((t + dt), (t1 + k3x2 * dt), (t2 + k3v2 * dt), (td1 + k3x1 * dt), (td2 + k3v1 * dt));



        td1 = td1 + (k1x1 + k2x1 * 2 + k3x1 * 2 + k4x1) / 6.0 * dt;

        t1 = t1 + (k1x2 + k2x2 * 2 + k3x2 * 2 + k4x2) / 6.0 * dt;

        td2 = td2 + (k1v1 + k2v1 * 2 + k3v1 * 2 + k4v1) / 6.0 * dt;

        t2 = t2 + (k1v2 + k2v2 * 2 + k3v2 * 2 + k4v2) / 6.0 * dt;



        x1 = L1 * sin(t1);

        x2 = x1 + L2 * sin(t2);

        y1 = -L1 * cos(t1);

        y2 = y1 - L2 * cos(t2);



        cout << t << "\t" << x1 << "\t" << y1 << "\t" << x2 << "\t" << y2 << "\t" << t1 << "\t" << t2 << "\t" << td1 << "\t" << td2 << endl;
        resultsFile << t << "\t" << x2 << "\t" << y2 << endl;

    }



    return 0;

}
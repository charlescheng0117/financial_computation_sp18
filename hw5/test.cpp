#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

bool greater_double (double d1, double d2) {
    return d1 > d2;
}


int main(int argc, char *argv[]) {
    
    vector<double> v = { 57.440422, 56.572552, 55.717796, 54.875954, 54.046831, 53.230235 };

    vector<double>::iterator low, up;

    low = lower_bound(v.begin(), v.end(), 57.440422, greater_double);
    up = upper_bound(v.begin(), v.end(), 57.440422, greater_double);

    cout << "lower_bound at position " << (low- v.begin()) << '\n';
    cout << "upper_bound at position " << (up - v.begin()) << '\n';

    


    return 0;
}

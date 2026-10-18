

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iostream>

/*
   g++ -std=c++17 -I /home/user/boost_1_85_0 main.cpp -o ublas_example
   g++ -std=c++17 -I 
*/

int main() {
    using namespace boost::numeric::ublas;
    vector<double> v(1024);
    for (unsigned i = 0; i < v.size(); ++i) {
        v(i) = static_cast<double>(i);
    }
    std::cout << v << std::endl;
    return 0;
}
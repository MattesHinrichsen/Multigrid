#include <iostream>
#include <chrono>
#include <cmath>

class Timer {
    std::chrono::_V2::system_clock::time_point startTime;
    std::chrono::_V2::system_clock::time_point endTime;
    static int counter;
public:
    Timer() {
        startTime = std::chrono::high_resolution_clock::now();
        counter++;
    }
    ~Timer() {
        endTime = std::chrono::high_resolution_clock::now();
        auto differenz = endTime-startTime;
        counter--;

        auto minutes = std::chrono::duration_cast< std::chrono::minutes >( differenz );
        differenz -= minutes;

        auto seconds = std::chrono::duration_cast< std::chrono::seconds >( differenz );
        differenz -= seconds;

        auto milliseconds = std::chrono::duration_cast< std::chrono::milliseconds >( differenz );
        differenz -= milliseconds;

        auto microseconds = std::chrono::duration_cast< std::chrono::microseconds >( differenz );
        differenz -= microseconds;

        auto nanoseconds = std::chrono::duration_cast< std::chrono::nanoseconds >( differenz );
        std::cout << "Timer " << counter << " finished in: ";
        std::cout << minutes.count() << " Minutes, "
                    << seconds.count() << " Seconds, "
                    << milliseconds.count() << " Milliseconds, "
                    << microseconds.count() << " Microseconds, "
                    << nanoseconds.count() << " Nanoseconds." << std::endl;

    }
};   int Timer::counter = 1;


// int main() {
    
//     Timer a;
//     double x = 0;
//     int less;
//     std::cin >> less;
//     for (int i = 0; i<less; i++) {
//         x = std::sqrt(i);
//     }
//     std::cout << x << std::endl;

    
// }





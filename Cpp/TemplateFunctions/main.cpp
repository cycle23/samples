#include <iostream>
#include <array>
#include <vector>
#include <forward_list>
using namespace std;

template<typename T>
void f(const T& param) {
    cout << param << endl;
}

template<typename T>
void refParam(T& param) {
}

template<typename T>
void constRefParam(const T& param) {
}

template<typename T>
void universalRef(T&& param) {

}

template<typename T>
void valueRef(T param) {

}

template<typename T>
class TD;




// return size of an array as a compile-time constant. (The
// array parameter has no name, because we care only about
// the number of elements it contains.)
template<typename T, size_t N>
constexpr size_t arraySize(T(&)[N]) noexcept
{
    return N;
}

template<typename T>
std::initializer_list<T> createInitList()
{
    return { 1, 2, 3};
}

template<typename T>
class Reset {

    void resetIt() {
        vector<T> v;
        auto resetV =
                [&v](const initializer_list<T> &newValue) { v = newValue; };
        resetV({1, 2, 3});
    }
};

void authenticateUser()
{
}

template<typename Container, typename Index>
auto authAndAccess(Container&& c, Index i) // universal ref, could use decltype(auto) here for C++14
    -> decltype(forward<Container>(c)[i]) // C++11 requirement
{
    authenticateUser();
    return forward<Container>(c)[i];
}


int main(int argc, char **argv)
{
    int x = 27;
    f(x);
    const int cx = x;
    const int& rx = x;
    refParam(x);  // T = int, param = int&
    refParam(cx); // T = const int, param = const int&
    refParam(rx); // T = const int, param = const int&
    constRefParam(x); // T = int, param = const int&
    constRefParam(cx);// T = int, param = const int&
    constRefParam(rx);// T = int, param = const int&

    universalRef(x);  // x is lvalue, T is int&, param int&
    universalRef(cx); // cx lvalue, T/param const int&
    universalRef(rx);// rx lvalue, T/param const int&
    universalRef(27);// 27 is rvalue, so T is int, param's type is int&&

    valueRef(x);  // simply T/param = int, no const or ref for all 3
    valueRef(cx);
    valueRef(rx);

    int keyVals[] = { 1, 3, 5, 7, 9, 11 };

    int mappedVals[arraySize(keyVals)];
    array<int, arraySize(keyVals)> mappedVals2;

    TD<decltype(universalRef(27))> type27;

    return 0;
}


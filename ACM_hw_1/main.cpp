#include <algorithm>
#include <functional>
#include <array>
#include <iostream>
#include <iterator>

//template< class RandomIt, class Compare >
//RandomIt MyPartition(RandomIt first, RandomIt last, Compare comp) {
//    auto pivot = *std::next(first, std::distance(first,last)/2);
//    for(RandomIt i = first; i != last; i++)
//        if (comp(*i, pivot)) {
//            std::iter_swap(i, first);
//            first++;
//        }
//    return first;
//}

template<class ForwardIt, class UnaryPredicate>
ForwardIt MyPartition(ForwardIt first, ForwardIt last, UnaryPredicate p)
{
    first = std::find_if_not(first, last, p);
    if (first == last) return first;

    for(ForwardIt i = std::next(first); i != last; ++i){
        if(p(*i)){
            std::iter_swap(i, first);
            ++first;
        }
    }
    return first;
}

template< class RandomIt, class Compare >
void MyNthElement( RandomIt first, RandomIt nth, RandomIt last, Compare comp )
{
    if(first == last) return;
    auto pivot = *std::next(first, std::distance(first,last)/2);
    RandomIt q = MyPartition(first, last, [pivot, comp] (const auto &a) { return comp(a, pivot);});
    if (q == nth)
        return;
    else
        if (comp(*q, *nth))
            MyNthElement(first, nth, q, comp);
        else
            MyNthElement(++q, nth, last, comp);
    return;
}

int main()
{
    std::array<int, 10> s{7, 3, 4, 2, 8, 1, 6, 9, 0, 5};

//    MyPartition(s.begin(), s.end(), [](auto a, auto b){ return a < b; });
//    std::copy(s.begin(), s.end(), std::ostream_iterator<int>(std::cout," "));
//    std::cout << std::endl;

    MyNthElement(s.begin(), s.begin() + s.size()/3, s.end(), [](auto a, auto b){ return a < b; });
    std::cout << s[s.size()/3];
    std::cout << std::endl;
}
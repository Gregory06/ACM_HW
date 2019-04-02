// Реализуйте следующие методы:

// 1. MyMerge, работающий аналогично std::merge
// https://en.cppreference.com/w/cpp/algorithm/merge
// ключевое слово "constexpr" можете исключить из объявления, если не хотите с этим разбираться
// Требуемая сложность O(N)
template< class InputIt1, class InputIt2, class OutputIt >
constexpr OutputIt MyMerge( InputIt1 first1, InputIt1 last1,
                          InputIt2 first2, InputIt2 last2,
                          OutputIt d_first )
{
    for(; first1 != last1; d_first++) {
        if (first2 == last2)
            return std::copy(first1, last1, d_first);
        if (*first1 < *first2) {
            *d_first = *first1;
            first1++;
        } else
        {
            *d_first = *first2;
            first2++;
        }
    }
    return std::copy(first2, last2, d_first);
}

// 2. MySort, как аналог std::sort, использующий для своей работы сортировку слиянием
// https://ru.cppreference.com/w/cpp/algorithm/sort
// Требуемая сложность O(NlogN)
template< class RandomIt, class Compare >
void MySort( RandomIt first, RandomIt last, Compare comp )
{
    auto pivot = *first;
    if (first != last) {
        RandomIt separator = MyPartition(first, last, [pivot, comp] (const auto &a) { return comp(a, pivot); });
        MySort(first, separator++, comp);
        MySort(separator, last, comp);
    }
}

// 3. MyPartition, работающий так же, как и std::partition (не обращайте внимание на параметр ExecutionPolicy)
// https://ru.cppreference.com/w/cpp/algorithm/partition
// RandomIt - некоторый удобный для вас тип итераторов, я бы не рекомендовал запариваться над его выбором, используйте Bidirectional
// Требуемая сложность O(N)
template< class RandomIt, class UnaryPredicate >
RandomIt MyPartition(RandomIt first, RandomIt last, UnaryPredicate p) {
    for(RandomIt i = first; i != last; i++)
        if (p(*i)) {
            std::iter_swap(i, first);
            first++;
        }
    return first;
}

// 4. MyNthElement, работающий аналогично std::nth_element 
// https://en.cppreference.com/w/cpp/algorithm/nth_element
// RandomIt - некоторый удобный для вас тип итераторов, я бы не рекомендовал запариваться, используйте Bidirection
// Требуемая сложность O(N), внутренная реализация - линейный алгоритм, рассказанный на семинаре.
template< class RandomIt, class Compare >
void MyNthElement( RandomIt first, RandomIt nth, RandomIt last, Compare comp )
{
    if(first == last) return;
    auto pivot = *std::next(first, std::distance(first,last)/2);
    RandomIt q = MyPartition(first, last, [pivot, comp] (const auto &a) { return comp(a, pivot); });
    if (q == nth)
        return;
    else
    if (comp(*q, *nth))
        MyNthElement(first, nth, q, comp);
    else
        MyNthElement(++q, nth, last, comp);
    return;
}

// В случае, если у вас нет достаточного знания языка С++, можете реализовать все методы в стиле C, самостоятельно 
// выбрав способы передачи и возврата аргументов, но семинаристы вам будут не очень благодарны :-)

// Решение - файл hw1_practic_group_surname_name.hpp, на почту amvshe4ka2019@gmail.com, тема письма должна иметь вид 
// hw1_practic_group_surname_name, например hw1_practic_ivanov_alexandr_573

// В данном случае специально не используем автоматические системы проверки для того, чтобы вы во-первых, 
// самостоятельно протестировали работу методов, и во-вторых, мы попробуем дать фидбек по стилю кода

// Если вы умеете пользоваться такими фреймворками как catch или google_test, то делайте это :-)
// Ваш код будет компилироваться с ключами g++-8 -Wall -Werror -Wextra -Wpedantic -std=c++17 -O2 -g -fsanitize=address,leak,undefined
// В случае использования C сами укажите желаемый стандарт языка и версию компилятора

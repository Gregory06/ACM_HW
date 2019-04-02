#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <complex>
#include <algorithm>
#include <iomanip>
#include <complex>
#include <initializer_list>

const double PI = std::acos(-1);
const double eps = 0.0000001;

// Реализуте пропущенные методы
// в качестве Т будет использоваться std::complex<float> // <double> // <long double>
namespace FFT {
    // Возвращает корень степени degree из 1
    template <typename T>
    T GetRoot(size_t degree);

    // Выполняет преобразование фурье квадратичным алгоритмом для вектора произвольной длины
    template <typename T>
    std::vector<T> FourierTransform(const std::vector<T>& data);

    // Выполняет обратное преобразование квадратичным алгоритмом для вектора фиксированной длины
    template <typename T>
    std::vector<T> InverseFourierTransform(const std::vector<T>& data);

    // Добивает вектор в конце нулями до длины expected_length,
    // выбрасывает std::runtime_error если expected_length < data.size()
    template <typename T>
    std::vector<T> AddPadding(const std::vector<T>& data, size_t expected_length);

    // Быстрое преобразование Фурье для вектора длины 2^k
    template <typename T>
    std::vector<T> FastFourierTransform(const std::vector<T>& data);

    // Обратное быстрое преобразование Фурье для вектора длины 2^k
    template <typename T>
    std::vector<T> FastInverseFourierTransform(const std::vector<T>& data);
} // namespace FFT


// Операции над многочленами с помощью ффт
template <typename T>
class Polynomial {
public:
    explicit Polynomial(const std::vector<T>& coefficients)
        : coefficients_(coefficients), degree_(coefficients.size()) {
    }

    // Чтобы можно было написать Polynomial<std::complex<double>> p({1, 2, 1})
    // И получить представление многочлена 1 + 2x + x^2
    Polynomial(const std::initializer_list<T>& coefficients)
        : coefficients_(coefficients.begin(), coefficients.end()), degree_(coefficients.size()) {
    }

    Polynomial(const Polynomial&) = default;
    Polynomial(Polynomial&&) noexcept = default;

//    Polynomial& operator=(const Polynomial) = default;
//    Polynomial& operator=(Polynomial&&) noexcept = default;

    size_t GetDegree() const {
        return degree_;
    }

    // Реализуйте следующие методы
    bool operator==(const Polynomial& other);

    Polynomial& clean();

    Polynomial& operator+=(const Polynomial& other);

    Polynomial& operator-=(const Polynomial& other);

    // Возведение в степень pow с помощью комбинации FFT и индийского возведения в степень
    Polynomial& operator^=(size_t pow);

    // С Использованием быстрого преобразования фурье
    Polynomial& operator*=(const Polynomial& other);

    friend std::ostream& operator<<(std::ostream& ostr, const Polynomial& a);

    // Используя предыдущее
    Polynomial operator+(const Polynomial& other);

    Polynomial operator-(const Polynomial& other);

    Polynomial operator*(const Polynomial& other);

    Polynomial operator^(size_t pow);

    // И еще один, унарный минус
    friend Polynomial operator-(const Polynomial& other);

private:
    std::vector<T> coefficients_;
    size_t degree_;
};


// Напишите какие-то тесты, демонстрирующие корректность написанных вами методов



// Задачи, решаемые с помощью умножения многочленов
// Если вы напишете решение, работающее для произольных строк над ascii кодировкой - укажете это и
// возможно, получите небольшой бонусный балл
namespace SubstringMatching {

    // Метод принимает две строки str и pattern, возвращает индексов,
    // указывающих начала подстрок str, совпадающих с pattern
    // Можете считать, что str и pattern состоят из символов 0 и 1
    std::vector<size_t> FindSubstrings(const std::string& str, const std::string& pattern);

    // Аналогично предыдущему, но теперь pattern может содержать символ '?', на месте которого
    // можно поставить любой символ алфавита
    std::vector<size_t> FindMatches(const std::string& str, const std::string& pattern);

} // namespace SubstringMatcher



// Напишите какие-то тесты, демонстрирующие корректность написанных вами методов
// Обратите внимание, что при проведении прямого и обратного преобразования для многочлена с
// целыми коэффициентами, новые коэффициенты могут выйти дробными и трубующими округления.
template <typename T>
T FFT::GetRoot(size_t degree) {
    return T(cos(2*PI / degree), sin(2*PI / degree));
}

template <typename T>
std::vector<T> FFT::FourierTransform(const std::vector <T> &data) {
    size_t n(data.size());
    std::vector<T> result(n);
    T kth_component(0);
    std::complex<double> root(GetRoot<T>(n));
    for (size_t i = 0; i < n; i++) {
        kth_component = 0;
        for (size_t j = 0; j < n; j++) {
            kth_component += pow(root, i*j) * data[j];
        }
        result[i] = kth_component;
    }
    return result;
}

template <typename T>
std::vector<T> FFT::InverseFourierTransform(const std::vector <T> &data) {
    size_t n(data.size());
    std::vector<T> result(n);
    T kth_component(0);
    T root(GetRoot<T>(n));
    for (size_t i = 0; i < n; i++) {
        kth_component = 0;
        for (size_t j = 0; j < n; j++)
            kth_component += pow(root, (int) -(i*j)) * data[j];
        result[i] = kth_component / T(n);
    }
    return result;
}

template <typename T>
std::vector<T> FFT::FastFourierTransform(const std::vector<T>& data){
    size_t n(data.size());
    std::vector<T> result(n);
    if (n == 1)
        return data;

    std::vector<T> first_half(n/2), second_half(n/2);
    for (size_t i = 0; i < n/2 ; i++) {
        first_half[i] = data[2*i];
        second_half[i] = data[2*i+1];
    }

    first_half = FastFourierTransform(first_half);
    second_half = FastFourierTransform(second_half);

    T w(1), root(GetRoot<T>(n));
    for (size_t i = 0; i < n/2; i++) {
        result[i] = first_half[i] + w * second_half[i];
        result[i+n/2] = first_half[i] - w * second_half[i];

        w *= root;
    }

    return result;
}

template <typename T>
std::vector<T> FFT::FastInverseFourierTransform(const std::vector<T>& data){
    size_t n(data.size());
    std::vector<T> result(n);
    if (n == 1)
        return data;

    std::vector<T> first_half(n/2), second_half(n/2);
    for (size_t i = 0; i < n/2; i++) {
        first_half[i] = data[2*i];
        second_half[i] = data[2*i+1];
    }

    first_half = FastInverseFourierTransform(first_half);
    second_half = FastInverseFourierTransform(second_half);

    T w(1), root(pow(GetRoot<T>(n), n-1));
    for (size_t i = 0; i < n/2; i++) {
        result[i] = (first_half[i] + w * second_half[i]) / T(2,0);
        result[i+n/2] = (first_half[i] - w * second_half[i]) / T(2,0);

        w *= root;
    }

    return result;
}

template <typename T>
std::vector<T> FFT::AddPadding(const std::vector<T>& data, size_t expected_length) {
    std::vector<T> result(data);
    for (size_t i = data.size(); i < expected_length; i++)
        result.push_back(T(0,0));
    return result;
}

template <typename T>
Polynomial<T>& Polynomial<T>::clean() {
    for (size_t i = coefficients_.size()-1; i != -1; i--) {
        if (abs(coefficients_[i]) < eps) {
            coefficients_.pop_back();
            continue;
        }

        if (coefficients_[i].real() < eps)
            coefficients_[i] = T(0, coefficients_[i].imag());
        if (coefficients_[i].imag() < eps)
            coefficients_[i] = T(coefficients_[i].real(), 0);

    }
    degree_ = coefficients_.size();

    return *this;
}

template <typename T>
Polynomial<T>& Polynomial<T>::operator*=(const Polynomial<T>& other) {
    size_t expected_length = (size_t) pow(2, (size_t) std::ceil(log(this->degree_ + other.degree_)/log(2)));
    std::vector<T> added_this = FFT::AddPadding(this->coefficients_, expected_length);
    std::vector<T> added_other = FFT::AddPadding(other.coefficients_, expected_length);
    std::vector<T> fft_this = FFT::FourierTransform(added_this);
    std::vector<T> fft_other = FFT::FourierTransform(added_other);

//    for (int i = 0; i < fft_this.size(); i++)
//        std::cout << fft_this[i] << ' ';
//    std::cout << std::endl;
//    for (int i = 0; i < fft_other.size(); i++)
//        std::cout << fft_other[i] << ' ';
//    std::cout << std::endl;

    for (int i = 0; i < expected_length; i++)
        fft_this[i] *= fft_other[i];

    fft_this = FFT::AddPadding(fft_this, expected_length);
    std::vector<T> result = FFT::FastInverseFourierTransform(fft_this);

    coefficients_ = result;
    degree_ = result.size();

    this->clean();

    return *this;
}

template <typename T>
bool Polynomial<T>::operator==(const Polynomial<T>& other) {
    if (this->degree_ != other.degree_)
        return false;
    double difference = 0;
    for (size_t i = 0; i < other.coefficients_.size(); i++)
        difference += abs(this->coefficients_[i] - other.coefficients_[i]);

    return difference < eps;
}

template <typename T>
Polynomial<T>& Polynomial<T>::operator+=(const Polynomial& other) {
    if (degree_ < other.degree_)
        coefficients_.resize(other.degree_);
    for (size_t i = 0; i < other.degree_; i++) {
        coefficients_[i] += other.coefficients_[i];
    }
    degree_ = coefficients_.size();

    return *this;
}

template <typename T>
Polynomial<T>& Polynomial<T>::operator-=(const Polynomial& other) {
    if (degree_ < other.degree_)
        coefficients_.resize(other.degree_);
    for (size_t i = 0; i < other.degree_; i++) {
        coefficients_[i] -= other.coefficients_[i];
    }
    degree_ = coefficients_.size();

    return *this;
}

template <typename T>
Polynomial<T>& Polynomial<T>::operator^=(size_t pow) {
    Polynomial<T> helper({1});
    for (int j = sizeof(size_t) * 8 - 1; j > 0; --j) {
        if ((pow >> j) & 1) {
            helper *= *this;
            helper *= helper;
        } else {
            helper *= helper;
        }
    }

    if (pow % 2)
        helper *= *this;

    coefficients_ = helper.coefficients_;
    degree_ = helper.degree_;

    this->clean();

    return *this;
}

template <typename T>
Polynomial<T> Polynomial<T>::operator+(const Polynomial<T>& other) {
    Polynomial<T> result (*this);
    result += other;

    return result;
}

template <typename T>
Polynomial<T> Polynomial<T>::operator-(const Polynomial<T>& other) {
    Polynomial<T> result (*this);
    result -= other;

    return result;
}

template <typename T>
Polynomial<T> Polynomial<T>::operator*(const Polynomial<T>& other) {
    Polynomial<T> result (*this);
    result *= other;

    return result;
}

template <typename T>
Polynomial<T> Polynomial<T>::operator^(size_t pow) {
    Polynomial<T> result (*this);
    result ^= pow;

    return result;
}

std::ostream& operator<<(std::ostream& ostr, const Polynomial<std::complex<double>>& a) {
    for (std::vector<std::complex<double>>::const_iterator i = a.coefficients_.begin(); i != a.coefficients_.end(); ++i)
        ostr << *i << ' ';
    return ostr;
}

int main() {

//    std::cout << pow(FFT::GetRoot<std::complex<double>>(4), -2);
//    Polynomial<std::complex<double>> p(std::vector<std::complex<double>>{1, 2, 1});
    std::vector<std::complex<double>> a = std::vector<std::complex<double>>{1,2,1,0,5,6,7,8}, b = std::vector<std::complex<double>>{4,8,4,12};
//    std::vector<std::complex<double>> a = FFT::AddPadding<std::complex<double>>(std::vector<std::complex<double>>{1,1}, pow(2, (size_t) std::ceil(log(4)/log(2))));
//    for (std::vector<std::complex<double>>::const_iterator i = a.begin(); i != a.end(); ++i)
//        std::cout << *i << ' ';
//    std::vector<std::complex<double>> b = FFT::FastFourierTransform(a);
    std::vector<std::complex<double>> c = FFT::FastInverseFourierTransform(a);

//    std::vector<std::complex<double>> c = FFT::FastInverseFourierTransform(std::vector<std::complex<double>>{1,1});
    Polynomial<std::complex<double>> p1({1, 2}), p2({1, 1});

    Polynomial<std::complex<double>> p3 = p1 ^ 3;

    std::cout << p3;
//    std::cout << std::endl;
//    for (std::vector<std::complex<double>>::const_iterator i = c.begin(); i != c.end(); ++i)
//        std::cout << *i << ' ';

//    c = FFT::FastInverseFourierTransform(FFT::FastFourierTransform(a));
//    std::cout << std::endl;
//    for (std::vector<std::complex<double>>::const_iterator i = c.begin(); i != c.end(); ++i)
//        std::cout << *i << ' ';
//    std::cout << std::endl;
//
//    c = FFT::InverseFourierTransform(a);
//    std::cout << std::endl;
//    for (std::vector<std::complex<double>>::const_iterator i = c.begin(); i != c.end(); ++i)
//        std::cout << *i << ' ';
//    std::cout << std::endl;

//    c = FFT::FastInverseFourierTransform(b);
//    std::cout << std::endl;
//    for (std::vector<std::complex<double>>::const_iterator i = c.begin(); i != c.end(); ++i)
//        std::cout << *i << ' ';
}
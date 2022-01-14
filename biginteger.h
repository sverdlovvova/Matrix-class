#include <iostream>

const int BASE = 1000'000'000;
const int BASE_CNT = 9;

class BigInteger;
class Rational;

BigInteger operator+(const BigInteger &first, const BigInteger &second);
BigInteger operator-(const BigInteger &first, const BigInteger &second);
BigInteger operator*(const BigInteger &first, const BigInteger &second);
BigInteger operator/(const BigInteger &first, const BigInteger &second);
BigInteger operator%(const BigInteger &first, const BigInteger &second);
bool operator==(const BigInteger &first, const BigInteger &second);
bool operator!=(const BigInteger &first, const BigInteger &second);
bool operator<(const BigInteger &first, const BigInteger &second);
bool operator>(const BigInteger &first, const BigInteger &second);
bool operator<=(const BigInteger &first, const BigInteger &second);
bool operator>=(const BigInteger &first, const BigInteger &second);
BigInteger GCD(BigInteger first, BigInteger second);

class BigInteger {
private:
    std::vector<int> digits;
    size_t size = 0;
    int sign = 1;

public:
    BigInteger() = default;

    BigInteger(const int &number) : digits(0), size(0), sign(1) {
        long long longNumber = number;
        if (number < 0) {
            sign = -1;
            longNumber *= -1;
        }
        while (longNumber > 0) {
            digits.push_back(int(longNumber % BASE));
            longNumber /= BASE;
        }
        size = digits.size();
    }

    std::string toString() const {
        std::string result;
        for (size_t i = 0; i < size; ++i) {
            int digit = digits[i];
            int cnt = 0;
            while (digit > 0) {
                digit /= 10;
                ++cnt;
            }
            digit = digits[i];
            while (digit > 0) {
                result += char(digit % 10 + '0');
                digit /= 10;
            }
            if (i != size - 1) {
                for (int j = cnt; j < BASE_CNT; ++j) {
                    result += '0';
                }
            }
        }
        if (sign == -1) {
            result += '-';
        }
        for (size_t i = 0; i < result.size() / 2; ++i) {
            std::swap(result[i], result[result.size() - i - 1]);
        }
        if (result.empty()) {
            result += '0';
        }
        return result;
    }

    BigInteger operator-() const {
        BigInteger result = *this;
        if (result.size > 0) result.sign *= -1;
        return result;
    }

    BigInteger &operator+=(const BigInteger &another) {
        BigInteger first = *this;
        BigInteger second = another;
        first.sign = 1;
        second.sign = 1;
        bool firstIsBigger = (first >= second);
        BigInteger& big = firstIsBigger ? first : second;
        BigInteger& little = firstIsBigger ? second : first;
        digits.resize(big.size);
        little.digits.resize(big.size);
        if (sign != another.sign) {
            int subtract = 0;
            for (size_t i = 0; i < big.size; ++i) {
                int digit = big.digits[i] - little.digits[i] - subtract;
                subtract = 0;
                if (digit < 0) {
                    digit += BASE;
                    subtract = 1;
                }
                digits[i] = digit;
            }
            sign = (firstIsBigger && sign == -1) || (!firstIsBigger && another.sign == -1) ? -1 : 1;
            while (!digits.empty() && digits.back() == 0) {
                digits.pop_back();
            }
            size = digits.size();
            if (size == 0) {
                sign = 1;
            }
            return *this;
        }
        int add = 0;
        for (size_t i = 0; i < digits.size(); ++i) {
            long long digit = 1LL * big.digits[i] + little.digits[i] + add;
            digits[i] = int(digit % BASE);
            add = int(digit / BASE);
        }
        if (add > 0) {
            digits.push_back(add);
        }
        size = digits.size();
        return *this;
    }

    BigInteger &operator-=(const BigInteger &another) {
        BigInteger copy = another;
        copy.sign *= -1;
        *this += copy;
        return *this;
    }

    BigInteger &operator*=(const BigInteger &another) {
        BigInteger result;
        result.digits.resize(size + another.size, 0);
        for (size_t i = 0; i < size; ++i) {
            int add = 0;
            for (size_t j = 0; j < another.size || add; ++j) {
                long long digitFromAnother = j < another.size ? another.digits[j] : 0;
                long long digit = result.digits[i + j] + 1LL * digits[i] * digitFromAnother + add;
                result.digits[i + j] = int(digit % BASE);
                add = int(digit / BASE);
            }
        }
        while (!result.digits.empty() && result.digits.back() == 0) {
            result.digits.pop_back();
        }
        result.size = result.digits.size();
        result.sign = sign * another.sign;
        if (result.size == 0) {
            result.sign = 1;
        }
        *this = result;
        return *this;
    }

    BigInteger &operator/=(const BigInteger &another) {
        BigInteger result = 0;
        BigInteger residue = 0;
        BigInteger absAnother = another;
        absAnother.sign = 1;
        for (size_t i = size - 1; i + 1 != 0; --i) {
            residue = residue * BASE;
            residue += digits[i];
            if (residue >= absAnother) {
                int leftQuotient = 0;
                int rightQuotient = BASE;
                while (rightQuotient - leftQuotient > 1) {
                    int quotient = leftQuotient + (rightQuotient - leftQuotient) / 2;
                    absAnother * quotient <= residue ? leftQuotient = quotient : rightQuotient = quotient;
                }
                result.digits.push_back(leftQuotient);
                residue -= absAnother * leftQuotient;
            } else if (!result.digits.empty()) {
                result.digits.push_back(0);
            }
        }
        result.size = result.digits.size();
        for (size_t i = 0; i < result.size / 2; ++i) {
            std::swap(result.digits[i], result.digits[result.size - i - 1]);
        }
        result.sign = sign * another.sign;
        if (result.size == 0) {
            result.sign = 1;
        }
        *this = result;
        return *this;
    }

    BigInteger &operator%=(const BigInteger &another) {
        BigInteger abs = *this;
        BigInteger absAnother = another;
        abs.sign = 1;
        absAnother.sign = 1;
        BigInteger result = abs - (abs / absAnother) * absAnother;
        result.sign = sign;
        if (result.size == 0) {
            result.sign = 1;
        }
        *this = result;
        return *this;
    }

    BigInteger& operator--() {
        if (*this == 0 || *this == 1) {
            *this -= 1;
            return *this;
        }
        size_t index = 0;
        while (index < size && ((sign == -1 && digits[index] == 9) || (sign == 1 && digits[index] == 0))) {
            digits[index] = sign == -1 ? 0 : 9;
            ++index;
        }
        if (index == size) {
            digits.push_back(1);
        } else {
            digits[index] += 1 - 2 * (sign == 1);
            if (digits[index] == 0 && index == size - 1) {
                digits.pop_back();
                --size;
            }
        }
        return *this;
    }

    BigInteger operator--(int) {
        *this -= 1;
        return *this;
    }

    BigInteger& operator++() {
        if (*this == 0 || *this == -1) {
            *this += 1;
            return *this;
        }
        size_t index = 0;
        while (index < size && ((sign == -1 && digits[index] == 0) || (sign == 1 && digits[index] == 9))) {
            digits[index] = sign == -1 ? 9 : 0;
            ++index;
        }
        if (index == size) {
            digits.push_back(1);
        } else {
            digits[index] += 1 - 2 * (sign == -1);
            if (digits[index] == 0 && index == size - 1) {
                digits.pop_back();
                --size;
            }
        }
        return *this;
    }

    BigInteger operator++(int) {
        *this += 1;
        return *this;
    }

    explicit operator bool() const {
        return size > 0;
    }

    BigInteger operator*=(const int &multiplier) {
        int add = 0;
        for (size_t i = 0; i < size; ++i) {
            long long cur = 1LL * digits[i] * multiplier + add;
            digits[i] = int(cur % BASE);
            add = int(cur / BASE);
        }
        digits.push_back(add);
        while (!digits.empty() && digits.back() == 0) {
            digits.pop_back();
        }
        size = digits.size();
        return *this;
    }

    BigInteger &operator/=(const int &divisor) {
        int residue = 0;
        for (size_t i = size - 1; i + 1 != 0; --i) {
            long long cur = digits[i] + 1LL * residue * BASE;
            digits[i] = int(cur / divisor);
            residue = int(cur % divisor);
        }
        while (!digits.empty() && digits.back() == 0) {
            digits.pop_back();
        }
        size = digits.size();
        return *this;
    }

    friend std::istream &operator>>(std::istream &in, BigInteger &number);
    friend std::ostream &operator<<(std::ostream &out, const BigInteger &number);
    friend bool operator==(const BigInteger &first, const BigInteger &second);
    friend bool operator<(const BigInteger &first, const BigInteger &second);
    friend class Rational;
    friend bool operator<(const Rational &first, const Rational &second);
    friend BigInteger GCD(BigInteger first, BigInteger second);
};

BigInteger operator+(const BigInteger &first, const BigInteger &second) {
    BigInteger copy = first;
    copy += second;
    return copy;
}

BigInteger operator-(const BigInteger &first, const BigInteger &second) {
    BigInteger copy = first;
    copy -= second;
    return copy;
}

BigInteger operator*(const BigInteger &first, const BigInteger &second) {
    BigInteger copy = first;
    copy *= second;
    return copy;
}

BigInteger operator/(const BigInteger &first, const BigInteger &second) {
    BigInteger copy = first;
    copy /= second;
    return copy;
}

BigInteger operator%(const BigInteger &first, const BigInteger &second) {
    BigInteger copy = first;
    copy %= second;
    return copy;
}

bool operator==(const BigInteger &first, const BigInteger &second) {
    if (first.sign != second.sign || first.size != second.size) {
        return false;
    }
    for (size_t i = 0; i < first.size; ++i) {
        if (first.digits[i] != second.digits[i]) {
            return false;
        }
    }
    return true;
}

bool operator!=(const BigInteger &first, const BigInteger &second) {
    return !(first == second);
}

bool operator<(const BigInteger &first, const BigInteger &second) {
    if (first.sign != second.sign) {
        return first.sign < second.sign;
    }
    if (first.size != second.size) {
        return (first.size < second.size && first.sign == 1) ||
        (first.size > second.size && first.sign == -1);
    }
    for (size_t i = first.size - 1; i + 1 != 0; --i) {
        if (first.digits[i] != second.digits[i]) {
            return (first.digits[i] < second.digits[i] && first.sign == 1) ||
            (first.digits[i] > second.digits[i] && first.sign == -1);
        }
    }
    return false;
}

bool operator>(const BigInteger &first, const BigInteger &second) {
    return second < first;
}

bool operator<=(const BigInteger &first, const BigInteger &second) {
    return !(first > second);
}

bool operator>=(const BigInteger &first, const BigInteger &second) {
    return !(first < second);
}

BigInteger GCD(BigInteger first, BigInteger second) {
    BigInteger powOf2 = 1;
    while (second && first) {
        bool isDividedFirst = first.digits[0] % 2 == 0;
        bool isDividedSecond = second.digits[0] % 2 == 0;
        if (isDividedFirst && isDividedSecond) {
            first /= 2;
            second /= 2;
            powOf2 *= 2;
        } else if (isDividedFirst) {
            first /= 2;
        } else if (isDividedSecond) {
            second /= 2;
        } else {
            first > second ? first -= second : second -= first;
        }
    }
    return first ? first * powOf2 : second * powOf2;
}

std::istream &operator>>(std::istream &in, BigInteger &number) {
    std::string input;
    in >> input;
    number.sign = 1 - 2 * (input[0] == '-');
    number.digits.clear();
    for (size_t i = input.size() - 1; i - (number.sign == -1) < input.size(); i -= BASE_CNT) {
        int digit = 0;
        size_t begin = i + 1 >= BASE_CNT ? i + 1 - BASE_CNT : 0;
        begin = std::max(begin, size_t(number.sign == -1));
        for (size_t j = begin; j <= i; ++j) {
            digit = digit * 10 + (input[j] - '0');
        }
        number.digits.push_back(int(digit));
    }
    number.size = number.digits.size();
    return in;
}

std::ostream &operator<<(std::ostream &out, const BigInteger &number) {
    //x0 + x1 * base + x2 * base^2 + x3 * base^3 + ...
    if (number.size == 0) {
        out << 0;
        return out;
    }
    if (number.sign == -1) {
        out << '-';
    }
    for (size_t i = number.size - 1; i + 1 != 0; --i) {
        int digit = number.digits[i];
        int cnt = 0;
        while (digit > 0) {
            digit /= 10;
            ++cnt;
        }
        if (cnt == 0) {
            cnt = 1;
        }
        if (i != number.size - 1) {
            for (int j = cnt; j < BASE_CNT; ++j) {
                out << '0';
            }
        }
        out << number.digits[i];
    }
    return out;
}
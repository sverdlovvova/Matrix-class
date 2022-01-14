#include <iostream>
#include "biginteger.h"

class Rational;
Rational operator+(const Rational &first, const Rational &second);
Rational operator-(const Rational &first, const Rational &second);
Rational operator*(const Rational &first, const Rational &second);
Rational operator/(const Rational &first, const Rational &second);
bool operator==(const Rational &first, const Rational &second);
bool operator!=(const Rational &first, const Rational &second);
bool operator>(const Rational &first, const Rational &second);
bool operator<=(const Rational &first, const Rational &second);
bool operator>=(const Rational &first, const Rational &second);

class Rational {
private:
    BigInteger numerator = 0;
    BigInteger denominator = 0;

public:
    Rational() = default;

    Rational(const int &number) : numerator(number), denominator(1) {}

    Rational(const BigInteger &number) : numerator(number), denominator(1) {}

    BigInteger getNumerator() const {
        return numerator;
    }

    BigInteger getDenominator() const {
        return denominator;
    }

    Rational operator-() const {
        Rational copy = *this;
        if (numerator.size > 0) {
            copy.numerator.sign *= -1;
        }
        return copy;
    }

    Rational &operator+=(const Rational &another) {
        numerator = (numerator * another.denominator + denominator * another.numerator);
        denominator = denominator * another.denominator;
        BigInteger copyN = numerator;
        BigInteger copyM = denominator;
        copyN.sign = 1;
        copyM.sign = 1;
        BigInteger gcd = (copyN > copyM ? GCD(copyN, copyM) : GCD(copyM, copyN));
        gcd.sign = 1;
        numerator /= gcd;
        denominator /= gcd;
        return *this;
    }

    Rational &operator-=(const Rational &another) {
        Rational copy = another;
        if (another.numerator.size > 0) copy.numerator.sign *= -1;
        *this += copy;
        return *this;
    }

    Rational &operator*=(const Rational &another) {
        numerator *= another.numerator;
        denominator *= another.denominator;
        BigInteger copyN = numerator;
        BigInteger copyM = denominator;
        copyN.sign = 1;
        copyM.sign = 1;
        BigInteger gcd = (copyN > copyM ? GCD(copyN, copyM) : GCD(copyM, copyN));
        gcd.sign = 1;
        numerator /= gcd;
        denominator /= gcd;
        return *this;
    }

    Rational &operator/=(const Rational &another) {
        Rational copy = another;
        std::swap(copy.numerator, copy.denominator);
        if (copy.denominator.sign == -1) {
            copy.numerator.sign = -1;
            copy.denominator.sign = 1;
        }
        *this *= copy;
        return *this;
    }

    std::string toString() const {
        std::string result = numerator.toString();
        if (denominator != 1) {
            result += '/';
            result += denominator.toString();
        }
        return result;
    }

    std::string asDecimal(size_t precision = 0) const {
        std::string result;
        if (numerator == 0) {
            result += '0';
            if (precision > 0) {
                result += '.';
                for (size_t i = 0; i < precision; ++i) {
                    result += '0';
                }
            }
            return result;
        }
        Rational z = numerator / denominator;
        if (z == 0 && numerator.sign == -1) {
            result += '-';
        }
        result += z.toString();
        if (precision > 0) {
            result += '.';
        } else {
            return result;
        }
        size_t pred = result.size();
        BigInteger now = numerator % denominator;
        now.sign = 1;
        BigInteger residue = 0;
        for (size_t i = 0; i < (precision + BASE_CNT - 1) / BASE_CNT; ++i) {
            now = now * BASE;
            if (now >= denominator) {
                int leftQuotient = 0;
                int rightQuotient = BASE;
                while (rightQuotient - leftQuotient > 1) {
                    int quotient = leftQuotient + (rightQuotient - leftQuotient) / 2;
                    denominator * quotient <= now ? leftQuotient = quotient : rightQuotient = quotient;
                }
                residue.digits.push_back(leftQuotient);
                now -= denominator * leftQuotient;
            } else {
                residue.digits.push_back(0);
            }
        }
        residue.size = residue.digits.size();
        int tmp = residue.digits[0];
        size_t cnt = 0;
        while (tmp > 0) {
            tmp /= 10;
            cnt++;
        }
        for (size_t i = 0; i < BASE_CNT - cnt; ++i) {
            result += '0';
        }
        for (size_t i = 0; i < residue.size / 2; ++i) {
            std::swap(residue.digits[i], residue.digits[residue.size - i - 1]);
        }
        residue.sign = 1;
        result += residue.toString();
        while (result.size() > pred + precision) {
            result.pop_back();
        }
        return result;
    }

    explicit operator double() const {
        size_t precision = 100;
        std::string result = (*this).asDecimal(precision);
        double ans = 0;
        for (char digit : result) {
            if (digit != '-' && digit != '.') {
                ans = ans * 10 + (digit - '0');
            }
        }
        for (size_t i = 0; i < precision; ++i) {
            ans /= 10;
        }
        if (numerator.sign == -1) {
            ans *= -1;
        }
        return ans;
    }
};

Rational operator+(const Rational &first, const Rational &second) {
    Rational copy = first;
    copy += second;
    return copy;
}

Rational operator-(const Rational &first, const Rational &second) {
    Rational copy = first;
    copy -= second;
    return copy;
}

Rational operator*(const Rational &first, const Rational &second) {
    Rational copy = first;
    copy *= second;
    return copy;
}

Rational operator/(const Rational &first, const Rational &second) {
    Rational copy = first;
    copy /= second;
    return copy;
}

bool operator==(const Rational &first, const Rational &second) {
    return first.getNumerator() == second.getNumerator() && first.getDenominator() == second.getDenominator();
}

bool operator!=(const Rational &first, const Rational &second) {
    return !(first == second);
}

bool operator<(const Rational &first, const Rational &second) {
    Rational difference = first;
    difference -= second;
    return difference.getNumerator() < 0;
}

bool operator>(const Rational &first, const Rational &second) {
    return (second < first);
}

bool operator<=(const Rational &first, const Rational &second) {
    return !(first > second);
}

bool operator>=(const Rational &first, const Rational &second) {
    return !(first < second);
}

std::istream &operator>>(std::istream &in, Rational &number) {
    BigInteger value;
    in >> value;
    number = Rational(value);
    return in;
}
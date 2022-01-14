int binPow(int number, size_t pow, size_t MOD) {
    if (pow == 0) {
        return 1;
    }
    if (pow & 1) {
        return int((1LL * number * binPow(number, pow - 1, MOD)) % MOD);
    } else {
        long long tmp = binPow(number, pow / 2, MOD);
        return int((tmp * tmp) % MOD);
    }
}

template<size_t V, size_t LeftBound, size_t RightBound>
struct SqrtHelper {
    static const size_t AVG = (LeftBound + RightBound) / 2;
    static const bool equal = (AVG * AVG <= V && (AVG + 1) * (AVG + 1) > V);
    static const bool side = (AVG * AVG < V);
    static const size_t value = SqrtHelper<V, equal ? AVG : side ? AVG + 1 : LeftBound, equal ? AVG : side ? RightBound : AVG - 1>::value;
};

template<size_t V, size_t LeftBound>
struct SqrtHelper<V, LeftBound, LeftBound> {
    static const size_t value = LeftBound;
};

template<size_t N>
struct Sqrt {
    static const size_t value = SqrtHelper<N, 0, N>::value;
};

template<size_t N, size_t Divisor>
struct IsPrimeHelper {
    static const bool value = (N % Divisor != 0 && IsPrimeHelper<N, Divisor - 1>::value);
};

template<size_t N>
struct IsPrimeHelper<N, 1> {
    static const bool value = true;
};

template<size_t N>
struct IsPrime {
    static const bool value = IsPrimeHelper<N, Sqrt<N>::value>::value;
};

template<size_t N>
class Residue {
private:
    int value = 0;

public:
    explicit Residue(int value_) : value(value_ % static_cast<int>(N)) {
        if (value < 0) {
            value += N;
        }
    }

    bool operator==(const Residue &another) const {
        return value == another.value;
    }

    bool operator!=(const Residue &another) const {
        return value != another.value;
    }

    Residue operator+(const Residue &another) const {
        return Residue((value + another.value) % N);
    }

    Residue operator-(const Residue &another) const {
        return Residue((value - another.value + N) % N);
    }

    Residue operator*(const Residue &another) const {
        return Residue((1LL * value * another.value) % N);
    }

    Residue operator/(const Residue &another) const {
        static_assert(IsPrime<N>::value);
        return Residue((1LL * value * binPow(another.value, N - 2, N)) % N);
    }

    Residue &operator+=(const Residue &another) {
        *this = *this + another;
        return *this;
    }

    Residue &operator-=(const Residue &another) {
        *this = *this - another;
        return *this;
    }

    Residue &operator*=(const Residue &another) {
        *this = *this * another;
        return *this;
    }

    Residue &operator/=(const Residue &another) {
        *this = *this / another;
        return *this;
    }

    explicit operator int() const {
        return value;
    }

};
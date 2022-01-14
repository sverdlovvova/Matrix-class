# Matrix class
Implementation of the class Matrix over the Residue and Rational fields.

### Residue class
The Residue<size_t N> class supports arithmetic operations (division is only for simple N, for composite N it gives a compilation error), a constructor from int and explicit conversions to int and back.

### BigInteger class
The BigInteger class supports long integers. The following operations have been implemented:
- Binary operators +, -, *, /, %
- Operators +=, -=, *=, /=, %=
- Unary minus, prefix and postfix increment and decrement
- Comparison operators <=, >=, <, >, ==, !=
- Input from a stream and output to a stream
- toString() method returning a string representation of a number
- Construction from int
- Convert to bool in conditional expressions.

### Rational class
The Rational class is implemented using the BigInteger class. The following operations have been implemented:
- Constructor from BigInteger and int
- Binary operators +, -, *, /
- Operators +=, -=, *=, /=
- Unary minus
- Comparison operators <=, >=, <, >, ==, !=
- Input from a stream and output to a stream
- toString() method that returns a string representation of a number of the form [minus]numerator/denominator
- The asDecimal(size_t precision = 0) method, which returns a number representation as a decimal fraction with precision decimal places
- Cast operator to double

### Matrix class
The class Matrix<size_t N, size_t M> implemented over the Residue and Rational fields. The following operations are supported:
- Default constructor that creates an identity matrix
- Default constructor that creates a matrix from vector<vector<int>> and from initializer_list<initializer_list<int>>
- Binary operators +, -, * (compilation error for matrices of inappropriate sizes)
- Comparison operators ==, !=
- Operators +=, -=
- Multiplication by a number
- *= operator for square matrices
- det() method returning matrix determinant
- transposed() method returning the transposed matrix
- rank() method, which returns the rank of a matrix.
- trace() method, which returns the trace of a matrix.
- inverted() method, which returns an inverse matrix.
- invert() method that inverts the given matrix.
- getRow(unsigned) and getColumn(unsigned) methods return the row and column of the matrix.
- [][] operator can be applied twice to a matrix.
- Square matrices can be declared with one template parameter SquareMatrix<size_t>

import numpy as np
import gmpy2
from math import lcm
import sympy as sp 
from sympy import primerange, symbols, diff, Poly, GF, gcd, factorint, legendre_symbol, Rational, floor, Matrix, div
from itertools import product
from multiprocessing import Pool

# Define x for polynomial operations
x = symbols('x')

# Define parameters
n = 45901
print(f"n = {n}")
d= 3
print(f"d = {d}")


# Function to compute m and generate polynomial f(x)
def polynomial_representation(n, d):
    m = 31 # Calculate the base m
    print(f"m = {m}")
    coefficients = []
    temp = n
    while temp > 0:
        coefficients.append(temp % m)
        temp //= m
    coefficients.reverse()

    # Construct polynomial
    polynomial = ' + '.join(
        f"{'' if coeff == 1 and power > 0 else coeff}{'x' if power > 0 else ''}{f'^{power}' if power > 1 else ''}"
        for power, coeff in zip(range(len(coefficients) - 1, -1, -1), coefficients) if coeff != 0
    )
    print(f"f(x) = {polynomial}")
    return m, coefficients

# Generate polynomial
m, coefficients = polynomial_representation(n, d)

# Polynomial and its derivative
f = sum(coeff * x**power for power, coeff in enumerate(reversed(coefficients)))
f_derivative = diff(f, x)

# Define bounds
B_1 = 42
B_2 = 97
k = 5 # Number of quadratic primes

# Rational base: primes up to B_1 and -1
rational_base = [-1] + list(primerange(2, B_1 + 1))
print(f"Rational Base (Primes ≤ {B_1}, including -1): {rational_base}")
print(f"Cardinality of Rational Base: {len(rational_base)}")

# Algebraic base: (p, r) where r is a root of f(x) mod p for primes ≤ B_2
algebraic_base = []
for p in primerange(2, B_2 + 1):
    roots = [r for r in range(p) if Poly(f, x, domain=GF(p)).eval(r) == 0]
    algebraic_base.extend((p, r) for r in roots)
print(f"Algebraic Base (Primes ≤ {B_2}): {algebraic_base}")
print(f"Cardinality of Algebraic Base: {len(algebraic_base)}")

# Quadratic characters: (q, s) for first k primes > B_2
quadratic_characters = set()
quadratic_primes = list(primerange(B_2 + 1, B_2 + 100))[:k]
for q in quadratic_primes:
    for s in range(q):
        f_mod = Poly(f, x, domain=GF(q)).eval(s)
        f_prime_mod = Poly(f_derivative, x, domain=GF(q)).eval(s)
        if f_mod == 0:
            if f_prime_mod != 0 or (f_prime_mod == 0 and (q, s) not in quadratic_characters):
                quadratic_characters.add((q, s))
print(f"Quadratic Characters (First {k} Primes > {B_2}): {sorted(quadratic_characters)}")
print(f"Cardinality of Quadratic Characters: {len(quadratic_characters)}")

# Define bounds for a and b
a_lower_bound = -80
a_upper_bound = 80
b_upper_bound = 10

# Function to find pairs (a, b) satisfying the conditions
def find_valid_pairs(algebraic_base, rational_base, f, m, d):
    valid_pairs = []
    for b in range(1, b_upper_bound + 1):  # b > 0
        for a in range(a_lower_bound, a_upper_bound + 1):  # a != 0
            if gcd(a, b) != 1 or a == 0:  # Ensure a and b are relatively prime
                continue

            # Compute N(a - bα) = b^d * f(a/b)
            # Ensure a/b is a rational number
            a_b_rational = Rational(a, b)
            N = b**d * Poly(f, x).eval(a_b_rational)

            # Check if N is an integer by using floor and comparing
            if floor(N) != N:  # This ensures that N is not very close to an integer
                raise ValueError(f"Computed N ({N}) is not an integer.")
            
            N = int(N)  # Convert N to an integer

            # Condition 1: All factors of N are in algebraic_base
            factors_N = factorint(N)
            if all((p, r) in algebraic_base for p in factors_N for r in range(p) if Poly(f, x, domain=GF(p)).eval(r) == 0):
                # Condition 2: All factors of a - b * m are in rational_base
                factors_a_b_m = factorint(a - b * m)
                if all(factor in rational_base for factor in factors_a_b_m):
                    valid_pairs.append((a, b))
    return valid_pairs

# Find valid pairs (a, b) satisfying the conditions
valid_pairs = find_valid_pairs(algebraic_base, rational_base, f, m, d)

# New functionality: Enumerate results, print factors and conditions
for index, (a, b) in enumerate(valid_pairs, start=1):
    print(f"\nPair #{index}: (a, b) = ({a}, {b})")

    # Compute and print a - bm and its factors
    a_minus_bm = a - b * m
    print(f"a - bm = {a_minus_bm}, Factors: {factorint(a_minus_bm)}")

    # Compute and print N(a - bα) and its factors
    # N = abs(b**d * Poly(f, x).eval(a / b))      #Old way of calculating N

    # Compute N(a - bα) = b^d * f(a/b)           
    N = b**d * f.subs(x, Rational(a, b))      # Added by moiz           

    # N = int(N)  # Convert N to an exact integer                                 # Added by moiz
    print(f"N(a - bα) = {N}, Factors: {factorint(int(N))}")

    # Check conditions for each factor of N
    for p, _ in factorint(N).items():
        for r in range(p):
            if Poly(f, x, domain=GF(p)).eval(r) == 0:
                if a % p == (b * r) % p:
                    print(f"(p, r) = ({p}, {r}) satisfies a ≡ br (mod p)")

    # Calculate and print the Legendre symbol for (a, b) and (q, s) pairs
    for q, s in quadratic_characters:
        legendre_val = legendre_symbol(a - b * s, q)
        if legendre_val == 1:
            print(f"Legendre Symbol ({a} - {b} * {s} / {q}) = 0")
        elif legendre_val == -1:
            print(f"Legendre Symbol ({a} - {b} * {s} / {q}) = 1")

# Output the results
print(f"\nValid Pairs (a, b): {valid_pairs}")
print(f"Number of Valid Pairs: {len(valid_pairs)}")

# Compute the mod 2 matrix
def compute_mod2_matrix_with_labels(valid_pairs, rational_base, algebraic_base, quadratic_characters, m, f, d):
    matrix = []
    column_labels = []

    # Define column labels
    column_labels.append("Sign(a - b * m)")  # First column label for the sign of a - b * m

    # Add labels for the rational base (p > 1)
    for p in rational_base:
        if p > 1:  # Exclude -1 for labeling
            column_labels.append(f"Rational Base: {p}")

    # Add labels for the algebraic base
    for p, r in algebraic_base:
        column_labels.append(f"Algebraic Base: (p={p}, r={r})")

    # Add labels for the quadratic characters
    for q, s in quadratic_characters:
        column_labels.append(f"Quadratic Character: (q={q}, s={s})")

    # Compute each row of the matrix
    for a, b in valid_pairs:
        # Initialize a row with zeros
        row = []

        # Determine the first element based on a - b * m
        a_minus_bm = a - b * m
        row.append(1 if a_minus_bm < 0 else 0)

        # Add the exponents in the prime factorization of |a - b * m|
        abs_a_minus_bm = abs(a_minus_bm)
        factors_a_minus_bm = factorint(abs_a_minus_bm)
        for p in rational_base:
            if p > 1:  # Only consider primes > 1
                row.append(factors_a_minus_bm.get(p, 0) % 2)  # Exponent mod 2

        # Compute |N(a - bα)| = b^d * f(a/b)
        N = abs(b**d * f.subs(x, Rational(a, b)))
        factors_N = factorint(N)

        # Add the exponents in the prime factorization of N for (p, r) in the algebraic base
        for p, r in algebraic_base:
            if p in factors_N and Poly(f, x, domain=GF(p)).eval(r) == 0 and a % p == (b * r) % p:
                row.append(factors_N[p] % 2)  # Exponent mod 2
            else:
                row.append(0)

        # Add 0 or 1 for each quadratic character based on the Legendre symbol
        for q, s in quadratic_characters:
            legendre_val = legendre_symbol(a - b * s, q)
            row.append(1 if legendre_val == -1 else 0)

        # Append the row to the matrix
        matrix.append(row)

    return np.array(matrix), column_labels

# Generate the matrix and column labels
mod2_matrix, column_labels = compute_mod2_matrix_with_labels(valid_pairs, rational_base, algebraic_base, quadratic_characters, m, f, d)

# Print the column labels
print("\nColumn Labels:")
for i, label in enumerate(column_labels, start=1):
    print(f"{i}: {label}")

# Print the matrix
print("\nMod 2 Matrix:")
for row in mod2_matrix:
    print(f"[{',' .join(map(str, row))}]")

# Print the number of rows and columns
num_rows = mod2_matrix.shape[0]  # Number of rows
num_columns = mod2_matrix.shape[1]  # Number of columns
print(f"\nNumber of Rows: {num_rows}")
print(f"Number of Columns: {num_columns}")



def left_kernel(matrix):
    # Convert input to a sympy Matrix for exact arithmetic
    M = Matrix(mod2_matrix).transpose()
    
    # Compute the null space (left kernel)
    null_space = M.nullspace()

    if not null_space:
        print("No non-trivial left kernel (null space is trivial).")
        return

    print("Left Kernel Vectors:")
    for vec in null_space:
        print([f"{frac.p}/{frac.q}" for frac in vec])
    for vec in null_space:
        denominators = [frac.q for frac in vec]
        lcm_den = lcm(*denominators)  # Unpack the list for lcm

        # Multiply each element by the LCM of denominators
        scaled_vec = [int(frac * lcm_den) for frac in vec]

        # Reduce modulo 2
        mod2_vec = [val % 2 for val in scaled_vec]

        print("\nMod 2 Vector:", mod2_vec)
        product = 1
        for i, val in enumerate(mod2_vec):
            if val == 1:
                a, b = valid_pairs[i]
                term = a - b * m  # Using M's symbol
                product *= term
                print(f"Entry {i+1}: (a, b) = ({a}, {b}), term = {term}")
        f_prime = diff(f, x)
        f_prime_m = f_prime.subs(x, m)

        final_product = product * (f_prime_m**2)
        sqrt_product = sp.sqrt(final_product)
        u = sqrt_product % n


        print("Product of (a - b*m):", product)
        print("Product multiplied by f'(m)^2:", final_product)
        print("Square root in Z:", sqrt_product)
        print("u (square root mod n):", u)
        
        product = Poly(1, x)
        for i, val in enumerate(mod2_vec):
            if val == 1:
                a, b = valid_pairs[i]
                term = Poly(a - b * x, x)
                product = div(product * term, f, domain='ZZ')[1]  # Reduce mod f
                print(f"Entry {i+1}: (a, b) = ({a}, {b}), term = {term.as_expr()} (mod f)")

        f_prime = diff(f, x)
        final_product = div(product * Poly(f_prime**2, x), f, domain='ZZ')[1]

        print("Product of (a - b*x) mod f:", product.as_expr())
        print("f'(x)^2:", f_prime**2)
        print("Final Product * f'(x)^2 mod f:", final_product.as_expr())

left_kernel(Matrix)

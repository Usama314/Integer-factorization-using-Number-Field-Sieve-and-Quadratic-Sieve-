import math
from math import lcm
from sympy import primerange, legendre_symbol, factorint, Matrix, gcd
def quadratic_sieve(n, N):
    # Step 1: Calculate B
    B = math.exp(math.sqrt(0.5 * math.log(n) * math.log(math.log(n))))
    B = math.floor(B)
    print(f"B: {B}") # Print the value of B
    
    # Step 2: Find all primes up to floor(B+1)
    primes = list(primerange(2, B + 1))
    
    # Step 3: Compute the Legendre symbol (n/p) for all primes > 2
    factor_base = [-1, 2]  # Start the factor base with -1 and 2
    for p in primes:
        if p > 2 and legendre_symbol(n, p) == 1:
            factor_base.append(p)
    
    # Step 4: Print the factor base and its cardinality
    print("Factor Base:", factor_base)
    print("Cardinality of Factor Base:", len(factor_base))
    
    # Step 5: Define the sieving interval using floor(sqrt(n))
    x_start = math.ceil(math.sqrt(n)) - N  # Starting value of x
    x_end = math.ceil(math.sqrt(n)) + N    # Ending value of x
    print(f"Sieving Interval: [{x_start}, {x_end}]")
    
    # Step 6: Find B-smooth numbers using x^2 - n
    b_smooth_numbers = []
    exponent_vectors = []  # To store exponent vectors mod 2
    x_values = []  # Store x values corresponding to smooth numbers
    for x in range(x_start, x_end + 1):
        value = x**2 - n
        factorization = factorint(abs(value))  # Factorize x^2 - n
        if value < 0:  # If x^2 - n is negative, include -1 in the factorization
            factorization[-1] = 1
        if all(p in factor_base for p in factorization):  # Check if all factors are in the factor base
            # Create the exponent vector mod 2
            vector = []
            for p in factor_base:
                exponent = factorization.get(p, 0)
                vector.append(exponent % 2)
            b_smooth_numbers.append((x, value, factorization))
            exponent_vectors.append(vector)
            x_values.append(x)
    
    # Step 7: Print x, x^2 - n, and its factorization
    print("\nB-smooth Numbers:")
    for x, value, factorization in b_smooth_numbers:
        print(f"x: {x}, x^2 - n: {value}, Factorization: {factorization}")
    
   # Step 8: Print the exponent matrix mod 2 with row labels for x values
    print("\nExponent Matrix Mod 2 :")
    for vector in exponent_vectors:
        print(f"{vector}")
        
        
       
        
        
        
   
    # Step 9: Linear Algebra - Find the kernel of the transpose
    def find_quadratic_relations(factorizations):
        """ Solve the system using linear algebra mod 2 and return null space vectors. """
        mat = Matrix(factorizations).transpose()
        null_space = mat.nullspace()
        if not null_space:
            print("No non-trivial left kernel (null space is trivial).")
            return []
        
        print("Left Kernel Vectors:")
        null_vectors = []
        for vec in null_space:
            print([f"{frac.p}/{frac.q}" for frac in vec])
            denominators = [frac.q for frac in vec]
            lcm_den = lcm(*denominators)  # Unpack the list for lcm

            # Multiply each element by the LCM of denominators
            scaled_vec = [int(frac * lcm_den) for frac in vec]

            # Reduce modulo 2
            reduced_null_space = [val % 2 for val in scaled_vec]
            null_vectors.append(reduced_null_space)
            print("\nMod 2 Vector:", reduced_null_space)
        return null_vectors
    
    null_vectors = find_quadratic_relations(exponent_vectors)
    print("\nQuadratic Relations (Null Vectors Reduced mod 2):")
    for i, nv in enumerate(null_vectors):
        print(f"Relation {i + 1}: {nv}")

    # Step 10: Factorization for all relations
    def find_factors_for_all_relations(smooth_values, null_vectors, n):
        factors = []
        for i, vector in enumerate(null_vectors):
            if not isinstance(vector, list):  # Check if vector is iterable
                print(f"Invalid vector: {vector}")
                continue
            selected_indices = [j for j, bit in enumerate(vector) if bit == 1]
            x_vals = [smooth_values[j][0] for j in selected_indices]
            x = math.prod(x_vals) % n  # x is reduced mod n
            
            # Compute y^2
            y_squared = 1
            y_factorization = {}
            for j in selected_indices:
                factorization = smooth_values[j][2]
                for p, exp in factorization.items():
                    y_factorization[p] = y_factorization.get(p, 0) + exp

            # Calculate y^2 as a product of all factors
            for p, exp in y_factorization.items():
                y_squared *= pow(p, exp)
            y_squared %= n  # Reduce y_squared mod n

            # Replace -1 with 1 in factorization output
            y_factorization_for_sqrt = {1 if p == -1 else p: exp for p, exp in y_factorization.items()}

            # Reduce y_factorization exponents by 2 (square root)
            y_factors_sqrt = {p: exp // 2 for p, exp in y_factorization_for_sqrt.items()}
            y = 1
            for p, exp in y_factors_sqrt.items():
                y *= pow(p, exp)
            y %= n  # y is reduced mod n

            # Print y^2 and the factorization
            factorization_str = " × ".join([f"{p}^{exp}" for p, exp in y_factorization.items()])
            sqrt_str = " × ".join([f"{p}^{exp}" for p, exp in y_factors_sqrt.items()])
            print(f"\nQuadratic Relation {i + 1}:")
            print(f"Selected Indices: {selected_indices}")
            print(f"x values for Relation: {x_vals}")
            print(f"x (mod n) = {x}")
            print(f"y^2 Factorization: {factorization_str}")
            print(f"y^2 = {y_squared} (mod n)")
            print(f"y = {sqrt_str} (mod n) = {y}")

            # Compute GCDs
            gcd1 = gcd(x - y, n)
            gcd2 = gcd(x + y, n)
            print(f"GCD(x - y, n): {gcd1}")
            print(f"GCD(x + y, n): {gcd2}")
            
            factors.append((gcd1, gcd2))
        return factors

    factors = find_factors_for_all_relations(b_smooth_numbers, null_vectors, n)
    print("\nFactors Found:")
    for f1, f2 in factors:
        print(f"{f1}, {f2}")
    if factors:
        return factors
    else:
        print("\nFailed to find factors.")
        return None

# Example Usage
n = 18446744073709551617 # Replace with the composite number you want to factorize
N = 20000  # Range for x values
result = quadratic_sieve(n, N)

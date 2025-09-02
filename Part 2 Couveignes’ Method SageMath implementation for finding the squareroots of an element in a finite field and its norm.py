# SageMath implementation for finding the squareroots of an element in a finite field and its norm
# Define the irreducible polynomial f over QQ
f = x^3 + 16*x^2 + 23*x + 21

# List of primes
primes = [1229, 1279, 1289, 1303]

# Loop over each prime
for p in primes:
    # Define the finite field GF(p^3) using modulus f reduced modulo p
    k.<x> = GF(p^3, modulus=f)
    print("Finite field defined by modulus f over GF(p):")
    print(k)
    
    # Define the element y in the finite field k 
    y = 954522564040241878802*x^2 + 1473437779284753373541 *x + 1465552076529211866174 
    print(f"Element y in GF({p}^3):")
    print(y)
    
    # Compute a square root of y in the finite field (if it exists)
    if y.is_square():
        a = y.sqrt()
        print(f"Square root of y in GF({p}^3): {a}")
        # Compute the second square root (since in a finite field of 0dd characteristic, if y has a sqrt, it has exactly two)
        b =-a
        print("Second square root b of y (i.e.,-a):")
        print(b)
        # Compute the norm of both square roots (norm from the fiel  extension to the base field GF(p))
        print("Norm of a:")
        print(a.norm())
        print("Norm of b:")
        print(b.norm())
    else:
        print(f"No square root exists for y in GF({p}^3).")

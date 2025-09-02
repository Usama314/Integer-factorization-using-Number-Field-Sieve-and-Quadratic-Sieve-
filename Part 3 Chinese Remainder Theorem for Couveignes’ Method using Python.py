#Chinese Remainder Theorem for Couveignes’ Method using Python
from sympy import mod_inverse, symbols, expand

# Step 1: Define the primes (moduli)
p = [
    1229,
    1279,
    1289,
    1303,
    1327
    ]

# Step 2: Polynomial coefficients modulo each prime (for a2 * x^2 + a1 * x + a0)
a2_list = [
    1202,
    149,
    800,
    594,
    564
          ] 
a1_list = [
   837,
  837,
  815,
  1001,
  472
          ]
a0_list = [
    526,
    610,
    1022,
    201,
    269
]

# Combine all coefficient lists
coeff_lists = [a2_list, a1_list, a0_list]

# Step 3: Compute total modulus M
M = 1
for pi in p:
    M *= pi

# Step 4: Compute Mi and yi
Mi_list = [M // pi for pi in p]
yi_list = [mod_inverse(Mi_list[i], p[i]) for i in range(len(p))]

# Step 5: CRT function
def crt(coeffs, Mi, yi, M):
    total = 0
    for i in range(len(coeffs)):
        total += coeffs[i] * Mi[i] * yi[i]
    return total % M
print("\nM =", M)


print("\nMi values:")
for i, Mi in enumerate(Mi_list, 1):
    print(f"M{i} = {Mi}")

print("\nyi values (modular inverses):")
for i, yi in enumerate(yi_list, 1):
    print(f"y{i} = {yi}")
# Step 6: Apply CRT and reduce each coefficient modulo M
crt_results = [crt(coeffs, Mi_list, yi_list, M) % M for coeffs in coeff_lists]
for idx, result in enumerate(crt_results):
    print(f"a{2-idx }_crt =", result)
# Step 7: Print each coefficient result - a_crt - M
for idx, result in enumerate(crt_results):
    print(f"a{2-idx }_crt - M =", result - M)

# Step 8: Build and display the final polynomial
x = symbols('x')
poly = sum(crt_results[i] * x**(2 - i) for i in range(len(crt_results)))
print("\nPolynomial r(x):")
print(expand(poly))
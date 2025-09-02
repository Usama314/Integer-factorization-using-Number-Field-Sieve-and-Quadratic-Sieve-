#Finding a Square Root in a Number Field Using SageMath
n = 45901
m = 31
f(x) = x^3 + 16*x^2 + 23*x + 21
y = 1892 # y is sqrt in Z
K.<x> = NumberField(f(x)) # f(x) is a polynomial

b = 954522564040241878802*x^2 + 1473437779284753373541*x + 1465552076529211866174 # a is gamma
    
# Check if 'b' is a square
if b.is_square():
  a = b.sqrt()
  print("Square Root of b:", a)
 
# Convert a to a polynomial expression
  a_poly = a.polynomial()
 
 # Evaluate at y = m
  a_m = a_poly(m)
  print(f"Value of a at x = {m}:", a_m)
 
# Compute a_m mod n
  x = a_m % n
  print(f"x mod {n}:", x)
  gcd_plus = gcd(x + y, n)
  gcd_minus = gcd(x- y, n)
 
  print(f"gcd(x + y, {n}):", gcd_plus)
  print(f"gcd(x-y, {n}):", gcd_minus)
 
else:
  print("a is not a square in K")
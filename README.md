# Integer-factorization-using-Number-Field-Sieve-and-Quadratic-Sieve-
This repository provides simple implementations of the Number Field Sieve (NFS) and the Quadratic Sieve (QS) for integer factorization. While the codes are not highly optimized, they serve as useful examples for understanding the fundamental principles. NFS is implemented using Python and SageMath, while QS is implemented in Python.
The code for the Number Field Sieve is written in two parts. First, Python is used to handle most of the computations, including finding the polynomial, constructing the factor base, and determining the elements $y \in \mathbb{Z}$ and $\beta \in \mathbb{Z}[\theta]$. Once the element $\beta \in \mathbb{Z}[\theta]$ is found, SageMath is employed to compute its square root in the ring $\mathbb{Z}[\theta]$. Finally, computing $\gcd(x \pm y, n)$ yields the factors of $n$.

#Finding a square root in a number field of odd degree.
Couveignes’ Method is Used to Find a Square Root in a Number Field of Odd Degree
Couveignes’ Method is implemented in three parts: the first part is used to find the norm of an element and its square root using SageMath; the second part finds the square root of an element in a finite field using SageMath; and the third part finds a square root in a number field using the Chinese Remainder Theorem using Python.

README


# Recursive Encoding Algorithm for Speilman Codes in Rust

The purpose of this project was to explore efficient bitwise operations and memory allocation using common structures in RUST. This Rust program demonstrates a recursive encoding algorithm for generating error-correcting Speilman codes. It uses random sparse matrices and bitwise operations for encoding.

The implementation is based on the 1996 paper "Linear-Time Encodable and Decodable Error-Correcting Codes" by Daniel A. Spielman. We define the code for any dimension $2^\ell$ with $\ell \geq \ell_0$. The code has a rate of $\frac{1}{4}$, meaning the encoding of a $2^\ell$-bit message is $2^{\ell+2}$-bits long. 

https://www.cs.yale.edu/homes/spielman/Research/ITsuperc.pdf. 

## Encoding Procedure
The encoding procedure is recursive. Let $E_\ell(m)$ denote the encoding function for length $2^\ell$ messages. The base case $E_{\ell_0}(m)$ is the encoding function of a random linear code with dimension $2^{\ell_0}$, rate $\frac{1}{4}$, and random generator matrix $G_0$. The recursive case uses random sparse binary matrices $A_t \in \{0,1\}^{2^{t-1} \times 2^t}$ with $2^{t-1}$ rows and $2^t$ columns. The $A_t$ matrices are sparse in that each column has exactly $g$ non-zero entries. The positions of the non-zero entries chosen independently at random.

The encoding algorithm $E_\ell(m)$ is as follows:

1. If $\ell$ is $\ell_0$, use the base case encoder defined by $G_0$.
2. Compute a checksum of the message, $x = A_\ell m$. $x$ is $2^{\ell-1}$ bits.
3. Recursively encode $x$ with $y = E_{\ell-1}(x)$. $y$ is $2^{\ell+1}$ bits.
4. Compute $z = A_{\ell+1}y$. $z$ is $2^\ell$ bits.
5. Return the concatenation $m\|y\|z$.

## Selection of Parameters

The program takes three command-line arguments: `l0`, `lmax`, and `g`, which determine the encoding parameters.  The parameters of the code are $g$, the column-density of the sparse matrices, $\ell_0$, the dimension of the base-case encoder, and $\ell_{max}$, the maximum message size. With these parameters, the matrices $G_0$ and $A_t$ are generated ahead of time and stored as parameters of the code for all $t = \ell_0 + 1,..., \ell_{max} + 1$.

This implementation was able to compute parameters $g = 25$, $\ell_0 = 12$, $\ell_{max} = 22$ efficiently on a desktop computer.

# Prerequisites

Before running the program, ensure you have Rust, Cargo and its dependencies installed on your system.

# Usage

To compile and run the program, use the following command:

```bash
cargo run --release <l0> <lmax> <g>
```
# Testing

Simple printing flag const PRINT was implemented and can be used that print to screen for small values of l0, lmax. 

# Non obvious decisions and other considerations

I experimented with a few data structures for the matrices and landed on implementing the sparce $A_t$ matrices as vectors of vectors of u32 elements that represent indices of nonzero matrix elements for each row.

I kept the base case G0 encoding matric as a vector of bitvecs, as $G_0$ is not a sparce matrix and using the above data structures would cause the memory usage to ballon with larger $\ell_0$, for example $\ell_0=15$ required about 12GB of RAM just to generate $G_0$.

To accomodate the two representations, two matrix-vector multiplication functions were implemented to accomodate the different matrix representation inputs, but have the same output type of bitvec.

The bitvec multiplication function for $G_0$ uses standard bitwise operations for dot products.

The indexing multiplication function sets a result bit to 0, then traverses the vector of u32 indices inherently knowing that all other entries are zero and the indexed entries are 1, checks to see if the corresponding bit in the bitvec message is 1 or 0, and flips the result bit everytime there is a 1.

The following design decisions were made given the scope of this project:
    - Basic command line input for running the program.
    - Bitvec crate was used to improve memory usage and operation speed.
    - The implementation of matrices $A_t$ using only indices of nonzero elements gave a huge boost even though u32 elements were used, as I assume $g$ will intentially be logarithmic or even constant when compared to message length $2^l$
    - Always prints elapsed time for $G_0$. generation, all At generation, and the actual encoding to help me (and you guys) with testing.
    - Does not output the encoded message anywhere, as this is memory management excersice, but could add a file output in the future.
    - The cosntruction of the $A_t$ could be sped up using rayon crate where multiple $A_t$ are generated at the same time.






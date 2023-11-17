extern crate rand;   
extern crate bitvec; 

use std::env;                // for command line inputs
use rand::Rng;               // Used to generate randomness for A, G0 matrices and message
use bitvec::prelude::BitVec; // Used for bit message, output encoding, interim vectors in recursion and G0 matrix.

// For debugging help only
const PRINT: bool = false;



////////////////// FUNCTIONS //////////////////

//////////////////
/// Description: Construct a randomized generator matrix G0 so that encoding has a rate of 1/4. Uses 
///              Vec<BitVec> data structure to reduce memory balooning becuase this matrix is not sparse 
///              and can potentially be very large. 
///  
/// Input:       l0: usize       - Dictates the dimension of the matrix. 
/// 
/// Output:      g0: Vec<BitVec> - Matrix with 2^l0 columns, 2^(l0+2) rows as a vector of BitVecs.
///
fn construct_g0_matrix(l0: usize) -> Vec<BitVec> {
    let mut rng = rand::thread_rng();
    let mut g0 = vec![BitVec::with_capacity(1 << l0); 1 << (l0 + 2)]; // allocates contingent memory

    // Pushes random 0 or 1 values to populate the Vec of BitVecs matrix
    for i in 0..(1 << (l0 + 2)) {
        for _j in 0..(1 << l0) {
            g0[i].push(rng.gen_range(0..=1) == 1);
        }
    }

    g0
}

//////////////////
/// Description: Construct a randomized sparse matrix of specified dimensions and column-density. 
///              Only stores non zero value indices for each row as usize elements. Since the
///              output matrix is sparse, the resulting storage size is lxg instead of l^2. Given
///              that g is relatively small to 2^l, we greatly reduce memory and time to compute 
///              matrix/vector multiplications.                  
///  
/// Inputs:      num_rows: usize - Specified number of rows for the output matrix.
///              num_cols: usize - Specified number of columns for the output matrix.
///              g: usize        - Column density of the sparce matrix output.
/// 
/// Output:      matrix: Vec<Vec<usize>> - Matrix representation that only stores non zero value indices 
///                                        in each row with only (2^l)*g usize elements in total.
/// 
fn construct_random_sparse_matrix(num_rows: usize, num_cols: usize, g: usize) -> Vec<Vec<usize>> {
    let mut matrix = vec![Vec::new(); num_rows];
    let mut rng = rand::thread_rng();

    //For every column 
    for col in 0..num_cols {
        //Create list of length g with random unqiue indices for non zero positions
        let mut nonzero_positions = Vec::with_capacity(g);
        while nonzero_positions.len() < g {
            let row = rng.gen_range(0..num_rows);
            if !nonzero_positions.contains(&row) {
                nonzero_positions.push(row);
            }
        }
        //Push the g nonzero positions to the appropriate row
        for &row in &nonzero_positions {
            matrix[row].push(col as usize);
        }
    }

    matrix
}


//////////////////
/// Description: Generate a vec of random sparse matrices with At in {0,1}^(2^t-1 x 2^t) where l0 < t <= lmax + 1                    
///  
/// Inputs:      num_rows: usize - Specified number of rows for the output matrix.
///              num_cols: usize - Specified number of columns for the output matrix.
///              g: usize        - Column density of the sparce matrix output.
/// 
/// Output:      matrix: Vec<Vec<usize>> - Only stores 1 coordinates as a usize elements representing 
///                                        the rows since the output matrix is sparse, resulting 
///                                        in lxg size instead of l^2. This greatly reduces memory and 
///                                        time to compute matrix/vector multiplications given g is small.
///
fn generate_all_sparse_matrices(l0: usize, lmax: usize, g: usize) -> Vec<Vec<Vec<usize>>> {
    let mut matrix_param: Vec<Vec<Vec<usize>>> = Vec::new();
    
    for t in (l0 + 1)..=(lmax + 1) {
        let matrix = construct_random_sparse_matrix(1 << (t - 1), 1 << t, g);
        matrix_param.push(matrix);
    }
    
    matrix_param
}


//////////////////
/// Description: Perform matrix-vector multiplication using bitwise operations on bitvec data structure 
///              representation of a matrix multiplied with bitvec vector. This function is only used in the base 
///              case encoding with G0 matrix.                    
///  
/// Inputs:      matrix: &Vec<BitVec> - pointer to dence encoding matrix using bitvec vectors.
///              vector: &BitVec      - pointer to message to be encoded represented by bitvec vector.
/// Output:      result: BitVec       - the resulting encoding as a bitvec vector.
/// 
/// 
fn matrix_vector_multiply_bitvec(matrix: &Vec<BitVec>, vector: &BitVec) -> BitVec {
    let num_rows = matrix.len();
    let num_cols = matrix[0].len();
    let mut result = BitVec::with_capacity(matrix.len());
    
    for i in 0..num_rows {
        let mut current_result_bit = false;
        for j in 0..num_cols {
            if matrix[i][j] && vector[j] {
                current_result_bit = !current_result_bit;
            }
        }
        result.push(current_result_bit);
    }
    
    result
}


//////////////////
/// Description: Perform matrix-vector multiplication using nonzero indices data structure representation 
///              of a matrix multiplied with bitvec vector. This function is used in all recursive matrix-vector
///              multiplications of the At's.                    
///  
/// Inputs:      matrix: &Vec<Vec<usize>> - pointer to sparce encoding matrix using nonzero index data structure
///              vector: &BitVec          - pointer to message to be encoded represented by bitvec vector.
/// Output:      result: BitVec           - the resulting encoding as a bitvec vector.
/// 
/// 
fn matrix_vector_multiply_indexing(matrix: &Vec<Vec<usize>>, vector: &BitVec) -> BitVec {
    let num_rows = matrix.len();
    let mut result = BitVec::with_capacity(num_rows);
    
    for i in 0..num_rows {
        let mut current_result_bit = false;
        for &coord in &matrix[i]{
            if vector[coord] {
                current_result_bit = !current_result_bit;
            }
        }
        result.push(current_result_bit);
    }
    
    result
}


//////////////////
/// Description: Perform recursive encoding algorithm E_l(m) as follows:
///                 1. If l is l0, use the base case encoder defined by G0.
///                 2. Compute a checksum of the message, x = A_l*m. (x is 2^{l-1} bits)
///                 3. Recursively encode x with y = E_{l−1}(x). (y is 2^{l+1} bits)
///                 4. Compute z = A_{l+1}*y. (z is 2^l bits)
///                 5. Return the concatenation m||y||z. (Total of 2^{l+2} bits) 
///                                
///  
/// Inputs:      message: &BitVec                 - pointer to message to be encoded represented by bitvec vector.
///              g0: &Vec<BitVec>                 - pointer to G0 matrix for base case.
///              matrix_param: &[Vec<Vec<usize>>] - 
///              l: usize                         - power of length of current message
///              l0: usize                        - power of minimal message length for recursion
/// Output:      result: BitVec                   - the resulting encoding of the message.
/// 
/// 
fn recursive_encode(message: &BitVec,
                    g0: &Vec<BitVec>,
                    matrix_param: &[Vec<Vec<usize>>],
                    l: usize,
                    l0: usize,
                    ) -> BitVec {

    //Base case, uses matrix vector multiply specialized for bitvec data structure representation
    if l == l0 {
        return matrix_vector_multiply_bitvec(g0, &message);

    //Recursive case, uses matrix vector multiply specialized for non zero position data structure representation
    } else {
        // Compute checksum x = A_{l} * message
        let x = matrix_vector_multiply_indexing(&matrix_param[l - l0 - 1], &message);

        // Recursive encoding of x as y = E_{l−1}(x)
        let y = recursive_encode(&x, g0, &matrix_param, l - 1, l0);

        // Compute z = A_{l+1} * y
        let z = matrix_vector_multiply_indexing(&matrix_param[l - l0], &y);
        
        // Create the result BitVec with size 2^{l+2}
        let mut result = BitVec::with_capacity(1 << l+2);

        // Extend the result BitVec with the bits from message, y, and z to avoid cloning
        result.extend_from_bitslice(&message);
        result.extend_from_bitslice(&y);
        result.extend_from_bitslice(&z);

        result
    }
}


//MAIN FUNCTION
fn main() {

    //Some basic command line user input checking  
    let args: Vec<String> = env::args().collect();

    if args.len() != 4 {
        eprintln!("Usage: {} <l0> <lmax> <g>", args[0]);
        std::process::exit(1);
    }

    //Assignment of inputs
    let l0: usize = args[1].parse().expect("Invalid l0 value");     // Base case dimension
    let lmax: usize = args[2].parse().expect("Invalid lmax value"); // Column-density of sparse matrices
    let g: usize = args[3].parse().expect("Invalid g value");       // Maximum message size (e.g., 2^lmax bits)

    //
    //Generate random message with length 2^lmax
    let mut rng = rand::thread_rng();
    let mut message = BitVec::with_capacity(1 << lmax);
    for _ in 0..(1 << lmax) {
        message.push(rng.gen::<bool>());
    }
    // for debugging
    if PRINT {
        println!("Message:"); 
        for entry in message.iter() {
            if *entry {
                print!("1");
            } else {
                print!("0");
            }
        }
        println!();
        println!();
    }
    

    // Create random generator matrix G0 so that encoding has rate 1/4 (i.e., 2^l0 columns, 2^(l0+2) rows)
    //Timing to help with testing
    use std::time::Instant;
    let now_g = Instant::now();
    let g0 = construct_g0_matrix(l0);
    let elapsed_g = now_g.elapsed();
    println!();
    println!("G0 Matrix Generated with elapsed: {:.2?}", elapsed_g);
    println!();
    //for debugging
    if PRINT {
        println!("Generated G0 Matrix:");
        for row in &g0 {
            for entry in row.iter() {
                if *entry {
                    print!("1 ");
                } else {
                    print!("0 ");
                }
            }
            println!();
        }
        println!();
        println!();

    }
    

    // Generate matrix parameters and store in matrix_param vector
    //Timing to help with testing
    let now_m = Instant::now();
    let matrix_param = generate_all_sparse_matrices(l0, lmax, g);
    let elapsed_m = now_m.elapsed();
    println!("A Matrices Generated with elapsed: {:.2?}", elapsed_m);
    println!();
    //for debugging
    if PRINT {
        for (t, matrix) in matrix_param.iter().enumerate() {
            println!("Matrix Parameter A({}):", l0 + t + 1);
            for row in matrix {
                for i in 0..(1<<(l0 + t + 1)) {
                    if row.contains(&(i as usize)) {
                        print!("1 ");
                    } else {
                        print!("0 ");
                    }
                }
                println!();
            }
            println!();
        }
        println!();
        println!();
    }
    
    
    //Invoke recursive encoding algorithm
    //Timing to help with testing
    let now = Instant::now();
    let encoding = recursive_encode(&message, &g0, &matrix_param, lmax, l0);
    let elapsed = now.elapsed();
    println!("Encoding finished elapsed: {:.2?}", elapsed);
    println!();

    // for debugging
    if PRINT {
        println!("Encoding: "); 
        for entry in encoding.iter() {
            if *entry {
                print!("1");
            } else {
                print!("0");
            }
        }
        println!();
        println!();
    }

}
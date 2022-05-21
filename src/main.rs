// main.rs
mod text;
use rand_distr::{Normal, Distribution};
use rand::{thread_rng, Rng};
use std::process;
use rand::prelude::SliceRandom;
use std::time::Instant;
use std::time::Duration;


fn main() {

    /*let mut p = 10;
    let mut n = 3;
    let mut p = 121;
    let mut n = 11;
    let mut B = alpha(n as f64);
    let mut E = 0.5;
    let m = ((1.0 + E) * (n as f64 + 1.0) * (p as f64).log2());
    let m = m.round() as u64;
    let mut s = vec![2,7,3];
    let mut a_vector = vec![
        vec![1,0,2],
        vec![8,3,2],
        vec![1,6,3],
        vec![5,5,7],
        vec![2,0,4],
        vec![2,3,7],
        vec![8,5,2],
        vec![1,0,8],
        vec![9,5,7],
        vec![0,2,2],
        vec![1,9,1],
        vec![0,8,2],
        vec![7,9,9],
        vec![9,4,6],
        vec![7,4,1],
        vec![6,0,6],
        vec![2,3,6],
        vec![7,2,6],
        vec![6,3,8],
        vec![8,0,5],
    ];
    let mut e_vector = vec![
        9,1,0,0,1,0,0,9,1,0,0,0,1,0,0,1,0,9,0,1
    ];
    //let mut e_vector = vec![
    //    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    //];
    let mut lst = Vec::new();
    for i in 0..a_vector.len() {
        let mut pair = pair {
            a: a_vector[i].clone(),
            b: (inner_product(&a_vector[i], &s, p) + e_vector[i]).rem_euclid(p)
        };
        lst.push(pair);
    }
    let public_key = public_key {
        m: m,
        n: n,
        p: p,
        list: lst
    };
    let private_key = private_key {
        p: p,
        s: s.clone()
    };*/
    

    let (public_key, private_key, e_vector) = create_system(5, 0.5);
    let p = public_key.p;

    let mut upper = (p/2)/2;
    let mut lower = p - upper;
    println!("acceptable error_sum range 1: [0, {}]", upper);
    println!("acceptable error_sum range 2: [{}, {}]", lower, p-1);

    // test bit 0
    let mut bit1 = 0; 
    let (pair, subset) = encrypt_bit(bit1, &public_key);
    let decrypted_bit = decrypt_bit(&pair, &private_key);
    
    let mut error_sum:u64 = 0;
    for index in &subset {
        error_sum += e_vector[*index as usize];
    }
    error_sum = error_sum.rem_euclid(p);
    println!("error_sum % p = {}", error_sum);
    if bit1 == decrypted_bit {
        println!("bit '0' decrypted correctly");
    } else {
        println!("bit '0' decrypted incorrectly");
    }

    // test bit 1
    let mut bit2 = 1; 
    let (pair, subset) = encrypt_bit(bit2, &public_key);
    let decrypted_bit = decrypt_bit(&pair, &private_key);
    let mut error_sum:u64 = 0;
    for index in &subset {
        error_sum += e_vector[*index as usize];
    }
    error_sum = error_sum.rem_euclid(p);
    println!("error_sum % p = {}", error_sum);
    if bit2 == decrypted_bit {
        println!("bit '1' decrypted correctly");
    } else {
        println!("bit '1' decrypted incorrectly");
    }


    /*let (public_key, private_key) = create_system(5, 0.5);
    let ciphertext = encrypt(String::from("one two three four five"), &public_key);
    let plaintext = decrypt(&ciphertext, &private_key);
    println!("plaintext = {}", plaintext);*/

    // Testing
    /*let messages = all_possible_messages(10);
    let mut iterations = 10;
    let mut av_correct = 0;
    let mut av_incorrect = 0;
    let mut total = 0;
    let (public_key, private_key, e_vector) = create_system(5, 0.5);
    for i in 0..iterations {
        let mut correct = 0;
        let mut incorrect = 0;
        total = 0;
        for message in &messages {
            let ciphertext = encrypt_vector(message, &public_key);
            let plaintext = decrypt_vector(&ciphertext, &private_key);
            if is_equal(message, &plaintext, 10) {
                correct += 1;
            } else {
                incorrect += 1;
            }
            total +=1 ; 
        }  
        println!("iteration {}: correct = {}/{}", i, correct, total);
        println!("iteration {}: incorrect = {}/{}\n", i, incorrect, total);
        av_correct += correct;
        av_incorrect += incorrect;
    }
    println!("average correct over {} tests = {}/{}", iterations, av_correct/iterations, total);
    println!("average incorrect over {} tests = {}/{}", iterations, av_incorrect/iterations, total);*/

}

fn create_system(n:u64, E:f64) -> (public_key, private_key, Vec<u64>) {
    let mut lower = n*n;
    let mut upper = 2*lower;
    let mut p = rand::thread_rng().gen_range(lower..upper) as u64;
    if p < 2 {  
        println!("Error: p must be >= 2");
        process::exit(1);
    }
    let m = ((1.0 + E) * (n as f64 + 1.0) * (p as f64).log2());
    let m = m.round() as u64;
    let mut B = alpha(n as f64);
    //println!("p = {}", p);
    //println!("m = {}", m);
    //let mut p = 4;
    //let mut m = 6;
    //let mut p = 5;
    //let mut m = 3;
    //let mut p = 6;
    //let mut m = 4;
    let private_key = private_key {
        p: p,
        s: draw_vector_from_Znp(n, p)
    };
    let (lst, e_vector) = get_list(m, n, p, B, &private_key.s);
    let public_key = public_key {
        m: m,
        n: n,
        p: p,
        //list: get_list(m, n, p, B, &private_key.s)
        list: lst
    };
    return (public_key, private_key, e_vector);
}

#[derive(Clone)]
pub struct pair {
    pub a: Vec<u64>,
    pub b: u64
}

pub struct private_key {
    pub s: Vec<u64>,
    pub p: u64
}

pub struct public_key {
    pub m: u64,
    pub n: u64,
    pub p: u64,
    pub list: Vec<pair>
}

/*pub fn encrypt(message:String, public_key:&public_key) -> Vec<pair> {
    let mut plaintext = string_to_binary_vec(message);
    return encrypt_vector(&plaintext, &public_key);
}*/

pub fn decrypt(ciphertext:&Vec<pair>, private_key:&private_key) -> String {
    let mut plaintext = decrypt_vector(ciphertext, &private_key);
    let mut plaintext = binary_vec_to_string(&plaintext);
    //let plaintext = String::from("");
    plaintext
}

fn encrypt_vector(message:&Vec<u64>, public_key:&public_key) -> Vec<pair> {
    let mut encrypted_bits = Vec::new();
    for bit in message {
        //let mut encrypted_bit = encrypt_bit(*bit, &public_key);
        let (pair, error_sum) = encrypt_bit(*bit, &public_key);
        //encrypted_bits.push(encrypted_bit);
        encrypted_bits.push(pair);
    }
    encrypted_bits
}

fn decrypt_vector(encrypted_bits:&Vec<pair>, private_key:&private_key) -> Vec<u64> {
    let mut decrypted_bits = Vec::new();
    for bit in encrypted_bits {
        let mut decrypted_bit = decrypt_bit(bit, &private_key);
        decrypted_bits.push(decrypted_bit);
    }
    decrypted_bits
}

fn decrypt_bit(encrypted_bit:&pair, private_key:&private_key) -> u64 {
    let s = private_key.s.clone();
    let p = private_key.p;
    let mut a = encrypted_bit.a.clone();
    let mut b = encrypted_bit.b.clone();
    let mut inner_product = inner_product(&a, &s, p);
    let mut point = (b as i64 - inner_product as i64).rem_euclid(p as i64) as u64;
    let mut q = (p as f64 / 2 as f64).floor() as u64;
    let mut circular_distance_to_0 = circular_distance(point as u64, 0, p);
    let mut circular_distance_to_q = circular_distance(point as u64, q, p);
    if circular_distance_to_0 < circular_distance_to_q {
        return 0;
    } else {
        return 1;
    }
}

fn circular_distance(point:u64, target:u64, modulus:u64) -> u64 {
    let diff = (point as i64 - target as i64).abs() as u64;
    return std::cmp::min(diff, modulus-diff);
}

fn encrypt_bit(bit:u64, public_key:&public_key) -> (pair, Vec<u64>) {
    let m = public_key.m;
    let n = public_key.n;
    let p = public_key.p;
    let list = public_key.list.clone();
    let mut subset = choose_subset(m); 
    let mut x = vec![0; n as usize];
    let mut y = 0 as u64;
    for index in subset.clone() {
        let mut a = list[index as usize].a.clone();
        let mut b = list[index as usize].b.clone();
        x = add_vector(&x, &a, p);
        y += b; 
        y = y.rem_euclid(p);
        
    }
    if bit == 0 {
        let pair = pair {
            a: x,
            b: y
        };
        return (pair, subset);
    } else {
        y += (p as f64 / 2 as f64).floor() as u64;
        y = y.rem_euclid(p);
        let pair = pair {
            a: x,
            b: y
        };
        return (pair, subset);
    }
}

fn add_vector(a:&Vec<u64>, b:&Vec<u64>, p:u64) -> Vec<u64> {
    let mut sum = Vec::new();
    for i in 0..a.len() {
        sum.push((a[i] + b[i]).rem_euclid(p));
    }
    sum
}

fn choose_subset(m:u64) -> Vec<u64> {
    let mut set = Vec::new();
    let mut num = Vec::new();
    while num.len() < m as usize { // create length m vector of 1s and 0s
        num.push(rand::thread_rng().gen_range(0..2) as u8); // 0 or 1
    }
    for index in 0..num.len() {
        if num[index]==1 {
            set.push((index) as u64);
        }
    }
    set
}

fn get_list(m:u64, n:u64, p:u64, B:f64, private_key:&Vec<u64>) -> (Vec<pair>, Vec<u64>) {
    let mut s = private_key.clone();
    let mut public_key = Vec::new();
    let mut e_vector = Vec::new();
    for i in 0..m {
        let mut a = draw_vector_from_Znp(n, p);
        let mut e = draw_element_from_X(p, B);
        e_vector.push(e);
        //let mut b = (inner_product(&a, &s, p) + draw_element_from_X(p, B)).rem_euclid(p); // reduce?
        let mut b = (inner_product(&a, &s, p) + e).rem_euclid(p); // reduce?
        let mut pair = pair {
            a: a, 
            b: b
        };
        public_key.push(pair);
    }
    (public_key, e_vector)
} 

fn inner_product(a:&Vec<u64>, b:&Vec<u64>, p:u64) -> u64 {
    let mut sum = 0;
    if a.len() != b.len() {
        println!("Error in inner product - lengths were different");
        process::exit(1);
    }
    for i in 0..a.len() {
        sum += a[i] * b[i];
    }
    return (sum).rem_euclid(p);
}

fn draw_vector_from_Znp(n:u64, p:u64) -> Vec<u64> {
    let mut vec = Vec::new();
    for i in 0..n {
        vec.push(rand::thread_rng().gen_range(0..p));
    }
    vec
}

fn draw_element_from_X(p:u64, B:f64) -> u64 {
    let mut mean = 0.0;
    let mut std_dev = B/(2 as f64*std::f64::consts::PI).sqrt(); // mean 0, standard deviation B/sqrt(2pi)
    let normal = Normal::new(mean, std_dev).unwrap(); // create normal distribution 
    let mut v = normal.sample(&mut rand::thread_rng()); // sample from normal distribution
    let mut w = (v).rem_euclid(1.0); // reduce result modulo 1
    let mut x = w * p as f64; // multiply by p
    let mut y = (x.round()).rem_euclid(p as f64); // round to nearest integer modulo p
    return y as u64;
}

// return B E R+
fn alpha(n:f64) -> f64 {
    return 1 as f64/(n.sqrt() * n.log2() * n.log2());
}


// returns vector of all possible binary vectors of width width
pub fn all_possible_messages(width:u64) -> Vec<Vec<u64>> {
    let mut vectors = Vec::new();
    for i in 0..2_u64.pow(width as u32) {
        let mut vec = vec(i, width);
        vectors.push(vec);
    }
    vectors
}
        
pub fn vec(u:u64, length:u64) -> Vec<u64> {
    let mut vec = vec![0 as u64; length as usize];
    let mut u = u.clone();
    for i in (0..vec.len()).rev() {
        if u%2 == 0 {
            vec[i] = 0 as u64;
        } else {
            vec[i] = 1 as u64;
        }
        u = u >> 1;// shift_right(&mut u);
    }
    vec
}

pub fn is_equal(a:&Vec<u64>, b:&Vec<u64>, width:u64) -> bool {
    for i in 0..width as usize {
        if a[i] != b[i] {
            return false;
        }
    }
    return true;
}
pub fn string_to_binary_vec(s:String) -> Vec<u64> {
    let mut binary_string = String::from("");
    for character in s.into_bytes() {
        if character == 32 {
            binary_string += &format!("00{:b}", character);
        } else {
            binary_string += &format!("0{:b}", character);
        }
    }
    let mut binary_vector = Vec::new();
    for ch in binary_string.chars() {
        if ch == '1' {
            binary_vector.push(1);
        } else {
            binary_vector.push(0);
        }
    }
    binary_vector
}

pub fn binary_vec_to_string(v:&Vec<u64>) -> String {
    let mut vec_of_chars = Vec::new();
    let mut binary_string = String::from("");
    for i in 0..v.len() {
        if v[i] == 0 {
            binary_string += "0";
        } else {
            binary_string += "1";
        }
        if binary_string.len() == 8 {
            vec_of_chars.push(binary_string.clone());
            binary_string = String::from("");
        }
    }
    let mut w = Vec::new();
    for i in 0..vec_of_chars.len() {
        let bin = vec_of_chars[i].clone();
        w.push(isize::from_str_radix(&vec_of_chars[i].clone(), 2).unwrap() as u8);
    }
    while w[w.len()-1] == 0 {
        w.pop();
    }
    let mut w = &w[..];
    //let st = std::str::from_utf8(&w).unwrap();
    //st.to_string()
    let st = String::from("");
    st
}

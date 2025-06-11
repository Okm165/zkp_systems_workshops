pub mod tests;

// Finite Field Arithmetic over GF(p)
#[derive(Debug, Clone)]
pub struct PrimeField {
    pub p: u64,
}

impl PrimeField {
    pub fn new(p: u64) -> Self {
        assert!(Self::is_prime(p), "p must be prime");
        Self { p }
    }

    pub fn add(&self, a: u64, b: u64) -> u64 {
        (a + b) % self.p
    }

    pub fn sub(&self, a: u64, b: u64) -> u64 {
        (a + self.p - b) % self.p
    }

    pub fn mul(&self, a: u64, b: u64) -> u64 {
        (a * b) % self.p
    }

    pub fn neg(&self, a: u64) -> u64 {
        (self.p - a % self.p) % self.p
    }

    pub fn inv(&self, a: u64) -> u64 {
        assert!(a != 0, "No inverse for 0");
        Self::mod_pow(a, self.p - 2, self.p)
    }

    pub fn div(&self, a: u64, b: u64) -> u64 {
        self.mul(a, self.inv(b))
    }

    pub fn is_prime(n: u64) -> bool {
        if n <= 1 {
            return false;
        }
        if n <= 3 {
            return true;
        }
        if n % 2 == 0 || n % 3 == 0 {
            return false;
        }
        let mut i = 5;
        while i * i <= n {
            if n % i == 0 || n % (i + 2) == 0 {
                return false;
            }
            i += 6;
        }
        true
    }

    pub fn mod_pow(mut base: u64, mut exp: u64, modulus: u64) -> u64 {
        let mut result = 1;
        base %= modulus;
        while exp > 0 {
            if exp % 2 == 1 {
                result = (result * base) % modulus;
            }
            base = (base * base) % modulus;
            exp /= 2;
        }
        result
    }
}

// Polynomial over GF(p)
#[derive(Debug, Clone)]
pub struct Polynomial {
    coeffs: Vec<u64>,
    field: PrimeField,
}

impl Polynomial {
    pub fn new(coeffs: Vec<u64>, field: PrimeField) -> Self {
        let coeffs = coeffs.into_iter().map(|c| c % field.p).collect();
        Self { coeffs, field }
    }

    pub fn evaluate(&self, x: u64) -> u64 {
        let mut result = 0;
        for (power, &coeff) in self.coeffs.iter().enumerate() {
            let x_pow = PrimeField::mod_pow(x, power as u64, self.field.p);
            let term = self.field.mul(coeff, x_pow);
            result = self.field.add(result, term);
        }
        result
    }
}

fn main() {
    let f = PrimeField::new(7);
    println!("5 + 3 mod 7 = {}", f.add(5, 3)); // Output: 1
    println!("3 / 4 mod 7 = {}", f.div(3, 4)); // Output: 6

    let poly = Polynomial::new(vec![1, 0, 1], f.clone()); // x^2 + 1
    println!("P(0) = {}", poly.evaluate(0)); // 1
    println!("P(2) = {}", poly.evaluate(2)); // 2^2 + 1 = 5
}

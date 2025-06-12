use proptest::prelude::*;

use crate::PrimeField;

// Define the field GF(7)
pub const PRIME: u64 = (1 << 31) - 1;
pub const FIELD: PrimeField = PrimeField { p: PRIME };

// Closure property for addition and multiplication:
// The result of a + b and a * b must remain in the field (i.e., less than p)
proptest! {
    #[test]
    fn closure_add_mul(a in 0..FIELD.p, b in 0..FIELD.p) {
        let sum = FIELD.add(a, b);
        let product = FIELD.mul(a, b);
        prop_assert!(sum < FIELD.p);
        prop_assert!(product < FIELD.p);
    }
}

// Associativity property:
// (a + b) + c == a + (b + c)
// (a * b) * c == a * (b * c)
proptest! {
    #[test]
    fn associativity_add_mul(a in 0..FIELD.p, b in 0..FIELD.p, c in 0..FIELD.p) {
        prop_assert_eq!(FIELD.add(FIELD.add(a, b), c), FIELD.add(a, FIELD.add(b, c)));
        prop_assert_eq!(FIELD.mul(FIELD.mul(a, b), c), FIELD.mul(a, FIELD.mul(b, c)));
    }
}

// Commutativity property:
// a + b == b + a
// a * b == b * a
proptest! {
    #[test]
    fn commutativity_add_mul(a in 0..FIELD.p, b in 0..FIELD.p) {
        prop_assert_eq!(FIELD.add(a, b), FIELD.add(b, a));
        prop_assert_eq!(FIELD.mul(a, b), FIELD.mul(b, a));
    }
}

// Identity elements:
// a + 0 == a (additive identity)
// a * 1 == a (multiplicative identity)
proptest! {
    #[test]
    fn identity_elements(a in 0..FIELD.p) {
        prop_assert_eq!(FIELD.add(a, 0), a);
        prop_assert_eq!(FIELD.mul(a, 1), a);
    }
}

// Additive inverse:
// a + (-a) == 0
proptest! {
    #[test]
    fn additive_inverse(a in 0..FIELD.p) {
        prop_assert_eq!(FIELD.add(a, FIELD.neg(a)), 0);
    }
}

// Multiplicative inverse:
// a * a^-1 == 1, for all a != 0
proptest! {
    #[test]
    fn multiplicative_inverse(a in 1..FIELD.p) {
        let inv = FIELD.inv(a);
        prop_assert_eq!(FIELD.mul(a, inv), 1);
    }
}

// Distributivity:
// a * (b + c) == a * b + a * c
proptest! {
    #[test]
    fn distributivity(a in 0..FIELD.p, b in 0..FIELD.p, c in 0..FIELD.p) {
        let left = FIELD.mul(a, FIELD.add(b, c));
        let right = FIELD.add(FIELD.mul(a, b), FIELD.mul(a, c));
        prop_assert_eq!(left, right);
    }
}

proptest! {
    #[test]
    #[should_panic]
    fn associativity_add_float(a in 0_f64..10_f64, b in 0_f64..10_f64, c in 0_f64..10_f64) {
        prop_assert_eq!((a + b) + c, a + (b + c));
    }
}

proptest! {
    #[test]
    #[should_panic]
    fn associativity_mul_float(a in 0_f64..10_f64, b in 0_f64..10_f64, c in 0_f64..10_f64) {
        prop_assert_eq!((a * b) * c, a * (b * c));
    }
}

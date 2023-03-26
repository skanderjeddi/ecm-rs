use std::{
    env::args,
    fmt::{Display, Formatter},
    time::Instant,
};

use log::{debug, info};

use num_bigint::{BigInt, RandBigInt};
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use primes::{PrimeSet, Sieve};
use rayon::prelude::{IntoParallelIterator, ParallelIterator};

#[derive(Debug, Clone, Eq, Hash, PartialEq)]
struct Point {
    x: BigInt,
    y: BigInt,
}

impl Display for Point {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let (x, y) = (self.x.clone(), self.y.clone());
        write!(f, "({x}, {y})")
    }
}

#[derive(Debug, Clone, Eq, Hash, PartialEq)]
struct EllipticCurveFactorizationMethod {
    n: BigInt,
    a: BigInt,
    b: BigInt,
}

impl Display for EllipticCurveFactorizationMethod {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let (n, a, b) = (self.n.clone(), self.a.clone(), self.b.clone());
        write!(f, "({a}, {b}) [{n}])")
    }
}

impl EllipticCurveFactorizationMethod {
    fn new(n: BigInt) -> Self {
        let mut rng = rand::thread_rng();
        let x0 = rng.gen_bigint_range(&BigInt::from(2u8), &n);
        let y0 = rng.gen_bigint_range(&BigInt::from(2u8), &n);
        let a = rng.gen_bigint_range(&BigInt::from(2u8), &n);
        let y0_squared = y0.pow(2) % n.clone();
        let x0_cubed = x0.pow(3) % n.clone();
        let a_x0 = (a.clone() * x0) % n.clone();
        let mut first_term: BigInt = y0_squared - x0_cubed;
        while first_term <= Zero::zero() {
            first_term += n.clone();
        }
        let mut b: BigInt = first_term.clone() - a_x0;
        while b <= Zero::zero() {
            b += n.clone();
        }
        Self { n, a, b }
    }

    fn get_point_on_curve(&self) -> Point {
        let mut rng = rand::thread_rng();
        let x = rng.gen_bigint_range(&BigInt::from(2u32), &self.n);
        let y_squared = (x.pow(3) + self.a.clone() * x.clone() + self.b.clone()) % self.n.clone();
        let y = y_squared.sqrt();
        Point { x, y }
    }

    #[inline(always)]
    fn modular_inverse(&self, p: BigInt) -> BigInt {
        let m0 = self.n.clone();
        if m0 == One::one() {
            return One::one();
        }
        let (mut a, mut m, mut x0, mut inv) = (p, m0.clone(), Zero::zero(), One::one());
        while a > One::one() {
            let (div, rem) = a.div_rem(&m);
            inv -= div * &x0;
            a = rem;
            std::mem::swap(&mut a, &mut m);
            std::mem::swap(&mut x0, &mut inv)
        }
        if inv < Zero::zero() {
            inv += m0
        }
        inv
    }

    #[inline(always)]
    fn add_points(&self, p: &Point, q: &Point) -> (bool, Point) {
        let (x1, y1) = (p.x.clone(), p.y.clone());
        let (x2, y2) = (q.x.clone(), q.y.clone());
        let x3: BigInt;
        let y3: BigInt;
        let mut numerator: BigInt;
        let mut denominator: BigInt;
        if x1 == x2 {
            numerator = BigInt::from(3u8) * x1.pow(2) + self.a.clone();
            denominator = BigInt::from(2u8) * y1.clone();
        } else {
            numerator = y2 - y1.clone();
            while numerator < Zero::zero() {
                numerator += self.n.clone();
            }
            denominator = x2.clone() - x1.clone();
            while denominator < Zero::zero() {
                denominator += self.n.clone();
            }
        }
        let x = numerator.gcd(&self.n);
        if x != One::one() {
            x3 = x.clone();
            y3 = self.n.clone() / x;
            return (true, Point { x: x3, y: y3 });
        }
        let x = denominator.gcd(&self.n);
        if x != One::one() {
            x3 = x.clone();
            y3 = self.n.clone() / x;
            return (true, Point { x: x3, y: y3 });
        }
        let slope = (numerator.clone() * self.modular_inverse(denominator.clone())).abs();
        let slope = slope % self.n.clone();
        let mut rx = slope.pow(2) - x1.clone() - x2;
        while rx < Zero::zero() {
            rx += self.n.clone();
        }
        rx %= self.n.clone();
        let mut second_factor = x1 - rx.clone();
        while second_factor < Zero::zero() {
            second_factor += self.n.clone();
        }
        let mut ry = slope * (second_factor) - y1;
        while ry < Zero::zero() {
            ry += self.n.clone();
        }
        ry %= self.n.clone();
        (false, Point { x: rx, y: ry })
    }
}

fn setup() {
    if std::env::var_os("RUST_LOG").is_none() {
        std::env::set_var("RUST_LOG", "debug");
    }
    pretty_env_logger::init();
    color_backtrace::install();
}

fn main() {
    setup();

    let n = args()
        .nth(1)
        .expect("Usage: ecm-rs <number>")
        .parse::<BigInt>()
        .expect("Failed to parse number");

    info!(
        "Factorizing N = {} (using {} threads)...",
        n,
        rayon::current_num_threads()
    );

    let friability_factor = 128;

    debug!(
        "Finding primes up to friability factor ({})...",
        friability_factor
    );
    let mut sieve = Sieve::new();
    let primes = sieve
        .iter()
        .take_while(|p| p <= &friability_factor)
        .collect::<Vec<_>>();
    debug!("Found primes: {:?}.", primes);
    debug!(
        "Finding exponent of primes for friability factor ({})...",
        friability_factor
    );
    let mut primes_exponents = Vec::with_capacity(primes.len());
    for prime in &primes {
        let mut e = 1;
        while prime.pow(e) <= friability_factor {
            e += 1;
        }
        primes_exponents.push(e - 1);
    }
    debug!("Found exponents: {:?}.", primes_exponents);

    let primes_and_exponents = primes
        .iter()
        .zip(primes_exponents.iter())
        .collect::<Vec<_>>();

    info!("Starting factorization...");

    let total_time = Instant::now();
    let number_of_curves = 128;
    loop {
        let factor = (0..number_of_curves)
            .into_par_iter()
            .map(|_| {
                
                EllipticCurveFactorizationMethod::new(n.clone())
            })
            .into_par_iter()
            .map(|ecm| {
                let p = ecm.get_point_on_curve();
                let p = p;
                let mut q = p.clone();
                for (&prime, &e) in primes_and_exponents.iter() {
                    for _ in 0..prime.pow(e - 1) {
                        let (found, r) = ecm.add_points(&p, &q);
                        if found {
                            return r.x;
                        }
                        q = r.clone();
                    }
                }
                Zero::zero()
            })
            .find_any(|x| x.clone() != Zero::zero());

        if let Some(factor) = factor {
            let time = total_time.elapsed();
            info!("Found factor {factor} in {time:?}",);
            break;
        }
    }
}

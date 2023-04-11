use std::{
    env::args,
    fmt::{Display, Formatter},
    time::Instant,
};

use log::{debug, info};

use num_bigint::{BigInt, RandBigInt};
use num_integer::Integer;
use num_traits::{One, Zero};
use primes::{PrimeSet, Sieve};
use rayon::{
    current_thread_index,
    prelude::{IntoParallelIterator, ParallelIterator},
};

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
        let b: BigInt = (y0.pow(2) - x0.pow(3) - a.clone() * x0) % n.clone();
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
        inv
    }

    #[inline(always)]
    fn add_points(&self, p: &Point, q: &Point) -> (bool, Point) {
        let (x1, y1) = (p.x.clone(), p.y.clone());
        let (x2, y2) = (q.x.clone(), q.y.clone());
        let x3: BigInt;
        let y3: BigInt;
        let numerator: BigInt;
        let denominator: BigInt;
        if x1 == x2 {
            numerator = BigInt::from(3u8) * x1.pow(2) + self.a.clone();
            denominator = BigInt::from(2u8) * y1.clone();
        } else {
            numerator = (y2 - y1.clone()) % self.n.clone();
            denominator = (x2.clone() - x1.clone()) % self.n.clone();
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
        let slope =
            (numerator.clone() * self.modular_inverse(denominator.clone())) % self.n.clone();
        let rx = (slope.pow(2) - x1.clone() - x2) % self.n.clone();
        let second_factor = (x1 - rx.clone()) % self.n.clone();
        let ry = slope * (second_factor) - y1;
        (false, Point { x: rx, y: ry })
    }

    fn point_multiplication(&self, p: &Point, k: &BigInt) -> (bool, Point) {
        let q = p.clone();
        let mut r = Point {
            x: Zero::zero(),
            y: Zero::zero(),
        };
        let mut found_factor = false;
        let mut i = k.bits() - 1;
        while i > 0 {
            r = self.add_points(&r, &r).1;
            if k.bit(i) {
                let (found, result) = self.add_points(&r, &q);
                if found {
                    found_factor = true;
                    return (found_factor, result);
                }
                r = result;
            }
            i -= 1;
        }
        (found_factor, r)
    }
}

fn setup() {
    if std::env::var_os("RUST_LOG").is_none() {
        std::env::set_var("RUST_LOG", "info");
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

    /* let num_digits = n.to_str_radix(10).len() as f64;
    let b1 = num_digits.sqrt() * (n.bits() as f64);
    let b1 = E.powf((b1.ln() * (b1.ln()).ln()) as f32) as u64; */

    let b1 = 5 * 7 * 11 * 17 * 23;

    info!(
        "Factorizing N = {} (using {} threads, B1 = {})...",
        n,
        rayon::current_num_threads(),
        b1
    );

    let total_time = Instant::now();

    debug!("Finding primes up to friability factor ({})...", b1);

    let mut sieve = Sieve::new();
    let primes_iter = sieve.iter();

    /* let mut product = 1;
    while product < b1 {
        let prime = primes_iter.next().unwrap();
        product *= prime;
        primes.push(prime);
    } */

    let primes: Vec<u64> = primes_iter.take_while(|p| p <= &b1).collect();

    let mut primes_exponents = Vec::with_capacity(primes.len());
    for prime in &primes {
        let mut e = 1;
        while prime.pow(e) <= b1 {
            e += 1;
        }
        primes_exponents.push(e);
    }

    let primes_raised = primes
        .iter()
        .zip(primes_exponents.iter())
        .map(|(&p, &e)| BigInt::from(p.pow(e)))
        .collect::<Vec<_>>();

    debug!("{:?}", primes_raised);

    debug!("Starting factorization...");

    // let number_of_curves = b1 / ((b1 as f64).ln() as u64);

    let number_of_curves = 8;

    let mut iteration = 1;
    loop {
        info!(
            "Generating {} curves (iteration #{iteration})...",
            number_of_curves
        );
        let curves = (0..number_of_curves)
            .into_par_iter()
            .map(|_| EllipticCurveFactorizationMethod::new(n.clone()))
            .collect::<Vec<_>>();
        iteration += 1;

        let factor = curves
            .into_par_iter()
            .map(|ecm| {
                let now = Instant::now();
                let p = ecm.get_point_on_curve();
                let p = p;
                for k in primes_raised.iter() {
                    let (found, r) = ecm.point_multiplication(&p, k);
                    if found {
                        let elaped = Instant::now() - now;
                        info!(
                            "Thread #{} found factor {} in {:?}",
                            current_thread_index().unwrap(),
                            r.x,
                            elaped
                        );
                        return (r.x, elaped);
                    }
                }
                /*for k in primes_raised.iter() {:?
                    for _ in 0..prime.pow(e - 1) {
                        let (found, r) = ecm.add_points(&p, &q);
                        if found {
                            return r.x;
                        }
                        q = r.clone();
                    }
                    while k < prime.pow(e - 1) {
                        let (found, r) = ecm.add_points(&p, &q);
                        if found {
                            return r.x;
                        }
                        p = r.clone();
                        q = r.clone();
                        k *= 2;
                    }
                    debug!("Computing kP for k = {:?} and P = {:?} on curve {}", k, p, &ecm);
                    let (found, r) = ecm.point_multiplication(&p, k);
                    if found {
                        return r.x;
                    }
                    p = r.clone();
                    debug!("Done. (kP = {:?} on curve {})", p, &ecm);
                }*/
                let elapsed = Instant::now() - now;
                info!(
                    "Curve done in {:?} (thread #{}).",
                    elapsed,
                    current_thread_index().unwrap()
                );
                (Zero::zero(), elapsed)
            })
            .find_any(|x| x.0.clone() != Zero::zero());

        if let Some(factor) = factor {
            let total_time = Instant::now() - total_time;
            info!(
                "Found factor {}, total time elapsed was {:?}",
                factor.0, total_time
            );
            info!("N = {} * {}", factor.0.clone(), n / factor.0);
            break;
        } else {
            debug!("Exhausted all curves, generating new ones...");
        }
    }
}

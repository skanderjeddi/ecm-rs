use std::time::Instant;

use log::{debug, info};

use num_bigint::{BigInt, BigUint, RandBigInt, ToBigInt, ToBigUint};
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

type ECMPoint = (BigUint, BigUint);

struct ECM {
    n: BigUint,
    a: BigUint,
    b: BigUint,
}

impl ECM {
    fn new(n: BigUint) -> Self {
        info!("Constructing ECM with parameter n = {}", n);
        let mut rng = rand::thread_rng();
        let x0 = rng.gen_biguint_range(&BigUint::from(2u32), &n);
        let y0 = rng.gen_biguint_range(&BigUint::from(2u32), &n);
        let a = rng.gen_biguint_range(&BigUint::from(2u32), &n);
        let y0s = (y0.clone() * y0.clone()) % n.clone();
        let x0c = (x0.clone() * x0.clone() * x0.clone()) % n.clone();
        let ax0 = (a.clone() * x0.clone()) % n.clone();
        let mut y0smx0c: BigInt =
            y0s.clone().to_bigint().unwrap() - x0c.clone().to_bigint().unwrap();
        while y0smx0c <= 0.to_bigint().unwrap() {
            y0smx0c += n.clone().to_bigint().unwrap();
        }
        let y0smx0c = y0smx0c.to_biguint().unwrap();
        let mut b: BigInt = y0smx0c.clone().to_bigint().unwrap() - ax0.clone().to_bigint().unwrap();
        while b <= 0.to_bigint().unwrap() {
            b += n.clone().to_bigint().unwrap();
        }
        let b = b.to_biguint().unwrap();
        debug!("x0 = {} [n]", x0);
        debug!("y0 = {} [n]", y0);
        debug!("a = {} [n]", a);
        debug!("b = (y_0)^2 - (x_0)^3 - a * x_0 = {} [n]", b);
        info!("ECM has equation y^2 = x^3 + {}x + {} (mod {})", a, b, n);
        Self { n, a, b }
    }

    fn get_point_on_curve(&self) -> ECMPoint {
        let mut rng = rand::thread_rng();
        let x = rng.gen_biguint_range(&BigUint::from(2u32), &self.n);
        let y_squared =
            (x.clone() * x.clone() * x.clone() + self.a.clone() * x.clone() + self.b.clone())
                % self.n.clone();
        let y = y_squared.sqrt();
        (x, y)
    }

    fn modular_inverse(&self, p: BigInt) -> BigInt {
        let now = Instant::now();
        let p0 = p.clone();
        debug!("Finding modular inverse of {} [{}]", p, self.n);
        let m0 = self.n.clone().to_bigint().unwrap();
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
        debug!(
            "Modular inverse of {} [{}] is {} (took {}ns)",
            p0,
            self.n,
            inv,
            now.elapsed().as_nanos()
        );
        inv
    }

    fn add_points(&self, p: &ECMPoint, q: &ECMPoint) -> (bool, ECMPoint) {
        debug!(
            "Adding points P = ({}, {}) and Q = ({}, {})",
            p.0, p.1, q.0, q.1
        );
        let (x1, y1) = p;
        let (x2, y2) = q;
        let x3: BigUint;
        let y3: BigUint;
        let numerator: BigInt;
        let denominator: BigInt;
        if x1 == x2 {
            numerator = (3.to_biguint().unwrap() * x1.clone() * x1.clone() + self.a.clone())
                .to_bigint()
                .unwrap();
            denominator = (2.to_biguint().unwrap() * y1.clone()).to_bigint().unwrap();
        } else {
            let mut y2c = y2.clone().to_bigint().unwrap();
            let y1c = y1.clone().to_bigint().unwrap();
            while y2c < y1c {
                debug!("y2 < y1, adding n to y2");
                y2c += self.n.clone().to_bigint().unwrap();
            }
            numerator = (y2c - y1c).to_bigint().unwrap();
            let mut x2c = x2.clone().to_bigint().unwrap();
            let x1c = x1.clone().to_bigint().unwrap();
            while x2c < x1c {
                debug!("y2 < y1, adding n to y2");
                x2c += self.n.clone().to_bigint().unwrap();
            }
            denominator = (x2c - x1c).to_bigint().unwrap();
        }
        let x = numerator.gcd(&self.n.to_bigint().unwrap());
        if x != 1.to_bigint().unwrap() {
            x3 = x.to_biguint().unwrap();
            y3 = self.n.clone() / x.to_biguint().unwrap();
            return (true, (x3, y3));
        }
        let x = denominator.gcd(&self.n.to_bigint().unwrap());
        if x != 1.to_bigint().unwrap() {
            x3 = x.to_biguint().unwrap();
            y3 = self.n.clone() / x.to_biguint().unwrap();
            return (true, (x3, y3));
        }
        let slope = (numerator.clone()
            * self
                .modular_inverse(denominator.clone())
                .clone()
                .to_bigint()
                .unwrap())
        .abs();
        let slope = slope.to_biguint().unwrap() % self.n.clone();
        debug!("numerator = {} [n]", numerator);
        debug!("denominator = {} [n]", denominator);
        debug!("slope = {} [n]", slope);
        // Rx = (slope^2 - x1 - x2) mod n
        let mut slope_squared = slope.clone() * slope.clone();
        while slope_squared < x1.clone() + x2.clone() {
            slope_squared += self.n.clone();
            debug!("adding n to slope_squared");
        }
        x3 = (slope_squared - x1.clone() - x2.clone()) % self.n.clone();
        debug!("x3 = {} [n]", x3);
        // Ry = (slope * (x1 - Rx) - y1) mod n
        let mut sx1 = (slope.clone() * x1.clone()) % self.n.clone();
        let sx3 = (slope.clone() * x3.clone()) % self.n.clone();
        while sx1 < sx3.clone() + y1.clone() {
            sx1 += self.n.clone();
            debug!("adding n to sx1");
        }
        y3 = (sx1 - sx3.clone() - y1.clone()) % self.n.clone();
        debug!("y3 = {} [n]", y3);
        return (false, (x3, y3));
    }
}

fn main() {
    if std::env::var_os("RUST_LOG").is_none() {
        std::env::set_var("RUST_LOG", "info");
    }
    pretty_env_logger::init();
    color_backtrace::install();
    let prime1 = "7543049".parse::<BigUint>().unwrap();
    let prime2 = "546341".parse::<BigUint>().unwrap();
    let n = prime1 * prime2;
    let now = Instant::now();
    let ecm = ECM::new(n);
    let p = ecm.get_point_on_curve();
    info!("P(x, y) = {:?}", p);
    let p = p.clone();
    let mut q = p.clone();
    let b = 2 * 3 * 5 * 7 * 9 * 11 * 13 * 17 * 19 * 23;
    info!("Using B = {}", b);
    for k in 2..=b {
        let (found, r) = ecm.add_points(&p, &q);
        if found {
            info!("Found factor: {} (for {k}P)", r.0);
            break;
        }
        debug!("{k}P = {:?}", r);
        q = r.clone();
    }
    info!("Took {}ms", now.elapsed().as_millis());
}

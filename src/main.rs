use std::{
    fmt::{Display, Formatter},
    time::Instant,
};

use colored::Colorize;
use log::info;

use num_bigint::{BigInt, RandBigInt};
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

#[derive(Debug, Clone)]
struct ECMPoint {
    x: BigInt,
    y: BigInt,
}

impl Display for ECMPoint {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let (x, y) = (self.x.clone(), self.y.clone());
        write!(
            f,
            "({}, {})",
            x.to_string().bright_red(),
            y.to_string().bright_red()
        )
    }
}

struct ECM {
    n: BigInt,
    a: BigInt,
    b: BigInt,
}

impl ECM {
    fn new(n: BigInt) -> Self {
        let mut rng = rand::thread_rng();
        let x0 = rng.gen_bigint_range(&BigInt::from(2u8), &n);
        let y0 = rng.gen_bigint_range(&BigInt::from(2u8), &n);
        let a = rng.gen_bigint_range(&BigInt::from(2u8), &n);
        let y0s = y0.pow(2) % n.clone();
        let x0c = x0.pow(3) % n.clone();
        let ax0 = (a.clone() * x0.clone()) % n.clone();
        let mut y0smx0c: BigInt = y0s.clone() - x0c.clone();
        while y0smx0c <= Zero::zero() {
            y0smx0c += n.clone();
        }
        let mut b: BigInt = y0smx0c.clone() - ax0.clone();
        while b <= Zero::zero() {
            b += n.clone();
        }
        info!(
            "y^2 = x^3 + {}x + {} [{}]",
            a.to_string().bright_red(),
            b.to_string().bright_red(),
            n.to_string().green()
        );
        Self { n, a, b }
    }

    fn get_point_on_curve(&self) -> ECMPoint {
        let mut rng = rand::thread_rng();
        let x = rng.gen_bigint_range(&BigInt::from(2u32), &self.n);
        let y_squared = (x.pow(3) + self.a.clone() * x.clone() + self.b.clone()) % self.n.clone();
        let y = y_squared.sqrt();
        ECMPoint { x, y }
    }

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

    fn add_points(&self, p: &ECMPoint, q: &ECMPoint) -> (bool, ECMPoint) {
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
            numerator = y2.clone() - y1.clone();
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
            y3 = self.n.clone() / x.clone();
            return (true, ECMPoint { x: x3, y: y3 });
        }
        let x = denominator.gcd(&self.n);
        if x != One::one() {
            x3 = x.clone();
            y3 = self.n.clone() / x.clone();
            return (true, ECMPoint { x: x3, y: y3 });
        }
        let slope = (numerator.clone() * self.modular_inverse(denominator.clone()).clone()).abs();
        let slope = slope % self.n.clone();
        let mut rx = slope.pow(2) - x1.clone() - x2.clone();
        while rx < Zero::zero() {
            rx += self.n.clone();
        }
        rx = rx % self.n.clone();
        let mut x1mrx = x1.clone() - rx.clone();
        while x1mrx < Zero::zero() {
            x1mrx += self.n.clone();
        }
        let mut ry = slope.clone() * (x1mrx) - y1.clone();
        while ry < Zero::zero() {
            ry += self.n.clone();
        }
        ry = ry % self.n.clone();
        return (false, ECMPoint { x: rx, y: ry });
    }
}

fn main() {
    if std::env::var_os("RUST_LOG").is_none() {
        std::env::set_var("RUST_LOG", "info");
    }
    pretty_env_logger::init();
    color_backtrace::install();
    let prime1 = "4094231".parse::<BigInt>().unwrap();
    let prime2 = "97114109".parse::<BigInt>().unwrap();
    let n = prime1.clone() * prime2.clone();
    info!(
        "{} = {} * {} = {}",
        "N".green(),
        prime1.to_string().blue(),
        prime2.to_string().yellow(),
        n.to_string().green()
    );
    let mut i = 1;
    let total_time = Instant::now();
    loop {
        let mut factor_found = false;
        info!("Iteration {}", i.to_string().underline());
        let ecm = ECM::new(n.clone());
        let p = ecm.get_point_on_curve();
        let p = p.clone();
        let mut q = p.clone();
        let b = 2 * 3 * 5 * 7 * 11 * 13 * 17;
        info!(
            "({}, {}) = (({}, {}), {})",
            "P".cyan(),
            "B".magenta(),
            p.x.to_string().cyan(),
            p.y.to_string().cyan(),
            b.to_string().bright_magenta()
        );
        for _ in 2..=b {
            let (found, r) = ecm.add_points(&p, &q);
            if found {
                let factor = r.x;
                if factor == prime1 {
                    info!("Found factor: {}", factor.to_string().blue());
                    factor_found = true;
                    break;
                } else {
                    info!("Found factor: {}", factor.to_string().yellow());
                    factor_found = true;
                    break;
                }
            }
            q = r.clone();
        }
        i += 1;
        if !factor_found {
            info!("No factors found with this pair (P, B)");
            continue;
        } else {
            break;
        }
    }
    info!(
        "Total time: {}ms",
        total_time.elapsed().as_millis().to_string().bright_green()
    );
}

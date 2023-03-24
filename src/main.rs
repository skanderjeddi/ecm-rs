use std::time::Instant;

use num_bigint::{BigInt, BigUint, RandBigInt, ToBigInt};
use log::{info, debug};

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
        Self { n, a, b }
    }

    fn get_point_on_curve(self) -> ECMPoint {
        let mut rng = rand::thread_rng();
        let x = rng.gen_biguint_range(&BigUint::from(2u32), &self.n);
        let y_squared =
            (x.clone() * x.clone() * x.clone() + self.a.clone() * x.clone() + self.b.clone())
                % self.n.clone();
        let y = y_squared.sqrt();
        (x, y)
    }
}

fn main() {
    if std::env::var_os("RUST_LOG").is_none() {
        std::env::set_var("RUST_LOG", "debug");
    }
    pretty_env_logger::init();
    let now = Instant::now();
    let ecm = ECM::new(
        "9470658081841584084540266594803402758422182767579119501969461691039581334408451838581104318520589571"
            .parse::<BigUint>()
            .unwrap(),
    );
    let p = ecm.get_point_on_curve();
    info!("P(x, y) = {:?}", p);
    info!("Time elapsed: {}ns", now.elapsed().as_nanos());
}

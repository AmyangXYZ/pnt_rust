use crate::gnss;
use ndarray::{Array1, Array2, Axis};

pub struct Satellite {
    pub id: u8,
    pub name: String,
    pub states: Vec<gnss::State>,
}

impl Satellite {
    pub fn new(id: u8, name: String) -> Self {
        Self {
            id,
            name,
            states: vec![],
        }
    }

    pub fn propagate(
        &mut self,
        start: std::time::SystemTime,
        duration: std::time::Duration,
        step: std::time::Duration,
        ephemeris_data: &[gnss::NavRecord],
    ) -> usize {
        // let end = start + duration;
        // let mut current = start;

        let gps_times: Array1<f64> = Array1::from_iter(
            (0..((duration.as_secs_f64() / step.as_secs_f64()) as usize))
                .map(|i| gnss::calculate_gps_time(start + step * i as u32) as f64),
        );

        // Convert ephemeris_data to Array2<f64>
        let ephem_array: Array2<f64> =
            Array2::from_shape_fn((16, ephemeris_data.len()), |(i, j)| match i {
                0 => ephemeris_data[j].sqrt_a,
                1 => ephemeris_data[j].eccentricity,
                2 => ephemeris_data[j].i0,
                3 => ephemeris_data[j].omega0,
                4 => ephemeris_data[j].omega,
                5 => ephemeris_data[j].m0,
                6 => ephemeris_data[j].toe,
                7 => ephemeris_data[j].delta_n,
                8 => ephemeris_data[j].omega_dot,
                9 => ephemeris_data[j].idot,
                10 => ephemeris_data[j].cus,
                11 => ephemeris_data[j].cuc,
                12 => ephemeris_data[j].crs,
                13 => ephemeris_data[j].crc,
                14 => ephemeris_data[j].cis,
                15 => ephemeris_data[j].cic,
                _ => unreachable!(),
            });

        let time_diff =
            (&ephem_array.row(6).to_owned().insert_axis(Axis(1)) - &gps_times).mapv(f64::abs);
        let closest_time_idx = time_diff
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .map(|(index, _)| index)
            .unwrap();

        let ephem = ephem_array.select(Axis(1), &[closest_time_idx]);
        println!("{:?}", ephem);
        // let a = ephem.row(0).mapv(|x| x.powi(2));
        // let e = ephem.row(1);
        // let i0 = ephem.row(2);
        // let omega0 = ephem.row(3);
        // let omega = ephem.row(4);
        // let m0 = ephem.row(5);
        // let toe = ephem.row(6);
        // let delta_n = ephem.row(7);
        // let omega_dot = ephem.row(8);
        // let idot = ephem.row(9);
        // let cus = ephem.row(10);
        // let cuc = ephem.row(11);
        // let crs = ephem.row(12);
        // let crc = ephem.row(13);
        // let cis = ephem.row(14);
        // let cic = ephem.row(15);

        // let tk = &gps_times - &toe;
        // let half_week = 302400.0;
        // let tk = (&tk + half_week) % (2.0 * half_week) - half_week;

        // let mu_earth = 3.986005e14; // Earth's gravitational constant
        // let n0 = (mu_earth / a.mapv(|x| x.powf(3.0))).mapv(f64::sqrt);
        // let n = &n0 + &delta_n;
        // let m = &m0 + &n * &tk;
        // let e_solved = solve_kepler_robust(&m, &e.to_owned());

        // let nu = (&e_solved.mapv(|x| 1.0 - x.powi(2)).mapv(f64::sqrt) * e_solved.mapv(f64::sin))
        //     .zip_map(&(e_solved.mapv(f64::cos) - &e), |y, x| y.atan2(x));
        // let phi = &nu + &omega;
        // let r = &a * (1.0 - &e * e.mapv(f64::cos));

        // let delta_u = &cus * (2.0 * &phi).mapv(f64::sin) + &cuc * (2.0 * &phi).mapv(f64::cos);
        // let delta_r = &crs * (2.0 * &phi).mapv(f64::sin) + &crc * (2.0 * &phi).mapv(f64::cos);
        // let delta_i = &cis * (2.0 * &phi).mapv(f64::sin) + &cic * (2.0 * &phi).mapv(f64::cos);

        // let u = &phi + &delta_u;
        // let r = &r + &delta_r;
        // let i = &i0 + &delta_i + &idot * &tk;

        // let x = &r * u.mapv(f64::cos);
        // let y = &r * u.mapv(f64::sin);

        // let omega = &omega0 + (&omega_dot - gnss::OMEGA_EARTH) * &tk - gnss::OMEGA_EARTH * &toe;

        // let x = &x * omega.mapv(f64::cos) - &y * &i.mapv(f64::cos) * omega.mapv(f64::sin);
        // let y = &x * omega.mapv(f64::sin) + &y * &i.mapv(f64::cos) * omega.mapv(f64::cos);
        // let z = &y * i.mapv(f64::sin);

        // for idx in 0..gps_times.len() {
        //     let state = gnss::State {
        //         time: vec![gps_times[idx]],
        //         position: vec![gnss::ECEF::new(x[idx], y[idx], z[idx])],
        //     };
        //     self.states.push(state);
        // }

        self.states.len()
    }

    // Helper function to solve Kepler's equation
    pub fn solve_kepler_robust(m: &Array1<f64>, e: &Array1<f64>) -> Array1<f64> {
        let max_iter = 30;
        let tolerance = 1e-8;

        let mut e_array = m.clone();
        for _ in 0..max_iter {
            let e_next = m + e * e_array.mapv(f64::sin);
            if (&e_next - &e_array).mapv(f64::abs).sum() < tolerance {
                return e_next;
            }
            e_array = e_next;
        }
        e_array
    }
}

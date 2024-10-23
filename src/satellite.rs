use crate::gnss;
use chrono::{DateTime, Utc};
use ndarray::{Array1, Array2, ArrayView1};

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
        start: DateTime<Utc>,
        duration: std::time::Duration,
        step: std::time::Duration,
        ephemeris_data: &[gnss::NavRecord],
    ) -> usize {
        let gps_times: Array1<f64> = Array1::from_iter(
            (0..((duration.as_millis() / step.as_millis()) as usize))
                .map(|i| gnss::calculate_gps_time((start + step * i as u32).into()) / 1000.0),
        );

        let closest_indices = gps_times.mapv(|time| {
            ephemeris_data
                .iter()
                .enumerate()
                .min_by_key(|(_, record)| {
                    let time_diff = (record.gps_millis - time).abs();
                    (time_diff * 1000.0) as u64
                })
                .map(|(index, _)| index)
                .unwrap_or(0)
        });

        let ephem = Array2::from_shape_fn((16, gps_times.len()), |(param, time_idx)| {
            let nav_record = &ephemeris_data[closest_indices[time_idx]];
            match param {
                0 => nav_record.sqrt_a,
                1 => nav_record.eccentricity,
                2 => nav_record.i0,
                3 => nav_record.omega0,
                4 => nav_record.omega,
                5 => nav_record.m0,
                6 => nav_record.toe,
                7 => nav_record.delta_n,
                8 => nav_record.omega_dot,
                9 => nav_record.idot,
                10 => nav_record.cus,
                11 => nav_record.cuc,
                12 => nav_record.crs,
                13 => nav_record.crc,
                14 => nav_record.cis,
                15 => nav_record.cic,
                _ => unreachable!(),
            }
        });

        let a = ephem.row(0).mapv(|x| x.powi(2));
        let e = ephem.row(1);
        let i0 = ephem.row(2);
        let omega0 = ephem.row(3);
        let omega = ephem.row(4);
        let m0 = ephem.row(5);
        let toe = ephem.row(6);
        let delta_n = ephem.row(7);
        let omega_dot = ephem.row(8);
        let idot = ephem.row(9);
        let cus = ephem.row(10);
        let cuc = ephem.row(11);
        let crs = ephem.row(12);
        let crc = ephem.row(13);
        let cis = ephem.row(14);
        let cic = ephem.row(15);
        let tk = &gps_times - &toe;
        let half_week = 302400.0;
        let tk = (&tk + half_week) % (2.0 * half_week) - half_week;
        let n0 = a.mapv(|a_val| (gnss::MU_EARTH / a_val.powi(3)).sqrt());
        let n = &n0 + &delta_n;
        let m = &m0 + &n * &tk;
        let e_array = Self::solve_kepler_robust(&m.view(), &e);

        let sin_e = e_array.mapv(f64::sin);
        let cos_e = e_array.mapv(f64::cos);
        let sqrt_1_minus_e2 = (1.0 - &e * &e).mapv(f64::sqrt);
        let nu = (&sqrt_1_minus_e2 * &sin_e)
            .iter()
            .zip((&cos_e - &e).iter())
            .map(|(&y, &x)| y.atan2(x))
            .collect::<Array1<f64>>();
        let phi = &nu + &omega;

        // Radius and argument of latitude correction
        let r = &a * (1.0 - &e * &cos_e);

        let phi_2 = &phi * 2.0;
        let sin_2phi = phi_2.mapv(f64::sin);
        let cos_2phi = phi_2.mapv(f64::cos);

        let delta_u = &cus * &sin_2phi + &cuc * &cos_2phi;
        let delta_r = &crs * &sin_2phi + &crc * &cos_2phi;
        let delta_i = &cis * &sin_2phi + &cic * &cos_2phi;

        // Corrected radius and argument of latitude
        let u = &phi + &delta_u;
        let r = &r + &delta_r;
        let i = &i0 + &delta_i + &idot * &tk;

        // Position in orbital plane
        let cos_u = u.mapv(f64::cos);
        let sin_u = u.mapv(f64::sin);
        let x = &r * &cos_u;
        let y = &r * &sin_u;

        // Earth-rotation correction
        let omega = &omega0 + (&omega_dot - gnss::OMEGA_E_DOT) * &tk - gnss::OMEGA_E_DOT * &toe;
        let cos_omega = omega.mapv(f64::cos);
        let sin_omega = omega.mapv(f64::sin);
        let cos_i = i.mapv(f64::cos);
        let sin_i = i.mapv(f64::sin);

        let x_ecef = &x * &cos_omega - &y * &cos_i * &sin_omega;
        let y_ecef = &x * &sin_omega + &y * &cos_i * &cos_omega;
        let z_ecef = &y * &sin_i;

        // Store states
        self.states.clear();
        for idx in 0..gps_times.len() {
            let state = gnss::State {
                time: vec![gps_times[idx]],
                position: vec![gnss::ECEF::new(x_ecef[idx], y_ecef[idx], z_ecef[idx])],
            };
            self.states.push(state);
        }
        println!("{:?}", self.states[0].position[0]);
        self.states.len()
    }

    fn solve_kepler_robust(m: &ArrayView1<f64>, e: &ArrayView1<f64>) -> Array1<f64> {
        let max_iter = 30;
        let tolerance = 1e-8;

        let mut e_array = m.to_owned();
        for _ in 0..max_iter {
            let e_next = m + e * &e_array.mapv(f64::sin);
            if (&e_next - &e_array).mapv(|x| x.abs()).sum() < tolerance {
                return e_next;
            }
            e_array = e_next;
        }
        e_array
    }
}

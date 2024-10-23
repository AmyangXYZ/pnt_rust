use chrono::{DateTime, TimeZone, Utc};
use std::fs::File;
use std::io::{BufRead, BufReader};

pub const OMEGA_E_DOT: f64 = 7.2921151467e-5; // WGS-84 earth rotation rate, rad/s
pub const MU_EARTH: f64 = 398600.5e9; // Earth's gravitational constant
pub const C_LIGHT: f64 = 299792458.0; // Speed of light, m/s

#[derive(Debug, PartialEq, Default, Clone, Copy)]
pub struct ECEF {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl ECEF {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
    pub fn to_lla(&self) -> LLA {
        LLA {
            latitude: 0.0,
            longitude: 0.0,
            altitude: 0.0,
        }
    }
}

pub struct LLA {
    pub latitude: f64,
    pub longitude: f64,
    pub altitude: f64,
}

impl LLA {
    pub fn new(latitude: f64, longitude: f64, altitude: f64) -> Self {
        Self {
            latitude,
            longitude,
            altitude,
        }
    }
    pub fn to_ecef(&self) -> ECEF {
        ECEF {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }
}

pub struct State {
    pub time: Vec<f64>,
    pub position: Vec<ECEF>,
}

impl State {
    pub fn new() -> Self {
        Self {
            time: vec![0.0],
            position: vec![ECEF::new(0.0, 0.0, 0.0)],
        }
    }
}

/// Calculate GPS time: milliseconds since GPS epoch (Jan 6, 1980) plus leap seconds
pub fn calculate_gps_time(time: std::time::SystemTime) -> f64 {
    let utc_time: DateTime<Utc> = time.into();
    let gps_epoch: DateTime<Utc> = Utc.with_ymd_and_hms(1980, 1, 6, 0, 0, 0).unwrap();
    let leap_seconds = 18.0; // As of 2024
    ((utc_time - gps_epoch).num_microseconds().unwrap() as f64 / 1e6 + leap_seconds) * 1000.0
}

#[derive(Debug, PartialEq, Default, Clone, Copy)]
pub struct NavRecord {
    pub sat_id: u8,
    pub epoch: (i32, i32, i32, i32, i32, i32),
    pub gps_millis: f64,
    pub sv_clock_bias: f64,
    pub sv_clock_drift: f64,
    pub sv_clock_drift_rate: f64,
    pub iode: f64,
    pub crs: f64,
    pub delta_n: f64,
    pub m0: f64,
    pub cuc: f64,
    pub eccentricity: f64,
    pub cus: f64,
    pub sqrt_a: f64,
    pub toe: f64,
    pub cic: f64,
    pub omega0: f64,
    pub cis: f64,
    pub i0: f64,
    pub crc: f64,
    pub omega: f64,
    pub omega_dot: f64,
    pub idot: f64,
    pub codes_on_l2_channel: f64,
    pub gps_week: f64,
    pub l2_p_data_flag: f64,
    pub sv_accuracy: f64,
    pub sv_health: f64,
    pub tgd: f64,
    pub iodc: f64,
    pub transmission_time: f64,
    pub fit_interval: f64,
}

pub struct RinexNav {
    pub records: Vec<NavRecord>,
}

impl RinexNav {
    pub fn from_file(filename: &str) -> Self {
        let mut records = Vec::new();
        let file = File::open(filename).expect("Failed to open file");
        let reader = BufReader::new(file);
        let mut lines = reader.lines();

        // Skip header
        while let Some(Ok(line)) = lines.next() {
            if line.contains("END OF HEADER") {
                break;
            }
        }

        // Parse records
        while let Some(Ok(line)) = lines.next() {
            if line.len() < 79 {
                continue;
            }

            let sat_id = line[1..3].trim().parse().unwrap_or(0);
            let epoch = Self::parse_epoch(&line[3..23]);
            let gps_millis = Self::epoch_to_gps_millis(&epoch);
            let sv_clock_bias = Self::parse_float(&line[23..42]);
            let sv_clock_drift = Self::parse_float(&line[42..61]);
            let sv_clock_drift_rate = Self::parse_float(&line[61..80]);

            let mut record = NavRecord {
                sat_id,
                epoch,
                gps_millis,
                sv_clock_bias,
                sv_clock_drift,
                sv_clock_drift_rate,
                ..Default::default()
            };

            // Parse additional lines
            let mut line_count = 0;
            for _ in 0..7 {
                if let Some(Ok(data_line)) = lines.next() {
                    Self::parse_data_line(&mut record, &data_line, line_count);
                    line_count += 1;
                }
            }
            records.push(record);
        }
        Self { records }
    }

    fn parse_epoch(s: &str) -> (i32, i32, i32, i32, i32, i32) {
        let parts: Vec<&str> = s.split_whitespace().collect();
        (
            parts[0].parse().unwrap_or(0),
            parts[1].parse().unwrap_or(0),
            parts[2].parse().unwrap_or(0),
            parts[3].parse().unwrap_or(0),
            parts[4].parse().unwrap_or(0),
            parts[5].parse().unwrap_or(0),
        )
    }

    fn epoch_to_gps_millis(epoch: &(i32, i32, i32, i32, i32, i32)) -> f64 {
        let utc_time = Utc
            .with_ymd_and_hms(
                epoch.0,
                epoch.1 as u32,
                epoch.2 as u32,
                epoch.3 as u32,
                epoch.4 as u32,
                epoch.5 as u32,
            )
            .unwrap();
        calculate_gps_time(utc_time.into())
    }

    fn parse_float(s: &str) -> f64 {
        s.trim().replace('D', "E").parse().unwrap_or(0.0)
    }

    fn parse_data_line(record: &mut NavRecord, line: &str, line_number: usize) {
        let values: Vec<f64> = line[4..]
            .chars()
            .collect::<Vec<char>>()
            .chunks(19)
            .take(4)
            .map(|chunk| {
                chunk
                    .iter()
                    .collect::<String>()
                    .replace('D', "E")
                    .trim()
                    .parse()
                    .unwrap_or(0.0)
            })
            .collect();

        match line_number {
            0 => {
                record.iode = values.get(0).copied().unwrap_or(0.0);
                record.crs = values.get(1).copied().unwrap_or(0.0);
                record.delta_n = values.get(2).copied().unwrap_or(0.0);
                record.m0 = values.get(3).copied().unwrap_or(0.0);
            }
            1 => {
                record.cuc = values.get(0).copied().unwrap_or(0.0);
                record.eccentricity = values.get(1).copied().unwrap_or(0.0);
                record.cus = values.get(2).copied().unwrap_or(0.0);
                record.sqrt_a = values.get(3).copied().unwrap_or(0.0);
            }
            2 => {
                record.toe = values.get(0).copied().unwrap_or(0.0);
                record.cic = values.get(1).copied().unwrap_or(0.0);
                record.omega0 = values.get(2).copied().unwrap_or(0.0);
                record.cis = values.get(3).copied().unwrap_or(0.0);
            }
            3 => {
                record.i0 = values.get(0).copied().unwrap_or(0.0);
                record.crc = values.get(1).copied().unwrap_or(0.0);
                record.omega = values.get(2).copied().unwrap_or(0.0);
                record.omega_dot = values.get(3).copied().unwrap_or(0.0);
            }
            4 => {
                record.idot = values.get(0).copied().unwrap_or(0.0);
                record.codes_on_l2_channel = values.get(1).copied().unwrap_or(0.0);
                record.gps_week = values.get(2).copied().unwrap_or(0.0);
                record.l2_p_data_flag = values.get(3).copied().unwrap_or(0.0);
            }
            5 => {
                record.sv_accuracy = values.get(0).copied().unwrap_or(0.0);
                record.sv_health = values.get(1).copied().unwrap_or(0.0);
                record.tgd = values.get(2).copied().unwrap_or(0.0);
                record.iodc = values.get(3).copied().unwrap_or(0.0);
            }
            6 => {
                record.transmission_time = values.get(0).copied().unwrap_or(0.0);
                record.fit_interval = values.get(1).copied().unwrap_or(0.0);
            }
            _ => {}
        }
    }
}

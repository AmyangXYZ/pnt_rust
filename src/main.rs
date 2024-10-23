use chrono::{DateTime, Utc};
use pnt_rust::{gnss::RinexNav, satellite::Satellite};

fn main() {
    let sat_id: u8 = 17;
    let mut satellite = Satellite::new(sat_id, String::from("ISS"));
    let start = std::time::SystemTime::now();

    let duration = std::time::Duration::from_secs(1);
    let step = std::time::Duration::from_millis(1);

    let nav_data = RinexNav::from_file("constellation/GCGO00USA_R_20231630000_01D_GN.rnx");

    let ephemeris_data: Vec<_> = nav_data
        .records
        .clone()
        .into_iter()
        .filter(|record| record.sat_id == sat_id)
        .collect();

    println!("Total records: {}", nav_data.records.len());
    println!("Filtered records for {}: {}", sat_id, ephemeris_data.len());

    let begin_time = std::time::SystemTime::now();
    let n_states = satellite.propagate(start, duration, step, &ephemeris_data);
    let end_time = std::time::SystemTime::now();
    let execution_time = end_time.duration_since(begin_time).unwrap();

    let start_datetime: DateTime<Utc> = start.into();
    let end_datetime: DateTime<Utc> = (start + duration).into();

    println!(
        "Propagated {} states from {} to {} in {:.3} ms",
        n_states,
        start_datetime.format("%Y-%m-%d %H:%M:%S%.3f UTC"),
        end_datetime.format("%Y-%m-%d %H:%M:%S%.3f UTC"),
        execution_time.as_micros() as f64 / 1_000.0
    );
}

use std::io::{self, BufRead, BufReader, Write};
use std::collections::HashMap;
use std::sync::{Arc, Mutex};
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;

const CHUNK_SIZE: usize = 40_000; // Process 40K lines (10K reads) at a time

fn main() -> io::Result<()> {
    let stdin = io::stdin();
    let mut reader = BufReader::new(stdin);

    // Check if input is gzipped by reading magic bytes
    let mut magic = [0u8; 2];
    {
        let buf = reader.fill_buf()?;
        if buf.len() >= 2 {
            magic.copy_from_slice(&buf[0..2]);
        }
    }

    let q_counts = Arc::new(Mutex::new(HashMap::<u8, u64>::new()));
    let mut line_num = 0u64;
    let mut chunk = Vec::with_capacity(CHUNK_SIZE);

    // Process in chunks to limit memory usage
    if magic[0] == 0x1f && magic[1] == 0x8b {
        // Gzipped input
        let gz_reader = BufReader::new(MultiGzDecoder::new(reader));
        for line_result in gz_reader.lines() {
            let line = line_result?;
            if line_num % 4 == 3 {  // Quality line
                chunk.push(line);
                if chunk.len() >= CHUNK_SIZE / 4 {
                    process_chunk(&chunk, &q_counts);
                    chunk.clear();
                }
            }
            line_num += 1;
        }
    } else {
        // Plain text input
        for line_result in reader.lines() {
            let line = line_result?;
            if line_num % 4 == 3 {  // Quality line
                chunk.push(line);
                if chunk.len() >= CHUNK_SIZE / 4 {
                    process_chunk(&chunk, &q_counts);
                    chunk.clear();
                }
            }
            line_num += 1;
        }
    }

    // Process remaining chunk
    if !chunk.is_empty() {
        process_chunk(&chunk, &q_counts);
    }

    // Sort by Q score and output
    let q_counts = Arc::try_unwrap(q_counts).unwrap().into_inner().unwrap();
    let mut scores: Vec<_> = q_counts.iter().collect();
    scores.sort_by_key(|&(k, _)| k);

    let stdout = io::stdout();
    let mut handle = stdout.lock();

    for (q, count) in scores {
        writeln!(handle, "{}\t{}", q, count)?;
    }

    Ok(())
}

fn process_chunk(chunk: &[String], global_counts: &Arc<Mutex<HashMap<u8, u64>>>) {
    let local_counts: HashMap<u8, u64> = chunk
        .par_iter()
        .fold(
            || HashMap::new(),
            |mut acc, line| {
                for byte in line.bytes() {
                    let q_score = byte.saturating_sub(33);
                    *acc.entry(q_score).or_insert(0) += 1;
                }
                acc
            },
        )
        .reduce(
            || HashMap::new(),
            |mut acc, map| {
                for (q, count) in map {
                    *acc.entry(q).or_insert(0) += count;
                }
                acc
            },
        );

    // Merge into global counts
    let mut global = global_counts.lock().unwrap();
    for (q, count) in local_counts {
        *global.entry(q).or_insert(0) += count;
    }
}

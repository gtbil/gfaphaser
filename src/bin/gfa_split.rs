// gfa_split: split a GFA into per-component files.
//
// Components whose unique sequence length is >= --min-len get their own
// file; smaller components are concatenated into one small_components.gfa.
// A-lines are dropped. See README / source comments for details.

use gfaphaser::{cigar_overlap_len, DSU};
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process;

fn die(msg: &str) -> ! {
    eprintln!("error: {}", msg);
    process::exit(1);
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!(
            "usage: {} <input.gfa> [--out-prefix PREFIX] [--min-len BP]",
            args.get(0).map(|s| s.as_str()).unwrap_or("gfa_split")
        );
        process::exit(2);
    }

    let input_path = PathBuf::from(&args[1]);
    let mut out_prefix: Option<String> = None;
    let mut min_len: u64 = 100_000;

    let mut i = 2;
    while i < args.len() {
        match args[i].as_str() {
            "--out-prefix" => {
                i += 1;
                if i >= args.len() {
                    die("--out-prefix needs a value");
                }
                out_prefix = Some(args[i].clone());
            }
            "--min-len" => {
                i += 1;
                if i >= args.len() {
                    die("--min-len needs a value");
                }
                min_len = args[i]
                    .parse()
                    .unwrap_or_else(|_| die("--min-len must be a non-negative integer"));
            }
            other => die(&format!("unknown argument: {}", other)),
        }
        i += 1;
    }

    let prefix: String = out_prefix.unwrap_or_else(|| {
        Path::new(&input_path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("out")
            .to_string()
    });
    let out_dir = input_path
        .parent()
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| PathBuf::from("."));

    eprintln!(
        "pass 1: scanning segments and links from {}",
        input_path.display()
    );

    let f = File::open(&input_path).unwrap_or_else(|e| die(&format!("cannot open input: {}", e)));
    let reader = BufReader::new(f);

    let mut name_to_id: HashMap<String, usize> = HashMap::new();
    let mut seg_lengths: Vec<u64> = Vec::new();
    let mut link_overlaps: Vec<(usize, usize, u64)> = Vec::new();

    let intern = |name: &str,
                  name_to_id: &mut HashMap<String, usize>,
                  seg_lengths: &mut Vec<u64>|
     -> usize {
        if let Some(&id) = name_to_id.get(name) {
            return id;
        }
        let id = seg_lengths.len();
        name_to_id.insert(name.to_string(), id);
        seg_lengths.push(0);
        id
    };

    for line in reader.lines() {
        let line = line.unwrap_or_else(|e| die(&format!("read error: {}", e)));
        if line.is_empty() {
            continue;
        }
        let mut fields = line.split('\t');
        let rec = fields.next().unwrap_or("");
        match rec {
            "S" => {
                let name = fields.next().unwrap_or("");
                let seq = fields.next().unwrap_or("");
                if name.is_empty() {
                    continue;
                }
                let id = intern(name, &mut name_to_id, &mut seg_lengths);
                let len: u64 = if seq != "*" {
                    seq.len() as u64
                } else {
                    let mut ln: u64 = 0;
                    for tag in fields {
                        if let Some(rest) = tag.strip_prefix("LN:i:") {
                            ln = rest.parse().unwrap_or(0);
                            break;
                        }
                    }
                    ln
                };
                seg_lengths[id] = len;
            }
            "L" => {
                let from = fields.next().unwrap_or("");
                let _from_orient = fields.next();
                let to = fields.next().unwrap_or("");
                let _to_orient = fields.next();
                let overlap = fields.next().unwrap_or("*");
                if from.is_empty() || to.is_empty() {
                    continue;
                }
                let a = intern(from, &mut name_to_id, &mut seg_lengths);
                let b = intern(to, &mut name_to_id, &mut seg_lengths);
                let ov = cigar_overlap_len(overlap);
                link_overlaps.push((a, b, ov));
            }
            _ => {}
        }
    }

    let n_nodes = seg_lengths.len();
    eprintln!("  {} segments, {} links", n_nodes, link_overlaps.len());

    let mut dsu = DSU::new(n_nodes);
    for &(a, b, _) in &link_overlaps {
        dsu.union(a, b);
    }

    let mut comp_seg_len: HashMap<usize, u64> = HashMap::new();
    let mut comp_overlap: HashMap<usize, u64> = HashMap::new();
    for id in 0..n_nodes {
        let r = dsu.find(id);
        *comp_seg_len.entry(r).or_insert(0) += seg_lengths[id];
    }
    for &(a, _b, ov) in &link_overlaps {
        let r = dsu.find(a);
        *comp_overlap.entry(r).or_insert(0) += ov;
    }

    let mut comps: Vec<(usize, u64)> = comp_seg_len
        .iter()
        .map(|(&root, &slen)| {
            let ov = *comp_overlap.get(&root).unwrap_or(&0);
            (root, slen.saturating_sub(ov))
        })
        .collect();
    comps.sort_by(|a, b| b.1.cmp(&a.1));

    let mut root_to_target: HashMap<usize, Option<usize>> = HashMap::new();
    let mut large_idx: usize = 0;
    let mut n_small: usize = 0;
    let mut small_uniq_total: u64 = 0;
    let mut large_uniq_total: u64 = 0;
    for (root, uniq) in &comps {
        if *uniq >= min_len {
            large_idx += 1;
            root_to_target.insert(*root, Some(large_idx));
            large_uniq_total += *uniq;
        } else {
            root_to_target.insert(*root, None);
            n_small += 1;
            small_uniq_total += *uniq;
        }
    }
    eprintln!(
        "  {} components total: {} large (>= {} bp, total {} bp unique), {} small (total {} bp unique)",
        comps.len(),
        large_idx,
        min_len,
        large_uniq_total,
        n_small,
        small_uniq_total
    );

    let pad = std::cmp::max(2, large_idx.to_string().len());
    let mut large_writers: Vec<BufWriter<File>> = Vec::with_capacity(large_idx);
    for k in 1..=large_idx {
        let path = out_dir.join(format!(
            "{}.component_{:0width$}.gfa",
            prefix,
            k,
            width = pad
        ));
        let f = File::create(&path)
            .unwrap_or_else(|e| die(&format!("cannot create {}: {}", path.display(), e)));
        large_writers.push(BufWriter::new(f));
    }
    let small_path = out_dir.join(format!("{}.small_components.gfa", prefix));
    let small_file = File::create(&small_path)
        .unwrap_or_else(|e| die(&format!("cannot create {}: {}", small_path.display(), e)));
    let mut small_writer = BufWriter::new(small_file);

    let mut wrote_any_small = false;

    eprintln!("pass 2: writing records");

    let f2 = File::open(&input_path)
        .unwrap_or_else(|e| die(&format!("cannot reopen input: {}", e)));
    let reader2 = BufReader::new(f2);

    let route = |seg_name: &str,
                 name_to_id: &HashMap<String, usize>,
                 dsu: &mut DSU,
                 root_to_target: &HashMap<usize, Option<usize>>|
     -> Option<Option<usize>> {
        let id = *name_to_id.get(seg_name)?;
        let root = dsu.find(id);
        root_to_target.get(&root).copied()
    };

    for line in reader2.lines() {
        let line = line.unwrap_or_else(|e| die(&format!("read error: {}", e)));
        if line.is_empty() {
            continue;
        }
        let rec_end = line.find('\t').unwrap_or(line.len());
        let rec = &line[..rec_end];

        if rec == "A" {
            continue;
        }

        match rec {
            "H" => {
                for w in &mut large_writers {
                    writeln!(w, "{}", line).unwrap();
                }
                writeln!(&mut small_writer, "{}", line).unwrap();
            }
            "S" => {
                let mut it = line.split('\t');
                it.next();
                let name = it.next().unwrap_or("");
                if let Some(target) = route(name, &name_to_id, &mut dsu, &root_to_target) {
                    match target {
                        Some(k) => writeln!(&mut large_writers[k - 1], "{}", line).unwrap(),
                        None => {
                            writeln!(&mut small_writer, "{}", line).unwrap();
                            wrote_any_small = true;
                        }
                    }
                }
            }
            "L" | "C" => {
                let mut it = line.split('\t');
                it.next();
                let name = it.next().unwrap_or("");
                if let Some(target) = route(name, &name_to_id, &mut dsu, &root_to_target) {
                    match target {
                        Some(k) => writeln!(&mut large_writers[k - 1], "{}", line).unwrap(),
                        None => {
                            writeln!(&mut small_writer, "{}", line).unwrap();
                            wrote_any_small = true;
                        }
                    }
                }
            }
            "P" => {
                let mut it = line.split('\t');
                it.next();
                let _pname = it.next();
                let seg_csv = it.next().unwrap_or("");
                let first = seg_csv.split(',').next().unwrap_or("");
                let first = first.trim_end_matches(|c: char| c == '+' || c == '-');
                if let Some(target) = route(first, &name_to_id, &mut dsu, &root_to_target) {
                    match target {
                        Some(k) => writeln!(&mut large_writers[k - 1], "{}", line).unwrap(),
                        None => {
                            writeln!(&mut small_writer, "{}", line).unwrap();
                            wrote_any_small = true;
                        }
                    }
                }
            }
            "W" => {
                let mut it = line.split('\t');
                it.next();
                let _ = it.next();
                let _ = it.next();
                let _ = it.next();
                let _ = it.next();
                let _ = it.next();
                let walk = it.next().unwrap_or("");
                let first = walk
                    .split(|c| c == '>' || c == '<')
                    .find(|s| !s.is_empty())
                    .unwrap_or("");
                if !first.is_empty() {
                    if let Some(target) = route(first, &name_to_id, &mut dsu, &root_to_target) {
                        match target {
                            Some(k) => writeln!(&mut large_writers[k - 1], "{}", line).unwrap(),
                            None => {
                                writeln!(&mut small_writer, "{}", line).unwrap();
                                wrote_any_small = true;
                            }
                        }
                    }
                }
            }
            _ => {
                for w in &mut large_writers {
                    writeln!(w, "{}", line).unwrap();
                }
                writeln!(&mut small_writer, "{}", line).unwrap();
            }
        }
    }

    for mut w in large_writers {
        w.flush().unwrap();
    }
    small_writer.flush().unwrap();
    drop(small_writer);

    if !wrote_any_small {
        let _ = std::fs::remove_file(&small_path);
        eprintln!("  no small components; not writing {}", small_path.display());
    } else {
        eprintln!("  wrote small components -> {}", small_path.display());
    }

    eprintln!(
        "done. wrote {} large component file(s) to {}",
        large_idx,
        out_dir.display()
    );
}

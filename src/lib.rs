//! Shared library for GFA tools.
//!
//! Provides:
//! - GFA parsing (segments with lengths/sequences, links with overlap CIGARs)
//! - A bidirected graph representation suitable for walking
//! - Union-find for connected-component decomposition
//! - Common types for walks and step orientations
//!
//! Conventions:
//! - "Side" of a node: 0 = the `-` end (left), 1 = the `+` end (right).
//!   A walk *enters* a node through one side and *exits* through the other.
//!   In GFA orientation terms, a node visited as "+" means entering from
//!   side 0 (left) and exiting from side 1 (right); a node visited as "-"
//!   means entering from side 1 and exiting from side 0.
//! - An L-line `L A oa B ob ovl` says: the side of A indicated by `oa`
//!   connects to the side of B indicated by `ob`. Specifically, the
//!   exit-side of A under orientation oa equals oa-as-side (+ -> 1, - -> 0),
//!   and B is entered through the entry-side under orientation ob
//!   (+ -> 0, - -> 1). We store the link symmetrically so that a walk can
//!   traverse it in either direction.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

// ---------------- Union-Find ----------------

pub struct DSU {
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl DSU {
    pub fn new(n: usize) -> Self {
        DSU {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }
    pub fn find(&mut self, x: usize) -> usize {
        let mut r = x;
        while self.parent[r] != r {
            r = self.parent[r];
        }
        let mut cur = x;
        while self.parent[cur] != r {
            let nxt = self.parent[cur];
            self.parent[cur] = r;
            cur = nxt;
        }
        r
    }
    pub fn union(&mut self, a: usize, b: usize) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra == rb {
            return;
        }
        let (smaller, larger) = if self.rank[ra] < self.rank[rb] {
            (ra, rb)
        } else {
            (rb, ra)
        };
        self.parent[smaller] = larger;
        if self.rank[ra] == self.rank[rb] {
            self.rank[larger] += 1;
        }
    }
}

// ---------------- CIGAR overlap ----------------

pub fn cigar_overlap_len(cigar: &str) -> u64 {
    if cigar == "*" || cigar.is_empty() {
        return 0;
    }
    let mut total: u64 = 0;
    let mut num: u64 = 0;
    let mut saw_digit = false;
    for c in cigar.bytes() {
        if c.is_ascii_digit() {
            num = num * 10 + (c - b'0') as u64;
            saw_digit = true;
        } else {
            if saw_digit {
                match c {
                    b'M' | b'=' | b'X' | b'I' | b'D' => total += num,
                    _ => {}
                }
            }
            num = 0;
            saw_digit = false;
        }
    }
    total
}

// ---------------- Orientation helpers ----------------

/// Side 0 = '-' end (left), 1 = '+' end (right).
pub type Side = u8;

/// For a segment visited with the given orientation char, return (entry_side, exit_side).
pub fn entry_exit_for_orient(orient: char) -> (Side, Side) {
    match orient {
        '+' => (0, 1), // enter from left, exit from right
        '-' => (1, 0),
        _ => (0, 1),
    }
}

/// Encode an oriented node (node_id, orientation char) as an exit half-edge: 2*id + exit_side.
#[inline]
pub fn halfedge(node: usize, side: Side) -> usize {
    node * 2 + side as usize
}

// ---------------- Bidirected Graph ----------------

#[derive(Debug, Clone)]
pub struct Segment {
    pub name: String,
    pub length: u64,
    /// Inline sequence if present; None if "*" with no LN tag, or if we chose not to keep it.
    pub seq: Option<Vec<u8>>,
}

#[derive(Debug, Clone)]
pub struct Edge {
    /// Neighbor node id.
    pub to: usize,
    /// Side of the neighbor we enter.
    pub to_side: Side,
    /// Overlap length (bp) consumed when traversing this edge.
    pub overlap: u64,
    /// Original CIGAR string from the L-line (e.g. "55M", or "*").
    /// Preserved verbatim for round-tripping into P-line Overlaps fields.
    pub cigar: String,
}

/// Bidirected graph.
/// `adj[halfedge(node, side)]` = list of edges you can take by exiting `node` via `side`.
pub struct Graph {
    pub segs: Vec<Segment>,
    pub name_to_id: HashMap<String, usize>,
    /// 2 * n_segments entries. Each entry is a list of outgoing edges from that (node, exit-side).
    pub adj: Vec<Vec<Edge>>,
}

impl Graph {
    pub fn n(&self) -> usize {
        self.segs.len()
    }

    /// Edges available when exiting `node` from `side`.
    pub fn neighbors(&self, node: usize, exit_side: Side) -> &[Edge] {
        &self.adj[halfedge(node, exit_side)]
    }

    /// A node is a "tip" on `side` if it has zero outgoing edges from that side.
    /// Walks starting at a tip can begin by exiting the *other* side; walks
    /// ending at a tip terminate having exited into the tip side.
    pub fn is_tip_on(&self, node: usize, side: Side) -> bool {
        self.adj[halfedge(node, side)].is_empty()
    }
}

// ---------------- GFA parsing ----------------

#[derive(Default)]
pub struct ParseOptions {
    /// Keep inline sequences (memory!). Off by default; gfa_split doesn't need them.
    pub keep_sequences: bool,
}

pub struct ParsedGfa {
    pub graph: Graph,
    /// Raw H-lines preserved in order.
    pub headers: Vec<String>,
}

pub fn parse_gfa<P: AsRef<Path>>(path: P, opts: &ParseOptions) -> std::io::Result<ParsedGfa> {
    let f = File::open(path.as_ref())?;
    let reader = BufReader::new(f);

    let mut name_to_id: HashMap<String, usize> = HashMap::new();
    let mut segs: Vec<Segment> = Vec::new();
    let mut headers: Vec<String> = Vec::new();
    // Collect raw links first; build adjacency after we've seen all segments.
    let mut raw_links: Vec<(usize, char, usize, char, u64, String)> = Vec::new();

    let intern = |name: &str,
                      name_to_id: &mut HashMap<String, usize>,
                      segs: &mut Vec<Segment>|
     -> usize {
        if let Some(&id) = name_to_id.get(name) {
            return id;
        }
        let id = segs.len();
        name_to_id.insert(name.to_string(), id);
        segs.push(Segment {
            name: name.to_string(),
            length: 0,
            seq: None,
        });
        id
    };

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let mut fields = line.split('\t');
        let rec = match fields.next() {
            Some(r) => r,
            None => continue,
        };
        match rec {
            "H" => headers.push(line.clone()),
            "S" => {
                let name = fields.next().unwrap_or("");
                let seq = fields.next().unwrap_or("");
                if name.is_empty() {
                    continue;
                }
                let id = intern(name, &mut name_to_id, &mut segs);
                let (len, kept) = if seq != "*" {
                    let l = seq.len() as u64;
                    let kept = if opts.keep_sequences {
                        Some(seq.as_bytes().to_vec())
                    } else {
                        None
                    };
                    (l, kept)
                } else {
                    let mut ln: u64 = 0;
                    for tag in fields.by_ref() {
                        if let Some(rest) = tag.strip_prefix("LN:i:") {
                            ln = rest.parse().unwrap_or(0);
                            break;
                        }
                    }
                    (ln, None)
                };
                segs[id].length = len;
                segs[id].seq = kept;
            }
            "L" => {
                let from = fields.next().unwrap_or("");
                let fo = fields.next().unwrap_or("+").chars().next().unwrap_or('+');
                let to = fields.next().unwrap_or("");
                let to_o = fields.next().unwrap_or("+").chars().next().unwrap_or('+');
                let overlap = fields.next().unwrap_or("*");
                if from.is_empty() || to.is_empty() {
                    continue;
                }
                let a = intern(from, &mut name_to_id, &mut segs);
                let b = intern(to, &mut name_to_id, &mut segs);
                let ov = cigar_overlap_len(overlap);
                raw_links.push((a, fo, b, to_o, ov, overlap.to_string()));
            }
            _ => {}
        }
    }

    let n = segs.len();
    let mut adj: Vec<Vec<Edge>> = vec![Vec::new(); n * 2];

    for (a, fo, b, to_o, ov, cigar) in &raw_links {
        let a_exit: Side = if *fo == '+' { 1 } else { 0 };
        let b_entry: Side = if *to_o == '+' { 0 } else { 1 };

        adj[halfedge(*a, a_exit)].push(Edge {
            to: *b,
            to_side: b_entry,
            overlap: *ov,
            cigar: cigar.clone(),
        });
        adj[halfedge(*b, b_entry)].push(Edge {
            to: *a,
            to_side: a_exit,
            overlap: *ov,
            cigar: cigar.clone(),
        });
    }

    Ok(ParsedGfa {
        graph: Graph {
            segs,
            name_to_id,
            adj,
        },
        headers,
    })
}

// ---------------- Walks ----------------

/// A step in a walk: a node visited with an orientation ('+' or '-').
/// The first step has no incoming edge; each subsequent step has an
/// overlap (bp) and CIGAR string with the previous step.
#[derive(Debug, Clone)]
pub struct Step {
    pub node: usize,
    pub orient: char, // '+' or '-'
    /// Overlap with the previous step in bp. 0 for the first step.
    pub overlap_in: u64,
    /// Original CIGAR string for the link from the previous step to this
    /// step (e.g. "55M"). Empty string for the first step.
    pub cigar_in: String,
}

#[derive(Debug, Clone)]
pub struct Walk {
    pub steps: Vec<Step>,
}

impl Walk {
    pub fn new() -> Self {
        Walk { steps: Vec::new() }
    }

    /// Total walked sequence length: sum of segment lengths minus overlaps.
    pub fn length_bp(&self, graph: &Graph) -> u64 {
        let mut total: u64 = 0;
        for s in &self.steps {
            total = total.saturating_add(graph.segs[s.node].length);
            total = total.saturating_sub(s.overlap_in);
        }
        total
    }

    /// Set of distinct nodes touched.
    pub fn node_set(&self) -> std::collections::HashSet<usize> {
        self.steps.iter().map(|s| s.node).collect()
    }

    /// Number of revisits (visits beyond the first to a given node).
    pub fn revisits(&self) -> usize {
        let mut seen = std::collections::HashSet::new();
        let mut r = 0;
        for s in &self.steps {
            if !seen.insert(s.node) {
                r += 1;
            }
        }
        r
    }

    /// Render in GFA W-line walk format: ">u1<u2>u3"
    pub fn to_walk_string(&self, graph: &Graph) -> String {
        let mut s = String::new();
        for st in &self.steps {
            s.push(if st.orient == '+' { '>' } else { '<' });
            s.push_str(&graph.segs[st.node].name);
        }
        s
    }

    /// Render as the SegmentNames field of a GFA P-line:
    /// "u1+,u2+,u3-". A path with one step renders as "u1+".
    pub fn to_path_segment_names(&self, graph: &Graph) -> String {
        let parts: Vec<String> = self
            .steps
            .iter()
            .map(|st| format!("{}{}", graph.segs[st.node].name, st.orient))
            .collect();
        parts.join(",")
    }

    /// Render as the Overlaps field of a GFA P-line: comma-separated CIGAR
    /// strings, one for each consecutive pair of steps. For a walk with N
    /// steps, returns N-1 CIGARs (so "55M,55M" for a 3-step walk). For a
    /// 1-step walk returns "*" (no overlaps field needed but the spec
    /// requires the column be present).
    pub fn to_path_overlaps(&self) -> String {
        if self.steps.len() < 2 {
            return "*".to_string();
        }
        let parts: Vec<String> = self
            .steps
            .iter()
            .skip(1)
            .map(|st| {
                if st.cigar_in.is_empty() {
                    "*".to_string()
                } else {
                    st.cigar_in.clone()
                }
            })
            .collect();
        parts.join(",")
    }
}

impl Default for Walk {
    fn default() -> Self {
        Walk::new()
    }
}

// ---------------- Connected components (using DSU on graph) ----------------

/// Compute connected components of the underlying undirected graph.
/// Returns (component_id_per_node, n_components, members[ci]).
pub fn connected_components(g: &Graph) -> (Vec<usize>, usize, Vec<Vec<usize>>) {
    let n = g.n();
    let mut dsu = DSU::new(n);
    for u in 0..n {
        for side in [0u8, 1u8] {
            for e in g.neighbors(u, side) {
                dsu.union(u, e.to);
            }
        }
    }
    let mut root_to_ci: HashMap<usize, usize> = HashMap::new();
    let mut comp_of = vec![0usize; n];
    let mut members: Vec<Vec<usize>> = Vec::new();
    for u in 0..n {
        let r = dsu.find(u);
        let ci = match root_to_ci.get(&r) {
            Some(&ci) => ci,
            None => {
                let ci = members.len();
                root_to_ci.insert(r, ci);
                members.push(Vec::new());
                ci
            }
        };
        comp_of[u] = ci;
        members[ci].push(u);
    }
    (comp_of, members.len(), members)
}

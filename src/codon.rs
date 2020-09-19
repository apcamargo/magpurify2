// This file is part of the magpurify2 package, available at:
// https://github.com/apcamargo/magpurify2
//
// Magpurify2 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.
//
// Contact: antoniop.camargo@gmail.com

use lazy_static::lazy_static;
use ndarray::{Array1, Array2};
use numpy::convert::ToPyArray;
use pyo3::{prelude::*, wrap_pyfunction};
use std::collections::{HashMap, HashSet};
use std::iter;

lazy_static! {
    static ref SYNONYMOUS_CODONS: HashMap<&'static str, Vec<&'static str>> = [
        ("CYS", vec!["TGT", "TGC"]),
        ("ASP", vec!["GAT", "GAC"]),
        ("SER", vec!["TCT", "TCG", "TCA", "TCC", "AGC", "AGT"]),
        ("GLN", vec!["CAA", "CAG"]),
        ("MET", vec!["ATG"]),
        ("ASN", vec!["AAC", "AAT"]),
        ("PRO", vec!["CCT", "CCG", "CCA", "CCC"]),
        ("LYS", vec!["AAG", "AAA"]),
        ("THR", vec!["ACC", "ACA", "ACG", "ACT"]),
        ("PHE", vec!["TTT", "TTC"]),
        ("ALA", vec!["GCA", "GCC", "GCG", "GCT"]),
        ("GLY", vec!["GGT", "GGG", "GGA", "GGC"]),
        ("ILE", vec!["ATC", "ATA", "ATT"]),
        ("LEU", vec!["TTA", "TTG", "CTC", "CTT", "CTG", "CTA"]),
        ("HIS", vec!["CAT", "CAC"]),
        ("ARG", vec!["CGA", "CGC", "CGG", "CGT", "AGG", "AGA"]),
        ("TRP", vec!["TGG"]),
        ("VAL", vec!["GTA", "GTC", "GTG", "GTT"]),
        ("GLU", vec!["GAG", "GAA"]),
        ("TYR", vec!["TAT", "TAC"]),
        ("STOP", vec!["TAG", "TGA", "TAA"])
    ]
    .iter()
    .cloned()
    .collect();
}

/// Counts the absolute codon frequencies from a list of in-frame ORF sequences. If the sequence
/// length if not a multiple of 3, the last nucleotides will be ignored.
fn count_codons(seq_list: Vec<&str>, use_pseudocount: bool) -> HashMap<&str, u32> {
    let pseudo: u32;
    if use_pseudocount {
        pseudo = 1;
    } else {
        pseudo = 0;
    }
    let mut codon_count: HashMap<&str, u32> = [
        ("TTT", pseudo),
        ("TTC", pseudo),
        ("TTA", pseudo),
        ("TTG", pseudo),
        ("CTT", pseudo),
        ("CTC", pseudo),
        ("CTA", pseudo),
        ("CTG", pseudo),
        ("ATT", pseudo),
        ("ATC", pseudo),
        ("ATA", pseudo),
        ("ATG", pseudo),
        ("GTT", pseudo),
        ("GTC", pseudo),
        ("GTA", pseudo),
        ("GTG", pseudo),
        ("TAT", pseudo),
        ("TAC", pseudo),
        ("TAA", pseudo),
        ("TAG", pseudo),
        ("CAT", pseudo),
        ("CAC", pseudo),
        ("CAA", pseudo),
        ("CAG", pseudo),
        ("AAT", pseudo),
        ("AAC", pseudo),
        ("AAA", pseudo),
        ("AAG", pseudo),
        ("GAT", pseudo),
        ("GAC", pseudo),
        ("GAA", pseudo),
        ("GAG", pseudo),
        ("TCT", pseudo),
        ("TCC", pseudo),
        ("TCA", pseudo),
        ("TCG", pseudo),
        ("CCT", pseudo),
        ("CCC", pseudo),
        ("CCA", pseudo),
        ("CCG", pseudo),
        ("ACT", pseudo),
        ("ACC", pseudo),
        ("ACA", pseudo),
        ("ACG", pseudo),
        ("GCT", pseudo),
        ("GCC", pseudo),
        ("GCA", pseudo),
        ("GCG", pseudo),
        ("TGT", pseudo),
        ("TGC", pseudo),
        ("TGA", pseudo),
        ("TGG", pseudo),
        ("CGT", pseudo),
        ("CGC", pseudo),
        ("CGA", pseudo),
        ("CGG", pseudo),
        ("AGT", pseudo),
        ("AGC", pseudo),
        ("AGA", pseudo),
        ("AGG", pseudo),
        ("GGT", pseudo),
        ("GGC", pseudo),
        ("GGA", pseudo),
        ("GGG", pseudo),
    ]
    .iter()
    .cloned()
    .collect();
    for seq in seq_list {
        let n_codons = seq.len() / 3;
        for i in 0..n_codons {
            let codon = &seq[i * 3..i * 3 + 3];
            if codon_count.contains_key(codon) {
                let count = codon_count.entry(codon).or_insert(0);
                *count += 1;
            }
        }
    }
    codon_count
}

/// Computed the Codon Adaptation Index (CAI) of a single in-frame ORF sequence. A codon usage index
/// is required as it stores the weight of each codon within a group of synonymous codons.
fn cai_for_gene(seq: &str, index: &HashMap<&str, f64>) -> f64 {
    let mut cai_value = 0.0;
    let mut cai_length = 0.0;
    let n_codons = seq.len() / 3;
    for i in 0..n_codons {
        let codon = &seq[i * 3..i * 3 + 3];
        if index.contains_key(codon) {
            if !vec!["ATG", "TGG"].contains(&codon) {
                cai_value += index.get(codon).unwrap().ln();
                cai_length += 1.0;
            }
        }
    }
    return (cai_value / (cai_length - 1.0)).exp();
}

/// get_codon_index(seq_list, use_pseudocount=True)
/// --
///
/// Computes the codon usage index (RSCU) from a list of in-frame ORF sequences.
///
/// Parameters
/// ----------
/// seq_list : list
///    In-frame DNA sequences of the ORFs that will be used to compute the index.
/// use_pseudocount : bool, optional
///    Use a pseudocount of 1 when counting the codons. Default is True.
///
/// Returns
/// -------
/// dict
///    A dictionary where the keys are codons and the values correspond to the
///    the weight of the given codon within its group of synonymous codons.
#[pyfunction(use_pseudocount = "true")]
fn get_codon_index(py: Python, seq_list: Vec<&str>, use_pseudocount: bool) -> PyObject {
    let mut index = HashMap::new();
    let codon_count = count_codons(seq_list, use_pseudocount);
    for aa in SYNONYMOUS_CODONS.keys() {
        let mut total = 0.0;
        let mut rcsu = vec![];
        let codons = SYNONYMOUS_CODONS.get(aa).unwrap();
        for codon in codons {
            total += *codon_count.get(codon).unwrap() as f64;
        }
        let denominator = total / codons.len() as f64;
        if total == 0.0 {
            rcsu = iter::repeat(0.0).take(codons.len()).collect();
        } else {
            for codon in codons {
                rcsu.push(*codon_count.get(codon).unwrap() as f64 / denominator);
            }
        }
        let rcsu_max = rcsu.iter().cloned().fold(0. / 0., f64::max);
        for (codon_index, codon) in codons.iter().enumerate() {
            if rcsu[codon_index] == 0.0 {
                index.insert(*codon, 0.0);
            } else {
                index.insert(*codon, rcsu[codon_index] / rcsu_max);
            }
        }
    }
    index.into_py(py)
}

/// get_cai(seq_list, use_pseudocount=True)
/// --
///
/// Computes Codon Adaptation Index (CAI) values for a list of in-frame ORF
/// sequences.
///
/// Parameters
/// ----------
/// seq_list : list
///    In-frame DNA sequences of the ORFs whose CAI will be computed.
/// codon_index : dict
///    A codon usage index generated with the `get_codon_index` function.
///
/// Returns
/// -------
/// ndarray
///    An array containing the CAI of each of the sequences in the input list.
#[pyfunction]
fn get_cai(py: Python, seq_list: Vec<&str>, codon_index: HashMap<&str, f64>) -> PyObject {
    let cai_array: Array1<f64> = seq_list
        .iter()
        .map(|seq| cai_for_gene(seq, &codon_index))
        .collect();
    cai_array
        .to_pyarray(Python::acquire_gil().python())
        .into_py(py)
}

/// get_codon_usage_profile(seq_list)
/// --
///
/// Computes five features that summarize the codon usage profile of an ORF
/// sequence (Hughes, A. L., & Langley, K. J. (2007)).
///
/// Parameters
/// ----------
/// seq_list : list
///    In-frame DNA sequences of the ORFs that will be used to compute the index.
///
/// Returns
/// -------
/// ndarray
///    An array containing the five descriptor features (pA4, pC4, pG4, pC2, and
///    pG2) for each input ORF.
#[pyfunction]
fn get_codon_usage_profile(py: Python, seq_list: Vec<&str>) -> PyObject {
    let n4_set: HashSet<&str> = vec![
        "GCT", "GCC", "GCA", "GCG", "CGT", "CGC", "CGA", "CGG", "GGT", "GGC", "GGA", "GGG", "CTT",
        "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "TCT", "TCC", "TCA", "TCG", "ACT", "ACC",
        "ACA", "ACG", "GTT", "GTC", "GTA", "GTG",
    ]
    .into_iter()
    .collect();
    let ag2_set: HashSet<&str> = vec![
        "AGA", "AGG", "CAA", "CAG", "GAA", "GAG", "TTA", "TTG", "AAA", "AAG", "AGT", "AGC", "TAT",
        "TAC",
    ]
    .into_iter()
    .collect();
    let tc2_set: HashSet<&str> = vec![
        "AAT", "AAC", "GAT", "GAC", "TGT", "TGC", "CAT", "CAC", "TTT", "TTC",
    ]
    .into_iter()
    .collect();

    let mut feature_vector = vec![];
    for seq in &seq_list {
        let mut n4 = 0;
        let mut ag2 = 0;
        let mut tc2 = 0;
        let mut pa4 = 0;
        let mut pc4 = 0;
        let mut pg4 = 0;
        let mut pc2 = 0;
        let mut pg2 = 0;
        let n_codons = seq.len() / 3;
        for i in 0..n_codons {
            let codon = &seq[i * 3..i * 3 + 3];
            let last_base = codon.chars().last().unwrap();
            if n4_set.contains(codon) {
                n4 += 1;
                match last_base {
                    'A' => pa4 += 1,
                    'C' => pc4 += 1,
                    'G' => pg4 += 1,
                    _ => (),
                }
            } else if ag2_set.contains(codon) {
                ag2 += 1;
                match last_base {
                    'G' => pg2 += 1,
                    _ => (),
                }
            } else if tc2_set.contains(codon) {
                tc2 += 1;
                match last_base {
                    'C' => pc2 += 1,
                    _ => (),
                }
            }
        }
        let pa4: f64 = (pa4 as f64 + 1.0) / (n4 as f64 + 1.0);
        let pc4: f64 = (pc4 as f64 + 1.0) / (n4 as f64 + 1.0);
        let pg4: f64 = (pg4 as f64 + 1.0) / (n4 as f64 + 1.0);
        let pc2: f64 = (pc2 as f64 + 1.0) / (tc2 as f64 + 1.0);
        let pg2: f64 = (pg2 as f64 + 1.0) / (ag2 as f64 + 1.0);
        feature_vector.extend(vec![pa4, pc4, pg4, pc2, pg2]);
    }
    let feature_vector = Array2::from_shape_vec((seq_list.len(), 5), feature_vector).unwrap();
    feature_vector
        .to_pyarray(Python::acquire_gil().python())
        .into_py(py)
}

#[pymodule]
fn _codon(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(get_codon_index))?;
    m.add_wrapped(wrap_pyfunction!(get_cai))?;
    m.add_wrapped(wrap_pyfunction!(get_codon_usage_profile))?;
    Ok(())
}

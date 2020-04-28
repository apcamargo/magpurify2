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

use ndarray::{Array2, Axis};
use numpy::convert::ToPyArray;
use pyo3::{prelude::*, wrap_pyfunction};
use rayon::{prelude::*, ThreadPoolBuilder};

static REV_COMPLEMENTARY_FOURMERS: [usize; 256] = [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 30, 11, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 23, 42, 43, 44, 7, 45, 46,
    47, 48, 49, 50, 51, 34, 52, 53, 54, 19, 55, 56, 57, 3, 58, 59, 60, 57, 61, 62, 63, 44, 64, 65,
    66, 30, 67, 68, 69, 14, 70, 71, 72, 54, 73, 74, 75, 41, 76, 77, 78, 26, 79, 80, 66, 10, 81, 82,
    83, 51, 84, 85, 86, 37, 87, 88, 75, 22, 89, 90, 63, 6, 91, 92, 93, 47, 94, 95, 83, 33, 96, 97,
    72, 18, 98, 99, 60, 2, 100, 101, 99, 56, 102, 103, 90, 43, 104, 105, 80, 29, 106, 107, 68, 13,
    108, 109, 97, 53, 110, 111, 88, 40, 112, 113, 77, 25, 114, 105, 65, 9, 115, 116, 95, 50, 117,
    118, 85, 36, 119, 111, 74, 21, 120, 103, 62, 5, 121, 122, 92, 46, 123, 116, 82, 32, 124, 109,
    71, 17, 125, 101, 59, 1, 126, 125, 98, 55, 127, 120, 89, 42, 128, 114, 79, 28, 129, 106, 67,
    12, 130, 124, 96, 52, 131, 119, 87, 39, 132, 112, 76, 24, 128, 104, 64, 8, 133, 123, 94, 49,
    134, 117, 84, 35, 131, 110, 73, 20, 127, 102, 61, 4, 135, 121, 91, 45, 133, 115, 81, 31, 130,
    108, 70, 16, 126, 100, 58, 0,
];

/// Counts the absolute frequencies of the k-mers of size `k` contained in a nucleotide sequence and
/// returns them in a vector. Only k-mers containing 'A', 'C', 'G', 'T', 'a', 'c', 'g' or 't' are
/// counted.
fn count_kmers(seq: &str, k: u32) -> std::vec::Vec<u32> {
    let mut counts = vec![0; 4_usize.pow(k)];
    let mut kmer = 0;
    let mut countdown = k - 1;
    let byteseq = seq.as_bytes();
    let mask = (1 << (2 * k)) - 1;
    for base in byteseq {
        match base {
            65 | 97 => kmer = ((kmer << 2) | 0) & mask,
            67 | 99 => kmer = ((kmer << 2) | 1) & mask,
            71 | 103 => kmer = ((kmer << 2) | 2) & mask,
            84 | 116 => kmer = ((kmer << 2) | 3) & mask,
            _ => countdown = k,
        }
        if countdown == 0 {
            counts[kmer] += 1;
        } else {
            countdown -= 1;
        }
    }
    counts
}

/// Counts the absolute frequencies of the canonical 4-mers in a nucleotide sequence and returns
/// them in a vector. Only k-mers containing 'A', 'C', 'G', 'T', 'a', 'c', 'g' or 't' are counted.
fn count_canonical_fourmers(seq: &str) -> std::vec::Vec<u32> {
    let mut canonical_fourmer_counts = vec![0; 136];
    let counts = count_kmers(seq, 4);
    for i in 0..256 {
        canonical_fourmer_counts[REV_COMPLEMENTARY_FOURMERS[i]] += counts[i];
    }
    canonical_fourmer_counts
}

/// get_tnf(seq_list, relative=True, threads=1)
/// --
///
/// Computes the 4-mer frequencies of a list of DNA sequences.
///
/// Parameters
/// ----------
/// seq_list : list
///    DNA sequence string from which k-mers will be counted.
/// relative : bool, optional
///    Return relative 4-mer frequencies. Default is True.
/// threads : int, optional
///    Number of threads to use for TNF computation. Default is 1.
///
/// Returns
/// -------
/// ndarray
///    An array containing the TNFs of each sequence in the input.
#[pyfunction(relative = "true", threads = "1")]
fn get_tnf(py: Python, seq_list: Vec<&str>, relative: bool, threads: usize) -> PyObject {
    let pool = ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap();
    let counts_array = pool.install(|| {
        Array2::from_shape_vec(
            (seq_list.len(), 136),
            seq_list
                .par_iter()
                .map(|seq| count_canonical_fourmers(seq))
                .flatten()
                .collect(),
        )
        .unwrap()
    });
    if relative {
        let counts_array = counts_array.mapv(|elem| elem as f32);
        let sums = counts_array.sum_axis(Axis(1));
        (counts_array.reversed_axes() / sums)
            .reversed_axes()
            .mapv(|elem| if elem.is_nan() { 0. } else { elem })
            .to_pyarray(Python::acquire_gil().python())
            .into_py(py)
    } else {
        counts_array
            .to_pyarray(Python::acquire_gil().python())
            .into_py(py)
    }
}

#[pymodule]
fn _tnf(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(get_tnf))?;
    Ok(())
}

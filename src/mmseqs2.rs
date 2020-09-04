use csv;
use pyo3::{prelude::*, wrap_pyfunction};
use std::collections::HashMap;
use std::error::Error;

/// Parses the content of the MMSeqs2 taxonomy output file and returns a dictionary where the keys
/// are the query ids and the values are their assigned taxids.
fn parse_mmseqs2_taxonomy(
    mmseqs2_taxonomy_output: &str,
) -> Result<HashMap<String, u32>, Box<dyn Error>> {
    let mut gene_taxonomy_dict = HashMap::new();
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(mmseqs2_taxonomy_output)?;
    let mut query: String;
    let mut taxid: u32;
    for result in rdr.records() {
        let record = result?;
        query = record[0].to_string();
        taxid = record[1].parse()?;
        gene_taxonomy_dict.insert(query, taxid);
    }
    Ok(gene_taxonomy_dict)
}

/// Parses the content of the MMSeqs2 alignment output file and returns two dictionaries. In the
/// first, the keys are the query ids and the values are 'alignment identity * bitscore'. In the
/// seconda, the keys are the query ids and the values are the alignment bitscores.
fn parse_mmseqs2_alignment(
    mmseqs2_alignment_output: &str,
) -> Result<(HashMap<String, f32>, HashMap<String, u32>), Box<dyn Error>> {
    let mut gene_identity_dict = HashMap::new();
    let mut gene_bitscore_dict = HashMap::new();
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_path(mmseqs2_alignment_output)?;
    let mut query: String;
    let mut fident: f32;
    let mut tcov: f32;
    let mut bitscore: u32;
    for result in rdr.records() {
        let record = result?;
        query = record[0].to_string();
        fident = record[1].parse()?;
        tcov = record[2].parse()?;
        bitscore = record[3].parse()?;
        gene_identity_dict.insert(query.clone(), fident * tcov);
        gene_bitscore_dict.insert(query.clone(), bitscore);
    }
    Ok((gene_identity_dict, gene_bitscore_dict))
}

/// get_mmseqs2(mmseqs2_taxonomy_output, mmseqs2_alignment_output)
/// --
///
/// Parses the taxonomy and alignment of a MMSeqs2 execution and returns a
/// nested dictionary containing the assigned taxonomy and alignment details for
/// each gene that had a match in the database search.
///
/// Parameters
/// ----------
/// mmseqs2_taxonomy_output : str
///    Path to the MMSeqs2 taxonomy output file.
/// mmseqs2_alignment_output : str
///    Path to the MMSeqs2 alignment output file.
///
/// Returns
/// -------
/// dictionary
///    A nested dictionary where the keys of the outer dictionary are the genome
///    names and the values are inner dictionaries where the keys are contig
///    names and the values are tuples in the following format:
///    (taxid, identity * target coverage, bitscore).
#[pyfunction]
fn get_mmseqs2(
    py: Python,
    mmseqs2_taxonomy_output: &str,
    mmseqs2_alignment_output: &str,
) -> PyObject {
    let gene_taxonomy_dict = parse_mmseqs2_taxonomy(mmseqs2_taxonomy_output).unwrap();
    let (gene_identity_dict, gene_bitscore_dict) =
        parse_mmseqs2_alignment(mmseqs2_alignment_output).unwrap();
    let mut mmseqs2_dict: HashMap<&str, HashMap<&str, (Vec<u32>, Vec<f32>, Vec<u32>)>> =
        HashMap::new();
    for (query, _) in &gene_identity_dict {
        let taxid = *gene_taxonomy_dict.get(query).unwrap();
        let identity = *gene_identity_dict.get(query).unwrap();
        let bitscore = *gene_bitscore_dict.get(query).unwrap();
        let v: Vec<&str> = query.rsplitn(3, "~").collect();
        if let [_, contig, genome] = &v[..] {
            let genome_dictionary = mmseqs2_dict.entry(genome).or_insert(HashMap::new());
            let contig_tuple = genome_dictionary
                .entry(contig)
                .or_insert((vec![], vec![], vec![]));
            contig_tuple.0.push(taxid);
            contig_tuple.1.push(identity);
            contig_tuple.2.push(bitscore);
        }
    }
    mmseqs2_dict.into_py(py)
}

#[pymodule]
fn _mmseqs2(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(get_mmseqs2))?;
    Ok(())
}

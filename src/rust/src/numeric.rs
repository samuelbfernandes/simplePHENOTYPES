use extendr_api::prelude::*;

// R represents NA_integer_ as i32::MIN.
const R_NA_INT: i32 = i32::MIN;

/// Numericalize a raw 0/1/2 dosage matrix into the user-facing coding.
///
/// `raw_dosage` is a flat integer vector of length `n_snp * n_samp` stored
/// **row-major** (i.e. all samples for SNP 0 come first, then all samples for
/// SNP 1, …).  Values must be 0 (hom allele-1), 1 (het), 2 (hom allele-2),
/// or `NA_integer_` (missing).
///
/// `flip[i]` = TRUE when allele-2 is the **major** allele for SNP i (so raw 2
/// should receive the major-allele code).  FALSE means allele-1 is major.
///
/// Returns an integer vector of the same length with values recoded per
/// `code_as` / `model` / `impute`.
///
/// @param raw_dosage Integer vector length n_snp * n_samp, row-major.
/// @param n_snp      Number of SNPs (rows).
/// @param n_samp     Number of samples (columns).
/// @param flip       Logical vector length n_snp.
/// @param code_as    "-101" (major=1, het=0, minor=-1) or "012" (major=2, het=1, minor=0).
/// @param model      "Add", "Dom", "Left", or "Right".
/// @param impute     "None", "Middle", "Minor", or "Major".
/// @return Integer vector length n_snp * n_samp.
/// @noRd
#[extendr]
pub fn numericalize_core(
    raw_dosage: &[i32],
    n_snp: i32,
    n_samp: i32,
    flip: Logicals,
    code_as: &str,
    model: &str,
    impute: &str,
) -> Vec<i32> {
    let n_snp = n_snp as usize;
    let n_samp = n_samp as usize;

    let (major_val, het_val, minor_val): (i32, i32, i32) = if code_as == "012" {
        (2, 1, 0)
    } else {
        (1, 0, -1) // default "-101"
    };

    let impute_val: Option<i32> = match impute {
        "Middle" => Some(het_val),
        "Minor"  => Some(minor_val),
        "Major"  => Some(major_val),
        _        => None, // "None"
    };

    let flip_vec: Vec<bool> = flip.into_iter().map(|v| v.is_true()).collect();
    let mut out = vec![R_NA_INT; n_snp * n_samp];

    for snp in 0..n_snp {
        let is_flipped = flip_vec.get(snp).copied().unwrap_or(false);
        let base = snp * n_samp;

        for samp in 0..n_samp {
            let raw = raw_dosage[base + samp];

            out[base + samp] = if raw == R_NA_INT {
                impute_val.unwrap_or(R_NA_INT)
            } else {
                // Map raw dosage → Add-model code
                let add_coded = match (raw, is_flipped) {
                    (0, false) => major_val,
                    (1, _)     => het_val,
                    (2, false) => minor_val,
                    (0, true)  => minor_val,
                    (2, true)  => major_val,
                    _          => R_NA_INT,
                };

                if add_coded == R_NA_INT {
                    impute_val.unwrap_or(R_NA_INT)
                } else {
                    // Apply genetic-model post-processing on top of Add coding.
                    // These match the transforms in the original numericalization():
                    //   Dom:   non-het → minor_val  (x1[x1 != 0] <- -1)
                    //   Left:  het     → minor_val  (x1[x1 == 0] <- -1)
                    //   Right: het     → major_val  (x1[x1 == 0] <- 1)
                    match model {
                        "Dom"   => if add_coded != het_val { minor_val } else { het_val },
                        "Left"  => if add_coded == het_val { minor_val } else { add_coded },
                        "Right" => if add_coded == het_val { major_val } else { add_coded },
                        _       => add_coded, // "Add" and unrecognised
                    }
                }
            };
        }
    }

    out
}

extendr_module! {
    mod numeric;
    fn numericalize_core;
}

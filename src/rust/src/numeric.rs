use extendr_api::prelude::*;

// R represents NA_integer_ as i32::MIN.
const R_NA_INT: i32 = i32::MIN;

/// Coding scheme: the values given to the major homozygote, the heterozygote
/// and the minor homozygote.
fn coding(code_as: &str) -> Result<(i32, i32, i32), String> {
    match code_as {
        "-101" => Ok((1, 0, -1)),
        "012" => Ok((2, 1, 0)),
        other => Err(format!(
            "code_as must be \"-101\" or \"012\", got {:?}",
            other
        )),
    }
}

/// Genetic-model post-transform applied on top of the Add coding.
#[derive(Clone, Copy)]
enum Model {
    Add,
    Dom,
    Left,
    Right,
}

fn parse_model(model: &str) -> Result<Model, String> {
    match model {
        "Add" => Ok(Model::Add),
        "Dom" => Ok(Model::Dom),
        "Left" => Ok(Model::Left),
        "Right" => Ok(Model::Right),
        other => Err(format!(
            "model must be \"Add\", \"Dom\", \"Left\" or \"Right\", got {:?}",
            other
        )),
    }
}

/// Convert a non-negative, non-NA R integer to a size.
fn size_arg(v: i32, what: &str) -> Result<usize, String> {
    usize::try_from(v).map_err(|_| {
        format!(
            "{} must be a non-negative whole number, got {}",
            what,
            if v == R_NA_INT {
                "NA".to_string()
            } else {
                v.to_string()
            }
        )
    })
}

/// The recoding itself, on validated inputs.
///
/// `raw` is row-major (`raw[snp * n_samp + samp]`); every value must be 0, 1,
/// 2 or `R_NA_INT`, and `flip` holds one flag per SNP. See [`numericalize_core`]
/// for the meaning of each argument. Nothing is defaulted: an unknown option,
/// a wrongly sized input or an out-of-domain call is an error, because each of
/// those used to produce a plausible but wrong coding silently (a short or NA
/// `flip` was read as FALSE, a call outside {0, 1, 2, NA} became missing, an
/// over-long `raw` was truncated, an unknown `code_as` / `model` / `impute` fell
/// back to `-101` / `Add` / `None`).
pub fn numericalize(
    raw: &[i32],
    n_snp: usize,
    n_samp: usize,
    flip: &[bool],
    code_as: &str,
    model: &str,
    impute: &str,
) -> Result<Vec<i32>, String> {
    let (major_val, het_val, minor_val) = coding(code_as)?;
    let model = parse_model(model)?;
    let impute_val: Option<i32> = match impute {
        "None" => None,
        "Middle" => Some(het_val),
        "Minor" => Some(minor_val),
        "Major" => Some(major_val),
        other => {
            return Err(format!(
                "impute must be \"None\", \"Middle\", \"Minor\" or \"Major\", got {:?}",
                other
            ))
        }
    };
    let n_cells = n_snp
        .checked_mul(n_samp)
        .ok_or_else(|| "n_snp * n_samp is too large".to_string())?;
    if raw.len() != n_cells {
        return Err(format!(
            "raw_dosage has {} values but n_snp * n_samp = {} * {} = {}",
            raw.len(),
            n_snp,
            n_samp,
            n_cells
        ));
    }
    if flip.len() != n_snp {
        return Err(format!(
            "flip has {} entries but there are {} SNPs",
            flip.len(),
            n_snp
        ));
    }
    if let Some(k) = raw
        .iter()
        .position(|&v| v != R_NA_INT && !(0..=2).contains(&v))
    {
        return Err(format!(
            "raw dosage values must be 0, 1, 2 or NA; found {} at position {}",
            raw[k],
            k + 1
        ));
    }

    let mut out = vec![R_NA_INT; n_cells];

    for (snp, &is_flipped) in flip.iter().enumerate() {
        let base = snp * n_samp;

        for samp in 0..n_samp {
            let raw_v = raw[base + samp];

            // The Add-model code comes from the observed dosage, or from the
            // imputation class when the call is missing (Major -> major_val,
            // Minor -> minor_val, Middle -> het_val); a missing call with
            // impute = "None" stays NA. The genetic-model transform below is then
            // applied to imputed calls too, so an imputed major homozygote is
            // coded exactly like an observed one under Dom/Left/Right. (Returning
            // the raw impute value here previously skipped the transform, so an
            // imputed genotype came out with the wrong Dom/Left/Right code.)
            let add_coded = if raw_v == R_NA_INT {
                impute_val.unwrap_or(R_NA_INT)
            } else {
                match (raw_v, is_flipped) {
                    (0, false) => major_val,
                    (1, _) => het_val,
                    (2, false) => minor_val,
                    (0, true) => minor_val,
                    (2, true) => major_val,
                    // unreachable: the domain was validated above
                    _ => R_NA_INT,
                }
            };

            out[base + samp] = if add_coded == R_NA_INT {
                R_NA_INT
            } else {
                // Apply genetic-model post-processing on top of Add coding.
                //   Dom:   non-het → minor_val  (x1[x1 != 0] <- -1)
                //   Left:  het     → minor_val  (x1[x1 == 0] <- -1)
                //   Right: het     → major_val  (x1[x1 == 0] <- 1)
                match model {
                    Model::Dom => {
                        if add_coded != het_val {
                            minor_val
                        } else {
                            het_val
                        }
                    }
                    Model::Left => {
                        if add_coded == het_val {
                            minor_val
                        } else {
                            add_coded
                        }
                    }
                    Model::Right => {
                        if add_coded == het_val {
                            major_val
                        } else {
                            add_coded
                        }
                    }
                    Model::Add => add_coded,
                }
            };
        }
    }

    Ok(out)
}

/// Numericalize a raw 0/1/2 dosage matrix into the user-facing coding.
///
/// `raw_dosage` is a flat integer vector of length `n_snp * n_samp` stored
/// **row-major** (i.e. all samples for SNP 0 come first, then all samples for
/// SNP 1, …).  Values must be 0 (hom allele-1), 1 (het), 2 (hom allele-2),
/// or `NA_integer_` (missing); anything else is an error.
///
/// `flip[i]` = TRUE when allele-2 is the **major** allele for SNP i (so raw 2
/// should receive the major-allele code).  FALSE means allele-1 is major.
/// `flip` must have exactly one non-missing entry per SNP.
///
/// Tie rule (MAF = 0.5): the R caller (`compute_flip()`) sets `flip` only when
/// raw 2 is *strictly* more frequent than raw 0, so a marker with equal counts
/// is not flipped and allele 1 (raw 0) is coded major. That is deterministic but
/// arbitrary: which allele is `allele 1` depends on the file format's own
/// convention. v1.3.0 differed for a marker with no heterozygote and exactly two
/// genotype classes (its `len == 2` branch coded the *major* class -1, the
/// opposite sign of every other marker, and re-coded missing calls); v2 codes
/// such markers like every other.
///
/// Returns an integer vector of the same length with values recoded per
/// `code_as` / `model` / `impute`. Returns `list(ok, err)` (feature
/// `result_list`); the R wrapper raises `err` as an R error.
///
/// @param raw_dosage Integer vector length n_snp * n_samp, row-major.
/// @param n_snp      Number of SNPs (rows).
/// @param n_samp     Number of samples (columns).
/// @param flip       Logical vector length n_snp, no NA.
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
) -> std::result::Result<Vec<i32>, String> {
    let n_snp = size_arg(n_snp, "n_snp")?;
    let n_samp = size_arg(n_samp, "n_samp")?;
    let mut flip_vec: Vec<bool> = Vec::with_capacity(flip.len());
    for (k, v) in flip.iter().enumerate() {
        if v.is_na() {
            return Err(format!("flip must not contain NA (entry {})", k + 1));
        }
        flip_vec.push(v.is_true());
    }
    numericalize(raw_dosage, n_snp, n_samp, &flip_vec, code_as, model, impute)
}

extendr_module! {
    mod numeric;
    fn numericalize_core;
}

#[cfg(test)]
mod tests {
    use super::*;

    const NA: i32 = R_NA_INT;

    fn run(
        raw: &[i32],
        n_snp: usize,
        n_samp: usize,
        flip: &[bool],
        c: &str,
        m: &str,
        i: &str,
    ) -> Vec<i32> {
        numericalize(raw, n_snp, n_samp, flip, c, m, i).unwrap()
    }

    #[test]
    fn additive_101_unflipped_and_flipped() {
        // one SNP, raw 0/1/2/NA
        let raw = [0, 1, 2, NA];
        assert_eq!(
            run(&raw, 1, 4, &[false], "-101", "Add", "None"),
            vec![1, 0, -1, NA]
        );
        assert_eq!(
            run(&raw, 1, 4, &[true], "-101", "Add", "None"),
            vec![-1, 0, 1, NA]
        );
        assert_eq!(
            run(&raw, 1, 4, &[false], "012", "Add", "None"),
            vec![2, 1, 0, NA]
        );
    }

    #[test]
    fn each_snp_uses_its_own_flip() {
        let raw = [0, 2, 0, 2];
        assert_eq!(
            run(&raw, 2, 2, &[false, true], "-101", "Add", "None"),
            vec![1, -1, -1, 1]
        );
    }

    #[test]
    fn full_valid_grid_matches_an_independent_mapper() {
        // 2 codings x 4 models x 4 imputations x 2 flips = 64 combinations,
        // every raw class (0, 1, 2, NA), against a tiny reference mapper.
        let raw = [0, 1, 2, NA];
        for code_as in ["-101", "012"] {
            let (maj, het, min) = if code_as == "012" {
                (2, 1, 0)
            } else {
                (1, 0, -1)
            };
            for model in ["Add", "Dom", "Left", "Right"] {
                for impute in ["None", "Middle", "Minor", "Major"] {
                    for flip in [false, true] {
                        let got = run(&raw, 1, 4, &[flip], code_as, model, impute);
                        let want: Vec<i32> = raw
                            .iter()
                            .map(|&r| {
                                let add = if r == NA {
                                    match impute {
                                        "Middle" => het,
                                        "Minor" => min,
                                        "Major" => maj,
                                        _ => NA,
                                    }
                                } else if r == 1 {
                                    het
                                } else if (r == 0) != flip {
                                    maj
                                } else {
                                    min
                                };
                                if add == NA {
                                    NA
                                } else {
                                    match model {
                                        "Dom" => {
                                            if add != het {
                                                min
                                            } else {
                                                het
                                            }
                                        }
                                        "Left" => {
                                            if add == het {
                                                min
                                            } else {
                                                add
                                            }
                                        }
                                        "Right" => {
                                            if add == het {
                                                maj
                                            } else {
                                                add
                                            }
                                        }
                                        _ => add,
                                    }
                                }
                            })
                            .collect();
                        assert_eq!(got, want, "{} {} {} flip={}", code_as, model, impute, flip);
                    }
                }
            }
        }
    }

    #[test]
    fn malformed_inputs_are_errors_not_defaults() {
        let raw = [0, 2];
        let n = |r: &[i32], s, m, f: &[bool], c: &str, mo: &str, i: &str| {
            numericalize(r, s, m, f, c, mo, i)
        };
        // short / long flip
        assert!(n(&raw, 2, 1, &[true], "-101", "Add", "None").is_err());
        assert!(n(&raw, 2, 1, &[true, false, true], "-101", "Add", "None").is_err());
        // raw length vs n_snp * n_samp (short and long)
        assert!(n(&raw, 2, 2, &[true, false], "-101", "Add", "None").is_err());
        assert!(n(&[0, 1, 2, 0, 1], 1, 2, &[false], "-101", "Add", "None").is_err());
        // out-of-domain call
        assert!(n(&[0, 3], 1, 2, &[false], "-101", "Add", "None").is_err());
        assert!(n(&[0, -1], 1, 2, &[false], "-101", "Add", "None").is_err());
        // unknown options
        assert!(n(&raw, 1, 2, &[false], "abc", "Add", "None").is_err());
        assert!(n(&raw, 1, 2, &[false], "-101", "add", "None").is_err());
        assert!(n(&raw, 1, 2, &[false], "-101", "Add", "none").is_err());
        // empty is fine
        assert_eq!(n(&[], 0, 0, &[], "-101", "Add", "None"), Ok(vec![]));
    }

    #[test]
    fn size_args_reject_na_and_negative() {
        assert_eq!(size_arg(3, "n"), Ok(3));
        assert!(size_arg(-1, "n").is_err());
        assert!(size_arg(NA, "n").is_err());
    }
}

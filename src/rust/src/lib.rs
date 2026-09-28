use extendr_api::prelude::*;

mod genome;
mod hash;
mod meiosis;
mod numeric;

extendr_module! {
    mod simplePHENOTYPES;
    use numeric;
    use meiosis;
    use hash;
}

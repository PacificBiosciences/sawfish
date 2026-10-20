use std::collections::HashSet;

#[derive(Debug)]
struct IgnoreMacros(HashSet<String>);

impl bindgen::callbacks::ParseCallbacks for IgnoreMacros {
    fn will_parse_macro(&self, name: &str) -> bindgen::callbacks::MacroParsingBehavior {
        if self.0.contains(name) {
            bindgen::callbacks::MacroParsingBehavior::Ignore
        } else {
            bindgen::callbacks::MacroParsingBehavior::Default
        }
    }
}

fn main() {
    let out_dir = cmake::Config::new("WFA2-lib")
        .cflag("-DCMAKE_BUILD_TYPE=Release")
        // As recommended by the README on master.
        .cflag("-DEXTRA_FLAGS=\"-ftree-vectorize -msse2 -mfpmath=sse -ftree-vectorizer-verbose=5\"")
        .build();
    println!("cargo:rustc-link-search=native={}/lib", out_dir.display());
    println!("cargo:rustc-link-search=native={}/lib64", out_dir.display());
    println!("cargo:rustc-link-lib=static=wfa2");

    let ignored_macros = IgnoreMacros(
        vec![
            "FP_INFINITE".into(),
            "FP_NAN".into(),
            "FP_NORMAL".into(),
            "FP_SUBNORMAL".into(),
            "FP_ZERO".into(),
            "IPPORT_RESERVED".into(),
        ]
        .into_iter()
        .collect(),
    );

    bindgen::Builder::default()
        .header("WFA2-lib/utils/commons.h")
        .header("WFA2-lib/wavefront/wfa.h")
        .clang_arg("--include-directory=WFA2-lib")
        .parse_callbacks(Box::new(ignored_macros))
        .generate()
        .expect("Unable to generate bindings")
        .write_to_file(out_dir.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}

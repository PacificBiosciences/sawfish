#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
#![allow(improper_ctypes)]
#![allow(deref_nullptr)]
// bindgen declares libc functions (memcmp, strlen, ...) with slightly different signatures
#![allow(unknown_lints)]
#![allow(suspicious_runtime_symbol_definitions)]

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

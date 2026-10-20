#[cxx::bridge(namespace = "spoa")]
pub mod ffi {
    /// Enumerates the possible alignment types
    #[repr(u32)]
    enum AlignmentType {
        /// Smith-Waterman alignment
        kSW,
        /// Needleman-Wunsch alignment
        kNW,
        /// Overlap alignment
        kOV,
    }

    unsafe extern "C++" {
        include!("spoa-sys/include/bindings.hpp");

        type Alignment;
        type AlignmentEngine;
        type AlignmentType;
        type Graph;

        fn create_graph() -> UniquePtr<Graph>;

        /// # Safety
        /// this function is unsafe because cxx says pointer arguments are unsafe
        unsafe fn clear(graph: Pin<&mut Graph>);

        /// # Safety
        /// this function is unsafe because cxx says pointer arguments are unsafe
        unsafe fn add_alignment(
            graph: Pin<&mut Graph>,
            alignment: &Alignment,
            sequence: *const c_char,
            sequence_len: u32,
            quality: *const c_char,
            quality_len: u32,
        );

        /// # Safety
        /// this function is unsafe because cxx says pointer arguments are unsafe
        unsafe fn add_alignment_noqual(
            graph: Pin<&mut Graph>,
            alignment: &Alignment,
            sequence: *const c_char,
            sequence_len: u32,
        );

        fn generate_consensus(graph: Pin<&mut Graph>) -> UniquePtr<CxxString>;
        fn generate_multiple_sequence_alignment(
            graph: Pin<&mut Graph>,
            include_consensus: bool,
        ) -> UniquePtr<CxxVector<CxxString>>;

        fn get_sequence_count(graph: &Graph) -> u32;

        fn create_alignment_engine(
            typ: AlignmentType,
            m: i8,
            n: i8,
            g: i8,
            e: i8,
            q: i8,
            c: i8,
        ) -> UniquePtr<AlignmentEngine>;

        fn create_alignment() -> UniquePtr<Alignment>;

        fn alignment_is_empty(alignment: &Alignment) -> bool;

        /// # Safety
        /// this function is unsafe because cxx says pointer arguments are unsafe
        unsafe fn align(
            alignment_engine: Pin<&mut AlignmentEngine>,
            sequence: *const c_char,
            sequence_len: u32,
            graph: &Graph,
            score: &mut i32,
        ) -> UniquePtr<Alignment>;

        /// # Safety
        /// this function is unsafe because cxx says pointer arguments are unsafe
        unsafe fn align_pruned(
            alignment_engine: Pin<&mut AlignmentEngine>,
            sequence: *const c_char,
            sequence_len: u32,
            graph: &Graph,
            min_score: i32,
            full_fallback: bool,
            score: &mut i32,
            reached: &mut bool,
        ) -> UniquePtr<Alignment>;

        /// # Safety
        /// this function is unsafe because cxx says pointer arguments are unsafe
        unsafe fn align_pruned_auto(
            alignment_engine: Pin<&mut AlignmentEngine>,
            sequence: *const c_char,
            sequence_len: u32,
            graph: &Graph,
            expected_fracs: *const f64,
            num_expected_fracs: u32,
            score: &mut i32,
        ) -> UniquePtr<Alignment>;

        fn max_possible_score(
            alignment_engine: Pin<&mut AlignmentEngine>,
            sequence_len: u32,
            graph: &Graph,
        ) -> i64;

        fn get_alignment_overlap_size(alignment: &Alignment) -> u32;

        fn get_alignment_clip_size(alignment: &Alignment) -> i32;
    }
}

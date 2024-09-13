/// Global sawfish version number
///
/// All client code should refer directly to this copy instead of using various possibly conflicting environment variables
pub const SAWFISH_VERSION: &str = env!("VERGEN_GIT_DESCRIBE");

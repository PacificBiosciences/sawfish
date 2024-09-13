pub use log::debug;

/// Print log message for either of two debug scenarios:
/// (1) The global debug log level has been activated
/// (2) A local debug flag has been enabled
///
/// The local debug flag is given as the first argument, and determines whether the debug message is directly printed to stdout.
///
/// # Examples
///
/// ```
/// debug_msg!(false, "At task foo {}", x); // prints debug log msg only if global user --debug flag is given
/// debug_msg!(true, "At task foo {}", x); // prints directly to stderr
/// ```
macro_rules! debug_msg {
    ($flag:expr, $($arg:tt)+) => {
        if $flag {
            eprintln!($($arg)+);
        } else {
            $crate::log_utils::debug!($($arg)+);
        }
    }
}

pub(crate) use debug_msg;

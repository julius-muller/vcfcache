//! Optional coarse spans. Disabled runs do not take per-record timestamps.
use std::{
    fs::OpenOptions,
    io::Write,
    path::PathBuf,
    sync::{Mutex, OnceLock},
    time::{Instant, SystemTime, UNIX_EPOCH},
};

static DIRECTORY: OnceLock<Option<PathBuf>> = OnceLock::new();
static LOCK: Mutex<()> = Mutex::new(());

pub struct Span {
    name: &'static str,
    start: Option<(Instant, u128)>,
    details: serde_json::Value,
}

impl Span {
    pub fn new(name: &'static str) -> Self {
        let enabled = DIRECTORY
            .get_or_init(|| std::env::var_os("VCFCACHE_PROFILE_DIR").map(PathBuf::from))
            .is_some();
        Self {
            name,
            start: enabled.then(|| {
                (
                    Instant::now(),
                    SystemTime::now()
                        .duration_since(UNIX_EPOCH)
                        .unwrap()
                        .as_nanos(),
                )
            }),
            details: serde_json::Value::Null,
        }
    }
    pub fn details(&mut self, value: serde_json::Value) {
        if self.start.is_some() {
            self.details = value;
        }
    }
}

impl Drop for Span {
    fn drop(&mut self) {
        let Some((start, epoch)) = self.start else {
            return;
        };
        let Some(directory) = DIRECTORY.get().and_then(Option::as_ref) else {
            return;
        };
        let record = serde_json::json!({
            "schema_version":1, "name":self.name, "pid":std::process::id(),
            "tid":format!("{:?}",std::thread::current().id()),
            "start_unix_ns":epoch, "duration_ns":start.elapsed().as_nanos(),
            "outcome":"scope_exit", "details":self.details
        });
        let _lock = LOCK.lock().unwrap_or_else(|error| error.into_inner());
        let result = (|| -> std::io::Result<()> {
            let mut file = OpenOptions::new()
                .append(true)
                .create(true)
                .open(directory.join(format!("rust-{}.jsonl", std::process::id())))?;
            writeln!(file, "{record}")
        })();
        if let Err(error) = result {
            eprintln!("PROFILE_WRITE_FAILED: {error}");
        }
    }
}

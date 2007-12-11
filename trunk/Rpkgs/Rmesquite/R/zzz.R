.First.lib <- function(lib,pkg) {
  require(rJava);
  .jpackage(pgk);
}

## Ideally we would want to clean up behind ourselves on exit, but at
## present there doesn't seem to be a good way to da this. Instead, we
## strongly recommend to initialize the JVM with forcing
## re-initialization, to startup with a clean environment.

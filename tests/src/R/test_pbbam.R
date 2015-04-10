
# usage: R [args] < test_pbbam.R --args <tests_scripts_dir> <PacBioBam_lib_dir> <test_data_dir>
args <- commandArgs(TRUE)
tests_path     <- args[1]
lib_path       <- args[2]
test_data_path <- args[3]

# load PacBioBAM lib & wrapper script
pbbam_libname <- paste(lib_path, "PacBioBam",   sep="/")
pbbam_wrapper <- paste(lib_path, "PacBioBam.R", sep="/")
dyn.load(paste(pbbam_libname, .Platform$dynlib.ext, sep=""))
source(pbbam_wrapper)
cacheMetaData(1)

# init test utils & run test cases
source(paste(tests_path, "utils.R", sep="/"))
run_test_suite(tests_path)

# print results & exit
results <- test_suite_results()
results$print_summary()
if (results$any_failed())
	quit(status=1)

# Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted (subject to the limitations in the
# disclaimer below) provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above
#    copyright notice, this list of conditions and the following
#    disclaimer in the documentation and/or other materials provided
#    with the distribution.
#
#  * Neither the name of Pacific Biosciences nor the names of its
#    contributors may be used to endorse or promote products derived
#    from this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
# GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
# BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
# USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
#
# Author: Derek Barnett

# main test suite runner
run_test_suite <- function(path) {
	test_files <- dir(path, "^test.*\\.[rR]$", full.names = TRUE)
	lapply(test_files, run_test_file)
	invisible()
}

run_test_file <- function(filename) {
	source(filename)
	invisible()
}

# main test case definition
#
# Example:
# test_case("name", {
#     ...tested code here...   
# })
#
test_case <- function(name, code) {
	test_case_runner(name, substitute(code))
	invisible()
}

# main assert definitions
#
# assertEqual(expected, actual)
# assertNotEqual(expected, actual)
# assertTrue(expr)
# assertFalse(expr)
#

assertEqual <- function(expected, actual) {
	assertHelper(identical(expected, actual), 
	             TRUE, 
				 expression_string("expected"), 
				 deparse(expected), 
				 expression_string("actual"), 
				 deparse(actual), 
				 "==")
}

assertAllEqual <- function(expected, actual) {
	assertHelper(all.equal(expected, actual), 
	             TRUE, 
				 expression_string("expected"), 
				 deparse(expected), 
				 expression_string("actual"), 
				 deparse(actual), 
				 "==")
}

assertNotEqual <- function(expected, actual) {
	assertHelper(identical(expected, actual), 
	             FALSE, 
				 expression_string("expected"), 
				 deparse(expected), 
				 expression_string("actual"), 
				 deparse(actual), 
				 "!=")
}

assertTrue <- function(expr) {
	assertHelper(as.vector(expr), 
	             TRUE, 
				 "TRUE", 
				 "TRUE", 
				 expression_string("expr"), 
				 deparse(expr), 
				 "==")
}

assertFalse <- function(expr) {
	assertHelper(as.vector(expr), 
	             FALSE, 
				 "FALSE", 
				 "FALSE", 
				 expression_string("expr"), 
				 deparse(expr), 
				 "==")
}

# TODO: (if needed) assertLessThan, assertGreaterThan, assertNull, etc

# ------------------------------- #
# internals 
# ------------------------------- #

expression_string <- function(name, env = parent.frame()) {
  subs <- do.call("substitute", list(as.name(name), env))
  paste0(deparse(subs, width.cutoff = 500), collapse = "\n")
}

assertHelper <- function(compare, 
	                     to, 
	                     expected_expr, 
						 expected_value, 
						 actual_expr, 
						 actual_value,
						 compare_type) 
{
	success <- identical(compare, to)
	
	result <- make_assert_result(success,
    				             expected_expr, 
                                 expected_value, 
 						         actual_expr, 
 						         actual_value,
 						         compare_type)
	
	
	# record result with testCaseCollector
	testCaseResults <- test_case_results()
	testCaseResults$add_result(result)
	invisible()
}

make_assert_result <- function(success,
	                           expected_expr, 
					           expected_value, 
					           actual_expr, 
					           actual_value,
					           compare_type)
{
	structure(list(
		success = success,
		expectedExpression = expected_expr,
		expectedValue = expected_value,
		actualExpression = actual_expr,
		actualValue = actual_value,
		compareType = compare_type	
	))					 	
}

TestCaseResults <- setRefClass("TestCaseResults",

	fields = list(
		test = "character",
		anyFailed = "logical",
		results = "list"
	),
	methods = list(
		initialize = function(...) {
			test <<- ""
			anyFailed <<- FALSE
		
			initFields(...)
		},
	
		start_test = function(name) {
			test <<- name
			results <<- list()
		},
		
		add_result = function(result) {
			if (!result$success) 
				anyFailed <<- TRUE
			results <<- c(results, list(result))
		},
		
		end_test = function() {
		
			# summarize test case results
			testOutput <- format_results()
		
			# report to test collector
			suiteResults <- test_suite_results()
			suiteResults$add_test_case_result(anyFailed, testOutput)
		
			# reset
			test <<- ""
			anyFailed <<- FALSE
			results <<- list()
		},
		
		format_results = function() {
			
			lines <- list()
			
			status <- "OK"
			if (anyFailed) 
				status <- "FAILED"
			
			headerLine <- paste("TestCase:", test, "...", status, sep=" ")
			lines <- c(lines, list(headerLine))
			
			for (result in results) {
				if (!result$success) {
					valueOfLabel <- paste(result$actualExpression, result$compareType, result$expectedExpression, sep=" ")
					valueOf  <- paste("  Value of:", valueOfLabel,         sep=" ")
					actual   <- paste("    Actual:", result$actualValue,   sep=" ")
					expected <- paste("  Expected:", result$expectedValue, sep=" ")
					lines <- c(lines, valueOf, actual, expected, "")
				}
			}
			invisible(lines)
		} 
	)
)

TestSuiteResults <- setRefClass("TestSuiteResults",

	fields = list(
		numTests = "integer",
		numFailed = "integer",
		results = "list"
	),
	methods = list(
		initialize = function(...) {
			numTests <<- 0L
			numFailed <<- 0L
			results <<- list()
		
			initFields(...)
		},
		
		add_test_case_result = function(testHasFailures, testOutput) {  #(results)
			numTests <<- numTests + 1L
			if (testHasFailures)
				numFailed <<- numFailed + 1L
			results <<- c(results, testOutput)
		},
		
		any_failed = function() {
			return (numFailed != 0L)
		},
		
		print_summary = function(...) {
				
			cat("\n")
			cat("-------------------------\n")
			cat("Tests Complete\n")
			cat("-------------------------\n")
			cat("\n")
			
			for (result in results) {
				cat(result)
				cat("\n")
			}
			cat("-------------------------\n")
			
			if (numFailed == 1L) {
				footer <- paste(numFailed, "test failed out of", numTests, sep=" ")
			} else {
				footer <- paste(numFailed, "tests failed out of", numTests, sep=" ")
			}
			cat(footer)
			cat("\n\n")
		}
	)
)

test_env = new.env()
test_env$testSuiteResults <- TestSuiteResults$new()
test_env$testCaseResults <- TestCaseResults$new()

test_suite_results <- function() {
	test_env$testSuiteResults
}

test_case_results <- function() {
	test_env$testCaseResults
}

test_case_runner <- function(name, code) {
	testCaseResults <- test_case_results()
	testCaseResults$start_test(name)
	eval(code, test_env)
	testCaseResults$end_test()
}

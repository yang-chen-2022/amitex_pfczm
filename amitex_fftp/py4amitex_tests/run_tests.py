# -*- coding:utf-8 -*-
from argparse import ArgumentParser
from contextlib import redirect_stdout
from importlib import import_module
import re
from os import putenv, getenv
from io import StringIO
import subprocess
from time import time
import traceback
from xml.sax.saxutils import escape as escapexml


import numpy as np

from py4amitex.simulation import getSimulationShellCommand
from py4amitex.output.amitexoutput import AmitexOutput

from testlist import TESTS


FAILCC = "\033[91m"
OKCC = "\033[92m"
ENDCC = "\033[0m"


def compare_mechanics(test, testresults):
    TOLERANCE = 1.0e-8
    # TODO: add option to silence AmitexOutput
    p4amout = StringIO()
    with redirect_stdout(p4amout):
        results = AmitexOutput(testresults)
        results.load_mean_values()
        ref = AmitexOutput(f"../resultats/{test}/reference/{test}")
        ref.load_mean_values()

    output = ""
    ok = True
    for quantity in ["stress", "strain", "niter"]:
        qref = ref.get_mean_values(quantity)
        qres = results.get_mean_values(quantity)
        norm = np.sqrt(np.sum(qref**2))
        error = np.sqrt(np.sum(qref - qres) ** 2)
        ok = error <= norm * TOLERANCE
        output += f"{quantity}: error rel = {error/norm:.7e} abs = {error:.7e}\n"

    return (ok, output, "")


def compare_diffusion(test, testresults):
    comp = subprocess.run(
        f"../validation/compare ../resultats/{test}/reference/{test} {testresults} -n {test} -s tests.log",
        shell=True,
        capture_output=True,
    )
    stdout = comp.stdout.decode()
    stderr = comp.stderr.decode()
    return (comp.returncode == 0 and (stdout.find("ERROR") == -1), stdout, stderr)


def compare(test, testresults, kind):
    if kind == "mechanics":
        return compare_mechanics(test, testresults)
    elif kind == "diffusion":
        return compare_diffusion(test, testresults)
    else:
        raise ValueError("kind for comparison is 'mechanics' or 'diffusion'")


def run_simulation(input, numberProcs):
    input.generateFiles()
    cmd = getSimulationShellCommand(input, numberProcs=numberProcs)
    comp = subprocess.run(cmd, shell=True, capture_output=True)
    return (comp.returncode == 0, comp.stdout.decode(), comp.stderr.decode())


def run_tests(regexp=""):
    testresults = []
    for testdesc in TESTS:
        test = testdesc["name"]
        typ = testdesc.get("type", "mechanics")
        numberProcs = testdesc.get("nprocs", 0)
        if (not regexp) or re.match(regexp, test):
            mod = import_module(f"tests.{test}")
            print(f"{test}... ", end="", flush=True)
            startTime = time()
            try:
                input = mod.make_input()
            except Exception as exc:
                duration = time() - startTime
                print(
                    f"{FAILCC}FAIL (in py4amitex){ENDCC} (executed in  {duration:.8g}s)"
                )
                testresults.append(
                    (
                        test,
                        False,
                        "",
                        "",
                        str(exc),
                        type(exc).__name__,
                        traceback.format_exc(),
                        duration,
                    )
                )
                continue

            input.resultsDir = f"resultats/{test}"
            ok, stdout, stderr = run_simulation(input, numberProcs)
            if ok:
                ok2, stdout2, stderr2 = compare(test, input.outputPrefix(), typ)
                ok = ok and ok2
                failureBody = stdout2
                stderr += stderr2
                message = "Comparaison failed"
                errorKind = "ValidationError"
            else:
                message = "Simulation failed"
                errorKind = "SimulationError"
                stdoutLines = stdout.splitlines()
                failureBody = "(last lines of stdout)...\n" + "\n".join(
                    stdoutLines[max(0, len(stdoutLines) - 10) : len(stdoutLines)]
                )
            duration = time() - startTime
            testresults.append(
                (test, ok, stdout, stderr, message, errorKind, failureBody, duration)
            )
            CC = OKCC if ok else FAILCC

            print(
                f"{CC}{'PASS' if ok else 'FAIL'}{ENDCC} (executed in  {duration:.6g}s)"
            )
    return testresults


def print_test_results(testresults, showOutput):
    nfail = 0
    for test, ok, stdout, stderr, _, _, failureBody, _ in testresults:
        if not ok:
            nfail += 1
            print(f"{FAILCC}Failed {test} {ENDCC}")
            print(failureBody)
            if showOutput:
                print(f"\nstdout{'-' * 40}\n{stdout}{'-' * 46}")
                if stderr:
                    print(f"\nstderr{'-' * 40}\n{stderr}{'-' * 46}")
            print()
    if nfail > 0:
        print(f"{FAILCC}{nfail} test(s) failed{ENDCC} out of {len(testresults)}")
    else:
        print(f"{OKCC}{len(testresults)} tests passed{ENDCC}")


def test_results_to_junit(testresults):

    totalDuration = sum(t[-1] for t in testresults)
    totalFails = sum(int(not t[1]) for t in testresults)
    result = f"""<?xml version="1.0" encoding="UTF-8"?>
<testsuites time="{totalDuration:.8g}" tests="{len(testresults)}" failures="{totalFails}" errors="0" skipped="0">
"""

    result += f"""<testsuite name="py4amitex" time="{totalDuration:.8g}" tests="{len(testresults)}" failures="{totalFails}" errors="0" skipped="0">
"""

    for (
        test,
        ok,
        stdout,
        stderr,
        message,
        errorKind,
        failureBody,
        duration,
    ) in testresults:
        result += f"""<testcase name="{test}" classname="py4amitex" time="{duration:.8g}" >
"""

        result += f"<system-out>{escapexml(stdout)}</system-out>\n"
        if stderr:
            result += f"<system-err>{escapexml(stderr)}</system-err>\n"
        if not ok:
            result += f"""<failure message="{escapexml(message)}" type="{errorKind}">{escapexml(failureBody)}</failure>
"""
        result += "</testcase>\n"

    result += "</testsuite>"
    result += "</testsuites>\n"
    return result


def set_environment():
    amitexDir = getenv("AMITEX_PATH", default="..")
    path = getenv("PATH")
    putenv("PATH", f"{amitexDir}/libAmitex/bin:{path}")
    ldpath = getenv("LD_LIBRARY_PATH")
    putenv("LD_LIBRARY_PATH", f"{amitexDir}/libAmitex/lib:{ldpath}")
    putenv("AMITEX_PATH", amitexDir)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-T", "--tests-regex", default="", help="Run tests matching regular expression"
    )
    parser.add_argument(
        "--junit-output",
        default="",
        help="Export test to result to XML file JUNIT_OUTPUT",
    )
    parser.add_argument(
        "--no-show-output",
        action='store_true',
        default=False,
        help="Do not show standard output/error of failed tests",
    )
    args = parser.parse_args()
    set_environment()
    testresults = run_tests(args.tests_regex)
    print()
    print_test_results(testresults, not args.no_show_output)
    if args.junit_output:
        with open(args.junit_output, "w", encoding="utf-8") as file:
            file.write(test_results_to_junit(testresults))

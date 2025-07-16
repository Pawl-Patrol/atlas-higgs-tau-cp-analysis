#!/usr/bin/env python3

import inquirer
import subprocess

# - selection of samples using multi-select
# - selection of debug (on/off)
# - selection of Events with Default 100
SAMPLES = {
    "cp-even\thadhad": "cp-even-hadhad",
    "cp-odd\thadhad": "cp-odd-hadhad",
    "cp-even\thadlep": "cp-even-hadlep",
    "cp-odd\thadlep": "cp-odd-hadlep",
    "cp-even\tlephad": "cp-even-lephad",
    "cp-odd\tlephad": "cp-odd-lephad",
    "cp-even\thadhad\tH200": "cp-even-hadhad-H200",
    "cp-even\thadhad\tH500": "cp-even-hadhad-H500",
    "cp-even\thadhad\tH1000": "cp-even-hadhad-H1000",
    "cp-even\thadhad\tH3000": "cp-even-hadhad-H3000",
}

questions = [
    inquirer.Checkbox(
        "samples",
        message="Select the samples you want to use",
        choices=list(SAMPLES.keys()),
    ),
    inquirer.Text(
        "events",
        message="How many events do you want to process?",
        default="100",
        validate=lambda _, x: x.isdigit() or x == "-1",
    ),
    inquirer.Confirm(
        "debug",
        message="Do you want to run in debug mode?",
        default=False,
    ),
]

if __name__ == "__main__":
    answers = inquirer.prompt(questions)
    events = int(answers["events"])
    debug = answers["debug"]

    subprocess.run(["make"], cwd="/srv/build")

    for sample in answers["samples"]:
        out = SAMPLES[sample]
        dir = "/samples/" + out
        cmd = ["ATestRun_eljob.py", "-c", dir, "-s", out, "-e", answers["events"]]
        if debug:
            cmd.append("--debug")
        subprocess.run(cmd, cwd="/srv/run")

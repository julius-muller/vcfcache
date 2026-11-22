"""Run a cache's annotation recipe on an input and nothing else.

This is the control of the cohort benchmark: the same recipe, substituted and
executed exactly as vcfcache does it, without cache lookup or output checks.
"""

import argparse
from pathlib import Path

from vcfcache.database.preserving import run_recipe
from vcfcache.database.workflow_manager import BcftoolsCommand, WorkflowManager


def main():
    """Substitute the recipe of one cache and run it once."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--params", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True)
    manager = WorkflowManager(
        input_file=args.input,
        output_dir=args.work,
        name="plain",
        anno_config_file=args.cache / "annotation.yaml",
        params_file=args.params,
        output_file=args.output,
        verbosity=1,
    )
    manager._create_work_dir(args.work, dirname="work")
    work = manager.work_dir

    def command(text):
        return BcftoolsCommand(text, manager.logger, work).run()

    run_recipe(manager, work, command, args.input, args.output)


if __name__ == "__main__":
    main()

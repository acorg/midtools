#!/bin/bash

# A Git pre-commit hook.
#
# To install:
#
#   $ cd .git/hooks
#   $ ln -s ../../bin/pre-commit.sh pre-commit

# Add our virtualenv bin to PATH. Git commit seems to muck with PATH :-(
if [ -n "$VIRTUAL_ENV" ]
then
    PATH="$VIRTUAL_ENV/bin:$PATH"
fi

ruff check

if [ $? -ne 0 ]
then
    echo "COMMIT FAILED: 'ruff check' did not run cleanly:" >&2
    exit 1
fi

make check

if [ $? -ne 0 ]
then
    echo "COMMIT FAILED: 'make check' did not run cleanly:" >&2
    exit 1
fi
